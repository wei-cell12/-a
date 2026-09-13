def volumetric_heat_capacity(cdec) :
    rho = 650.0 + 128.0 * cdec
    cp = 1450.0 + 2736.0 * cdec / (cdec + 1.0)
    return rho * cp


def conductivity(cdec) :
    return 0.21 + 0.38 * cdec / (cdec + 1.0)


def diffusivity(t_c, cdec) :
    if np.any(cdec <= 0.0):
        raise FloatingPointError("含水率必须为正，程序不使用截断掩盖数值错误")
    return 2.4e-3 * np.exp(-0.45 / cdec) * np.exp(-3850.0 / (t_c + 273.15))


def radial_geometry(cfg: Config) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if cfg.intervals < 4 or cfg.mesh_power < 1.0:
        raise ValueError("intervals至少为4，mesh_power不得小于1")
    xi = np.linspace(0.0, 1.0, cfg.intervals + 1)
    r = RADIUS_M * (1.0 - (1.0 - xi) ** cfg.mesh_power)
    faces = np.concatenate(([0.0], 0.5 * (r[:-1] + r[1:]), [RADIUS_M]))
    volumes = 0.5 * np.diff(faces**2)
    return r, faces, volumes


def spatial_operator(coefficient, boundary_coefficient: float,
                     boundary_value: float, geom: tuple[np.ndarray, ...]
                     ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    r, faces, volumes = geom
    face_values = (2.0 * coefficient[:-1] * coefficient[1:] /
                   np.maximum(coefficient[:-1] + coefficient[1:], 1.0e-300))
    conce = faces[1:-1] * face_values / np.diff(r)
    lwr = np.zeros_like(r)
    upr = np.zeros_like(r)
    lwr[1:] = conce / volumes[1:]
    upr[:-1] = conce / volumes[:-1]
    diaal = -(lwr + upr)
    diaal[-1] -= RADIUS_M * boundary_coefficient / volumes[-1]
    soce = np.zeros_like(r)
    soce[-1] = RADIUS_M * boundary_coefficient * boundary_value / volumes[-1]
    return lwr, diaal, upr, soce


def apply_operator(op: tuple[np.ndarray, ...], field) :
    lwr, diaal, upr, soce = op
    relt = diaal * field + soce
    relt[1:] += lwr[1:] * field[:-1]
    relt[:-1] += upr[:-1] * field[1:]
    return relt


def cn_linear_solve(old, old_op: tuple[np.ndarray, ...],
                    new_op: tuple[np.ndarray, ...], capacity,
                    dt_s: float, theta: float) -> tuple[np.ndarray, float]:
    lwr, diaal, upr, soce = new_op
    rhs = capacity * old + dt_s * ((1.0 - theta) * apply_operator(old_op, old)
                                   + theta * soce)
    maix = np.zeros((3, old.size))
    maix[0, 1:] = -dt_s * theta * upr[:-1]
    maix[1] = capacity - dt_s * theta * diaal
    maix[2, :-1] = -dt_s * theta * lwr[1:]
    aer = solve_banded((1, 1), maix, rhs, check_finite=False)
    reual = (capacity * (aer - old) - dt_s * (
        (1.0 - theta) * apply_operator(old_op, old)
        + theta * apply_operator(new_op, aer)))
    scaled = float(np.max(np.abs(reual) / np.maximum(capacity, 1.0)))
    return aer, scaled


def advance(t0, c0, time0_s: float, dt_s: float,
            measured, geom: tuple[np.ndarray, ...], cfg: Config,
            thermal_locked: bool = False
            ) -> tuple[np.ndarray, np.ndarray, int, float, float]:
    env0 = environment(time0_s, measured, cfg)
    env1 = environment(time0_s + dt_s, measured, cfg)
    cap0 = volumetric_heat_capacity(c0)
    heat0 = spatial_operator(conductivity(c0), HEAT_TRANSFER, env0[0], geom)
    mass0 = spatial_operator(diffusivity(t0, c0), MASS_TRANSFER, env0[1], geom)
    tg, cg = t0.copy(), c0.copy()
    max_residual_t = 0.0
    max_residual_c = 0.0
    for iteration in range(1, cfg.max_iterations + 1):
        if thermal_locked:
            tn = np.full_like(t0, cfg.post_temperature_c)
            residual_t = 0.0
        else:
            cap1 = volumetric_heat_capacity(cg)
            heat1 = spatial_operator(conductivity(cg), HEAT_TRANSFER, env1[0], geom)
            tn, residual_t = cn_linear_solve(
                t0, heat0, heat1, 0.5 * (cap0 + cap1), dt_s, cfg.theta)
        mass1 = spatial_operator(diffusivity(tn, cg), MASS_TRANSFER, env1[1], geom)
        cn, residual_c = cn_linear_solve(
            c0, mass0, mass1, np.ones_like(c0), dt_s, cfg.theta)
        max_residual_t = max(max_residual_t, residual_t)
        max_residual_c = max(max_residual_c, residual_c)
        error = max(float(np.max(np.abs(tn - tg))) / 50.0,
                    float(np.max(np.abs(cn - cg))) / INITIAL_C)
        if not np.all(np.isfinite(tn)) or not np.all(np.isfinite(cn)):
            raise FloatingPointError("出现NaN或Inf")
        if np.min(cn) <= 0.0:
            raise FloatingPointError("出现非正含水率，请减小时间步")
        if error < cfg.tolerance:
            return tn, cn, iteration, max_residual_t, max_residual_c
        tg, cg = tn, cn
    raise RuntimeError(f"Picard迭代未在{cfg.max_iterations}次内收敛，误差={error:.3e}")


def locate_event(t0, c0, time0_s: float, dt_s: float,
                 measured, geom: tuple[np.ndarray, ...], cfg: Config,
                 thermal_locked: bool
                 ) -> dict:
    low, high = 0.0, dt_s
    te = t0.copy()
    ce = c0.copy()
    while high - low > cfg.event_tolerance_s:
        middle = 0.5 * (low + high)
        tm, cm, _, _, _ = advance(
            t0, c0, time0_s, middle, measured, geom, cfg, thermal_locked)
        if float(np.max(cm)) < THRESHOLD:
            high, te, ce = middle, tm, cm
        else:
            low = middle
    if not np.isfinite(ce).all() or float(np.max(ce)) >= THRESHOLD:
        te, ce, _, _, _ = advance(
            t0, c0, time0_s, high, measured, geom, cfg, thermal_locked)
    return {"time_s": time0_s + high, "bracket_s": [time0_s + low, time0_s + high],
            "temperature_C": te, "moisture_kg_kg": ce,
            "max_C": float(np.max(ce)), "argmax_node": int(np.argmax(ce)),
            "average_C": area_average(ce, geom[2]), "surface_C": float(ce[-1])}


