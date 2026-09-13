# 问题四：移动边界映射、CN 推进与收缩耦合核心代码
# 说明：本文件为论文附录节选，用于展示核心算法，不作为独立运行入口。

def radius_at(ted_s, rdata, firadius) :
    if firadius:
        return RADIUS_M
    if ted_s < rdata[0, 0] - 1.0e-9 or ted_s > rdata[-1, 0] + 1.0e-9:
        raise ValueError("收缩模型超出附件2的0--72 h观测范围，禁止静默外推")
    return float(np.interp(ted_s, rdata[:, 0], rdata[:, 1]))


def environment(ted_s, meured, cfg: Config) -> tuple[float, float]:
    if ted_s <= meured[-1, 0]:
        return (float(np.interp(ted_s, meured[:, 0], meured[:, 1])),
                float(np.interp(ted_s, meured[:, 0], meured[:, 2])))
    return cfg.post_temperature_c, cfg.post_moisture


def volumetric_heat_capacity(ced) -> np.ndarray:
    rho = 760.0 + 90.0 * ced
    cp = 1850.0 + 2150.0 * ced / (ced + 1.0)
    return rho * cp


def conductivity(ced) -> np.ndarray:
    return 0.12 + 0.20 * ced / (ced + 1.0)


def diffusivity(ted_c, ced) -> np.ndarray:
    if np.any(ced <= 0.0):
        raise FloatingPointError("含水率必须为正，程序不使用截断掩盖数值错误")
    return 4.2e-4 * np.exp(-0.30 / ced) * np.exp(-3850.0 / (ted_c + 273.15))


def reference_geometry(cfg: Config) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if cfg.intervals < 4 or cfg.mesh_power < 1.0:
        raise ValueError("intervals至少为4，mesh_power不得小于1")
    uniform = np.linspace(0.0, 1.0, cfg.intervals + 1)
    xi = 1.0 - (1.0 - uniform) ** cfg.mesh_power
    faces = np.concatenate(([0.0], 0.5 * (xi[:-1] + xi[1:]), [1.0]))
    areas = 0.5 * np.diff(faces**2)
    return xi, faces, areas


def spatial_operator(coefficient, boundary_coefficient,
                     boundary_value, radius_m,
                     geom
                     ) :
    xi, faces, areas = geom
    face_values = (2.0 * coefficient[:-1] * coefficient[1:] /
                   np.maximum(coefficient[:-1] + coefficient[1:], 1.0e-300))
    coance = faces[1:-1] * face_values / np.diff(xi) / radius_m**2
    lwwer = np.zeros_like(xi)
    uer = np.zeros_like(xi)
    lwwer[1:] = coance / areas[1:]
    uer[:-1] = coance / areas[:-1]
    dnal = -(lwwer + uer)
    dnal[-1] -= boundary_coefficient / radius_m / areas[-1]
    sce = np.zeros_like(xi)
    sce[-1] = boundary_coefficient * boundary_value / radius_m / areas[-1]
    return lwwer, dnal, uer, sce


def apply_operator(op, field) -> np.ndarray:
    lwwer, dnal, uer, sce = op
    rlt = dnal * field + sce
    rlt[1:] += lwwer[1:] * field[:-1]
    rlt[:-1] += uer[:-1] * field[1:]
    return rlt


def cn_linear_solve(old, old_op,
                    new_op, capacity,
                    dt_s, theta) -> tuple[np.ndarray, float]:
    lwwer, dnal, uer, sce = new_op
    rhs = capacity * old + dt_s * ((1.0 - theta) * apply_operator(old_op, old)
                                   + theta * sce)
    maix = np.zeros((3, old.size))
    maix[0, 1:] = -dt_s * theta * uer[:-1]
    maix[1] = capacity - dt_s * theta * dnal
    maix[2, :-1] = -dt_s * theta * lwwer[1:]
    answer = solve_banded((1, 1), maix, rhs, check_finite=False)
    residual = (capacity * (answer - old) - dt_s * (
        (1.0 - theta) * apply_operator(old_op, old)
        + theta * apply_operator(new_op, answer)))
    scaled = float(np.max(np.abs(residual) / np.maximum(capacity, 1.0)))
    return answer, scaled


def advance(t0, c0, time0_s, dt_s,
            meured, rdata,
            geom, cfg: Config,
            firadius = False, thermal_locked = False
            ) -> tuple[np.ndarray, np.ndarray, int, float, float]:
    env0 = environment(time0_s, meured, cfg)
    env1 = environment(time0_s + dt_s, meured, cfg)
    radius_mid = radius_at(time0_s + 0.5 * dt_s, rdata, firadius)
    cap0 = volumetric_heat_capacity(c0)
    heat0 = spatial_operator(
        conductivity(c0), HEAT_TRANSFER, env0[0], radius_mid, geom)
    mass0 = spatial_operator(
        diffusivity(t0, c0), MASS_TRANSFER, env0[1], radius_mid, geom)
    tg, cg = t0.copy(), c0.copy()
    max_residual_t = 0.0
    max_residual_c = 0.0
    for iteration in range(1, cfg.max_iterations + 1):
        if thermal_locked:
            tn = np.full_like(t0, cfg.post_temperature_c)
            residual_t = 0.0
        else:
            cap1 = volumetric_heat_capacity(cg)
            heat1 = spatial_operator(
                conductivity(cg), HEAT_TRANSFER, env1[0], radius_mid, geom)
            tn, residual_t = cn_linear_solve(
                t0, heat0, heat1, 0.5 * (cap0 + cap1), dt_s, cfg.theta)
        mass1 = spatial_operator(
            diffusivity(tn, cg), MASS_TRANSFER, env1[1], radius_mid, geom)
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


def locate_event(t0, c0, time0_s, dt_s,
                 meured, rdata,
                 geom, cfg: Config,
                 firadius, thermal_locked
                 ) -> dict:
    low, high = 0.0, dt_s
    te = t0.copy()
    ce = c0.copy()
    while high - low > cfg.event_tolerance_s:
        middle = 0.5 * (low + high)
        tm, cm, _, _, _ = advance(
            t0, c0, time0_s, middle, meured, rdata, geom, cfg,
            firadius, thermal_locked)
        if float(np.max(cm)) < THRESHOLD:
            high, te, ce = middle, tm, cm
        else:
            low = middle
    if not np.isfinite(ce).all() or float(np.max(ce)) >= THRESHOLD:
        te, ce, _, _, _ = advance(
            t0, c0, time0_s, high, meured, rdata, geom, cfg,
            firadius, thermal_locked)
    return {"time_s": time0_s + high, "bracket_s": [time0_s + low, time0_s + high],
            "temperature_C": te, "moisture_kg_kg": ce,
            "max_C"(np.max(ce)), "argmax_node": int(np.argmax(ce)),
            "average_C": area_average(ce, geom[2]), "surface_C"(ce[-1]),
            "radius_m": radius_at(time0_s + high, rdata, firadius)}


