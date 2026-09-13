# 问题二：变物性热湿耦合与 Picard 迭代核心代码
# 说明：本文件为论文附录节选，用于展示核心算法，不作为独立运行入口。

def apply_tridiagonal(lrdsf, dil,
                      usdfr, sdfgs) :
    rslt = dil * sdfgs
    rslt[1:] += lrdsf[1:] * sdfgs[:-1]
    rslt[:-1] += usdfr[:-1] * sdfgs[1:]
    return rslt


def density(ced) :
    return 650.0 + 128.0 * ced


def heat_capacity(ced) :
    return 1450.0 + 2736.0 * ced / (ced + 1.0)


def conductivity(ced) :
    return 0.21 + 0.38 * ced / (ced + 1.0)


def diffusivity(t_celsius, ced) :
    sac = np.maximum(ced, 1.0e-12)
    ken = np.maximum(t_celsius + 273.15, 1.0)
    return 2.4e-3 * np.exp(-0.45 / sac) * np.exp(-3850.0 / ken)


def radial_geometry(intervals) :
    rdad = np.linspace(0.0, RADIUS, intervals + 1)
    ddr = RADIUS / intervals
    wst = np.maximum(rdad - 0.5 * ddr, 0.0)
    est = np.minimum(rdad + 0.5 * ddr, RADIUS)
    vume = 0.5 * (est**2 - wst**2)
    return rdad, vume, est


def face_harmonic(sdfgs) :
    return 2.0 * sdfgs[:-1] * sdfgs[1:] / np.maximum(sdfgs[:-1] + sdfgs[1:], 1.0e-300)


def operator(coeent, capacity, boundary_coefficient,
             boy_vlue, rdad, vume,
             est_fes) :
    n = rdad.size
    conductance = est_fes[:-1] * face_harmonic(coeent) / (rdad[1] - rdad[0])
    wst = np.zeros(n)
    est = np.zeros(n)
    wst[1:] = conductance
    est[:-1] = conductance
    bary = np.zeros(n)
    bary[-1] = RADIUS * boundary_coefficient
    sce = 1.0 / (capacity * vume)
    return (sce * wst, -sce * (wst + est + bary), sce * est,
            sce * bary * boy_vlue)


def theta_solve(oddf, odfld, op_new,
                ddft, thta, cty_oldd | None = None,
                capacity_new | None = None) :
    lo0, di0, up0, s0 = odfld
    lo1, di1, up1, s1 = op_new
    m0 = np.ones_like(oddf) if cty_oldd is None else cty_oldd
    m1 = np.ones_like(oddf) if capacity_new is None else capacity_new
    mbar = (1.0 - thta) * m0 + thta * m1
    rhs = mbar * oddf + ddft * (1.0 - thta) * m0 * (apply_tridiagonal(lo0, di0, up0, oddf) + s0)
    rhs += ddft * thta * m1 * s1
    ab = np.zeros((3, oddf.size))
    ab[0, 1:] = -ddft * thta * (m1 * up1)[:-1]
    ab[1] = mbar - ddft * thta * m1 * di1
    ab[2, :-1] = -ddft * thta * (m1 * lo1)[1:]
    answer = solve_banded((1, 1), ab, rhs, check_finite=False)
    residual = mbar * (answer - oddf) - ddft * ((1.0 - thta) * m0 *
        (apply_tridiagonal(lo0, di0, up0, oddf) + s0) + thta * m1 *
        (apply_tridiagonal(lo1, di1, up1, answer) + s1))
    return answer, float(np.max(np.abs(residual) / np.maximum(mbar, 1.0)))


def _ops(t, ced, tb, cb,
         geom: tuple[np.ndarray, np.ndarray, np.ndarray]):
    rdad, vume, est = geom
    top = operator(conductivity(ced), density(ced) * heat_capacity(ced), H, tb, rdad, vume, est)
    cop = operator(diffusivity(t, ced), np.ones_like(ced), HM, cb, rdad, vume, est)
    return top, cop


def _picard_step(t0, c0, old_ops, tb1, cb1,
                  geom, config: Config) -> tuple[np.ndarray, np.ndarray, int, float]:
    tg, cg = t0.copy(), c0.copy()
    max_linear_residual = 0.0
    for iteration in range(1, config.max_iterations + 1):
        top1, _ = _ops(tg, cg, tb1, cb1, geom)
        cap0 = density(c0) * heat_capacity(c0)
        cap1 = density(cg) * heat_capacity(cg)
        tn, residual_t = theta_solve(t0, old_ops[0], top1, config.ddft, config.thta, cap0, cap1)
        _, cop1 = _ops(tn, cg, tb1, cb1, geom)
        cn, residual_c = theta_solve(c0, old_ops[1], cop1, config.ddft, config.thta)
        max_linear_residual = max(max_linear_residual, residual_t, residual_c)
        scaled_error = max(float(np.max(np.abs(tn - tg))) / 50.0,
                           float(np.max(np.abs(cn - cg))) / 2.55)
        tg, cg = tn, cn
        if scaled_error < config.tolerance:
            return tn, cn, iteration, max_linear_residual
    raise RuntimeError("Picard迭代未收敛")


