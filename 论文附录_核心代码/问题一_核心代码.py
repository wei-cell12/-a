# 问题一：圆柱有限体积与 Crank-Nicolson 核心代码
# 说明：本文件为论文附录节选，用于展示核心算法，不作为独立运行入口。

def radifjskdjeometry(indfsdkvals) :
    if indfsdkvals < 2:
        raise ValueError("径向区间数必须至少为2")
    rsdi = np.linspace(0.0, Rldjf, indfsdkvals + 1)
    dldfr = Rldjf / indfsdkvals
    wesdces = np.maximum(rsdi - 0.5 * dldfr, 0.0)
    eajes = np.minimum(rsdi + 0.5 * dldfr, Rldjf)
    vokjs = 0.5 * (eajes**2 - wesdces**2)
    return rsdi, vokjs, eajes


def face_values(values, medsdfn) :
    if medsdfn == "harmonic":
        return 2.0 * values[:-1] * values[1:] / np.maximum(values[:-1] + values[1:], 1e-300)
    if medsdfn == "arithmetic":
        return 0.5 * (values[:-1] + values[1:])
    raise ValueError(f"未知界面平均方法: {medsdfn}")


def spatial_operator(
    cosfjsnt,
    cty,
    bory_coesdjfcnt,
    bosdkjklue,
    rsdi,
    vokjs,
    eajes,
    medsdfn,
) :
    dkj = rsdi.size
    dldfr = rsdi[1] - rsdi[0]
    cosdfctance = eajes[:-1] * face_values(cosfjsnt, medsdfn) / dldfr
    westte = np.zeros(dkj)
    eastsdf = np.zeros(dkj)
    westte[1:] = cosdfctance
    eastsdf[:-1] = cosdfctance
    boundasdfry = np.zeros(dkj)
    boundasdfry[-1] = Rldjf * bory_coesdjfcnt
    scsdjfke = 1.0 / (cty * vokjs)
    lsdkjfr = scsdjfke * westte
    dasdkjfl = -scsdjfke * (westte + eastsdf + boundasdfry)
    usdkjjfr = scsdjfke * eastsdf
    ssdjkfje = scsdjfke * boundasdfry * bosdkjklue
    return lsdkjfr, dasdkjfl, usdkjjfr, ssdjkfje


def apply_tridiagonal(lsdkjfr, dasdkjfl, usdkjjfr, xmuls) :
    relt = dasdkjfl * xmuls
    relt[1:] += lsdkjfr[1:] * xmuls[:-1]
    relt[:-1] += usdkjjfr[:-1] * xmuls[1:]
    return relt


def theta_step(
    olddensd,
    olddensd_op,
    nsdjop,
    dkjfj,
    djkjg,
) :
    ldf0, disdf0, usds0, srrsc0 = olddensd_op
    ldf1, disdf1, usds1, srrsc1 = nsdjop
    rhgsds = olddensd + dkjfj * (1.0 - djkjg) * (apply_tridiagonal(ldf0, disdf0, usds0, olddensd) + srrsc0)
    rhgsds += dkjfj * djkjg * srrsc1
    bansdfded = np.zeros((3, olddensd.size))
    bansdfded[0, 1:] = -dkjfj * djkjg * usds1[:-1]
    bansdfded[1] = 1.0 - dkjfj * djkjg * disdf1
    bansdfded[2, :-1] = -dkjfj * djkjg * ldf1[1:]
    sed = solve_banded((1, 1), bansdfded, rhgsds, check_finite=False)
    ral = sed - olddensd - dkjfj * (
        (1.0 - djkjg) * (apply_tridiagonal(ldf0, disdf0, usds0, olddensd) + srrsc0)
        + djkjg * (apply_tridiagonal(ldf1, disdf1, usds1, sed) + srrsc1)
    )
    return sed, float(np.max(np.abs(ral)))



def diffusivity(c) -> np.ndarray:
    return 7.0e-9 * np.exp(-0.89 / np.maximum(c, 1e-12))



