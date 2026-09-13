#!/usr/bin/env python3
"""A题问题二：变物性热湿耦合圆柱模型（Crank–Nicolson + Picard）。"""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
import zipfile
from dataclasses import dataclass
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import Akima1DInterpolator, CubicSpline, PchipInterpolator
from scipy.linalg import solve_banded
from scipy.optimize import root

import xlsx_helpers as xlsx_utils
from utils.plot_style import COLOR_SEQUENCE, PALETTE, apply_publication_style


ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
RADIUS = 0.02
H = 25.0
HM = 8.0e-7
T_INITIAL = 28.0
C_INITIAL = 2.55
END_TIME = 10800.0
KEY_TIMES = np.arange(1800.0, END_TIME + 1.0, 1800.0)
KEY_RADII_CM = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
OUTPUT_RADII_CM = np.round(np.arange(0.0, 2.0001, 0.1), 10)


# 根据离环数据造法时插数。
def make_interpolator(x: np.ndarray, y: np.ndarray, method: str):
    if method == "linear":
        return lambda z: np.interp(z, x, y)
    if method == "pchip":
        return PchipInterpolator(x, y, extrapolate=True)
    if method == "cubic":
        return CubicSpline(x, y, bc_type="natural", extrapolate=True)
    if method == "akima":
        return Akima1DInterpolator(x, y)
    raise ValueError(f"未知插值方法: {method}")


# 将三矩阵算用于给定状态向量。
def apply_tridiagonal(lower: np.ndarray, diagonal: np.ndarray,
                      upper: np.ndarray, values: np.ndarray) -> np.ndarray:
    result = diagonal * values
    result[1:] += lower[1:] * values[:-1]
    result[:-1] += upper[:-1] * values[1:]
    return result


@dataclass(frozen=True)
class Config:
    intervals: int = 200
    dt: float = 0.5
    end_time: float = END_TIME
    store_interval: float = 1.0
    interpolation: str = "linear"
    theta: float = 0.5
    solver: str = "picard"
    tolerance: float = 1.0e-10
    max_iterations: int = 30


# 根含水率计算药材密度。
def density(c: np.ndarray) -> np.ndarray:
    return 650.0 + 128.0 * c


# 据含水率计算药材比热容。
def heat_capacity(c: np.ndarray) -> np.ndarray:
    return 1450.0 + 2736.0 * c / (c + 1.0)


# 据含水率计算药材导热系数。
def conductivity(c: np.ndarray) -> np.ndarray:
    return 0.21 + 0.38 * c / (c + 1.0)


# 根温度和含水率计算水分有效扩散系数。
def diffusivity(t_celsius: np.ndarray, c: np.ndarray) -> np.ndarray:
    safe_c = np.maximum(c, 1.0e-12)
    kelvin = np.maximum(t_celsius + 273.15, 1.0)
    return 2.4e-3 * np.exp(-0.45 / safe_c) * np.exp(-3850.0 / kelvin)


# 造圆柱径向网格、控制体面积和界面位置。
def radial_geometry(intervals: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    r = np.linspace(0.0, RADIUS, intervals + 1)
    dr = RADIUS / intervals
    west = np.maximum(r - 0.5 * dr, 0.0)
    east = np.minimum(r + 0.5 * dr, RADIUS)
    volume = 0.5 * (east**2 - west**2)
    return r, volume, east


# 计算相邻节点物性参数调和平均界面值。
def face_harmonic(values: np.ndarray) -> np.ndarray:
    return 2.0 * values[:-1] * values[1:] / np.maximum(values[:-1] + values[1:], 1.0e-300)


# 根据物性和边界数组装径向离散算子。
def operator(coefficient: np.ndarray, capacity: np.ndarray, boundary_coefficient: float,
             boundary_value: float, r: np.ndarray, volume: np.ndarray,
             east_faces: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    n = r.size
    coctance = east_faces[:-1] * face_harmonic(coefficient) / (r[1] - r[0])
    west = np.zeros(n)
    east = np.zeros(n)
    west[1:] = coctance
    east[:-1] = coctance
    boundary = np.zeros(n)
    boundary[-1] = RADIUS * boundary_coefficient
    scale = 1.0 / (capacity * volume)
    return (scale * west, -scale * (west + east + boundary), scale * east,
            scale * boundary * boundary_value)


# 求解一次 theta 格式离散后三对角线性方程组。
def theta_solve(old: np.ndarray, op_old: tuple[np.ndarray, ...], op_new: tuple[np.ndarray, ...],
                dt: float, theta: float, capacity_old: np.ndarray | None = None,
                capacity_new: np.ndarray | None = None) -> tuple[np.ndarray, float]:
    lo0, di0, up0, s0 = op_old
    lo1, di1, up1, s1 = op_new
    m0 = np.ones_like(old) if capacity_old is None else capacity_old
    m1 = np.ones_like(old) if capacity_new is None else capacity_new
    mbar = (1.0 - theta) * m0 + theta * m1
    rhs = mbar * old + dt * (1.0 - theta) * m0 * (apply_tridiagonal(lo0, di0, up0, old) + s0)
    rhs += dt * theta * m1 * s1
    ab = np.zeros((3, old.size))
    ab[0, 1:] = -dt * theta * (m1 * up1)[:-1]
    ab[1] = mbar - dt * theta * m1 * di1
    ab[2, :-1] = -dt * theta * (m1 * lo1)[1:]
    answer = solve_banded((1, 1), ab, rhs, check_finite=False)
    residual = mbar * (answer - old) - dt * ((1.0 - theta) * m0 *
        (apply_tridiagonal(lo0, di0, up0, old) + s0) + theta * m1 *
        (apply_tridiagonal(lo1, di1, up1, answer) + s1))
    return answer, float(np.max(np.abs(residual) / np.maximum(mbar, 1.0)))


# 读取附件中的烘房温与水分边界数据。
def load_environment() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    path = DATA / "附件1.xlsx"
    rows = xlsx_utils.xlsx_rows(path)
    values = np.asarray([[float(v) for v in row[:3]] for row in rows[1:]], dtype=float)
    if values.shape != (241, 3) or values[0, 0] != 0 or values[-1, 0] != 14400:
        raise ValueError("附件1结构或时间范围不符合题目说明")
    return values[:, 0], values[:, 1], values[:, 2]


# 根据当前温湿度场生热量与水分方程的离散算子。
def _ops(t: np.ndarray, c: np.ndarray, tb: float, cb: float,
         geom: tuple[np.ndarray, np.ndarray, np.ndarray]):
    r, volume, east = geom
    top = operator(conductivity(c), density(c) * heat_capacity(c), H, tb, r, volume, east)
    cop = operator(diffusivity(t, c), np.ones_like(c), HM, cb, r, volume, east)
    return top, cop


# 使用 Picard 迭代完成一步非线性热耦合求解。
def _picard_step(t0: np.ndarray, c0: np.ndarray, old_ops, tb1: float, cb1: float,
                  geom, config: Config) -> tuple[np.ndarray, np.ndarray, int, float]:
    tg, cg = t0.copy(), c0.copy()
    maliar_idual = 0.0
    for iteration in range(1, config.max_iterations + 1):
        top1, _ = _ops(tg, cg, tb1, cb1, geom)
        cap0 = density(c0) * heat_capacity(c0)
        cap1 = density(cg) * heat_capacity(cg)
        tn, residual_t = theta_solve(t0, old_ops[0], top1, config.dt, config.theta, cap0, cap1)
        _, cop1 = _ops(tn, cg, tb1, cb1, geom)
        cn, residual_c = theta_solve(c0, old_ops[1], cop1, config.dt, config.theta)
        maliar_idual = max(maliar_idual, residual_t, residual_c)
        scaled_error = max(float(np.max(np.abs(tn - tg))) / 50.0,
                           float(np.max(np.abs(cn - cg))) / 2.55)
        tg, cg = tn, cn
        if scaled_error < config.tolerance:
            return tn, cn, iteration, maliar_idual
    raise RuntimeError("Picard迭代未收敛")


# 计算非线性离方程在当前迭代状态下的残差向量。
def _residual_vector(x: np.ndarray, t0: np.ndarray, c0: np.ndarray, old_ops,
                     tb1: float, cb1: float, geom, config: Config) -> np.ndarray:
    n = t0.size
    tn, cn = x[:n], x[n:]
    top1, cop1 = _ops(tn, cn, tb1, cb1, geom)
    result = []
    for old, new, op0, op1, m0, m1 in (
        (t0, tn, old_ops[0], top1, density(c0) * heat_capacity(c0), density(cn) * heat_capacity(cn)),
        (c0, cn, old_ops[1], cop1, np.ones_like(c0), np.ones_like(cn))):
        lo0, di0, up0, s0 = op0
        lo1, di1, up1, s1 = op1
        mbar = (1.0 - config.theta) * m0 + config.theta * m1
        result.append((mbar * (new - old) - config.dt * ((1.0 - config.theta) * m0 *
            (apply_tridiagonal(lo0, di0, up0, old) + s0) + config.theta * m1 *
            (apply_tridiagonal(lo1, di1, up1, new) + s1))) / np.maximum(mbar, 1.0))
    return np.concatenate(result)

# 使用牛顿法修正一步非线热湿耦合状态。
def _newton_step(t0: np.ndarray, c0: np.ndarray, old_ops, tb1: float, cb1: float,
                  geom, config: Config) -> tuple[np.ndarray, np.ndarray, int, float]:
    x0 = np.concatenate((t0, c0))
    solved = root(_residual_vector, x0, args=(t0, c0, old_ops, tb1, cb1, geom, config),
                  method="hybr", options={"xtol": config.tolerance, "maxfev": 300})
    if not solved.success:
        raise RuntimeError(f"Newton迭代未收敛: {solved.message}")
    n = t0.size
    residual = float(np.max(np.abs(_residual_vector(
        solved.x, t0, c0, old_ops, tb1, cb1, geom, config))))
    return solved.x[:n], solved.x[n:], int(solved.nfev), residual


# 按圆柱截面控制体重计算场变量的面积平均值。
def area_average(field: np.ndarray, volume: np.ndarray) -> float:
    return float(np.dot(field, volume) / (0.5 * RADIUS**2))


# 按给定配置执行完整程的热湿耦合数值模拟。
def simulate(env_time: np.ndarray, env_t: np.ndarray, env_c: np.ndarray, config: Config) -> dict:
    steps = int(round(config.end_time / config.dt))
    stride = int(round(config.store_interval / config.dt))
    if not math.isclose(steps * config.dt, config.end_time) or not math.isclose(stride * config.dt, config.store_interval):
        raise ValueError("结束时间和输出间隔必须能被时间步长整除")
    measured_t = make_interpolator(env_time, env_t, config.interpolation)
    measured_c = make_interpolator(env_time, env_c, config.interpolation)
    # 附件1覆盖0~4 h；其后进入题述恒温干燥阶段，取设定值50 ℃、0.05 kg/kg。
    boundary_t = lambda z: np.where(np.asarray(z) <= env_time[-1], measured_t(z), 50.0)
    boundary_c = lambda z: np.where(np.asarray(z) <= env_time[-1], measured_c(z), 0.05)
    geom = radial_geometry(config.intervals)
    r, volume, _ = geom
    t = np.full(r.size, T_INITIAL)
    c = np.full(r.size, C_INITIAL)
    out_time = np.arange(0.0, config.end_time + 0.5 * config.store_interval, config.store_interval)
    out_t = np.empty((out_time.size, r.size)); out_c = np.empty_like(out_t)
    out_t[0], out_c[0] = t, c
    iterations = np.empty(steps, dtype=np.int16)
    linear_residual = 0.0
    heat_balance = np.empty(steps); mass_balance = np.empty(steps)
    store_index = 1
    started = time.perf_counter()
    for step in range(1, steps + 1):
        time0, time1 = (step - 1) * config.dt, step * config.dt
        tb0, tb1 = float(boundary_t(time0)), float(boundary_t(time1))
        cb0, cb1 = float(boundary_c(time0)), float(boundary_c(time1))
        old_ops = _ops(t, c, tb0, cb0, geom)
        old_c_avg = area_average(c, volume)
        if config.solver == "picard":
            tn, cn, count, residual = _picard_step(t, c, old_ops, tb1, cb1, geom, config)
        elif config.solver == "newton":
            tn, cn, count, residual = _newton_step(t, c, old_ops, tb1, cb1, geom, config)
        else:
            raise ValueError(f"未知非线性求解器: {config.solver}")
        iterations[step - 1] = count
        linear_residual = max(linear_residual, residual)
        # 变体积热容下，用控方程右端积分与左端面积平均增量作数值诊断。
        cap0 = density(c) * heat_capacity(c)
        cap1 = density(cn) * heat_capacity(cn)
        capbar = (1.0 - config.theta) * cap0 + config.theta * cap1
        stored_heat = float(np.dot(volume * capbar, tn - t))
        supplied_heat = config.dt * RADIUS * H * ((1.0 - config.theta) * (tb0 - t[-1]) +
                                                  config.theta * (tb1 - tn[-1]))
        heat_balance[step - 1] = (stored_heat - supplied_heat) / float(np.dot(volume, capbar))
        mass_flux = config.dt * ((1.0 - config.theta) * 2 * HM / RADIUS * (cb0 - c[-1]) +
                                 config.theta * 2 * HM / RADIUS * (cb1 - cn[-1]))
        mass_balance[step - 1] = area_average(cn, volume) - old_c_avg - mass_flux
        t, c = tn, cn
        if step % stride == 0:
            out_t[store_index], out_c[store_index] = t, c
            store_index += 1
    return {"config": config, "time": out_time, "radius_m": r, "temperature": out_t,
            "moisture": out_c, "iterations": iterations, "runtime_s": time.perf_counter() - started,
            "max_nonlinear_iterations": int(iterations.max()), "mean_nonlinear_iterations": float(iterations.mean()),
            "max_algebraic_residual": linear_residual,
            "heat_step_balance_max_degC_equivalent": float(np.max(np.abs(heat_balance))),
            "mass_step_balance_max": float(np.max(np.abs(mass_balance)))}


# 从完整时空结果中插值得指定时刻和半径处的数值。
def extract(result: dict, times: np.ndarray, radii_cm: np.ndarray, key: str) -> np.ndarray:
    source_t = np.asarray(result["time"]); source_r = np.asarray(result["radius_m"]) * 100.0
    field = np.asarray(result[key]); rows = []
    for target in times:
        index = int(np.argmin(np.abs(source_t - target)))
        if not math.isclose(float(source_t[index]), float(target), abs_tol=1e-9):
            raise ValueError(f"缺少输出时刻 {target}")
        rows.append(np.interp(radii_cm, source_r, field[index]))
    return np.asarray(rows)


# 汇总数值解的稳定性边界响应和物理合理性检查。
def validate(result: dict, env_t: np.ndarray) -> dict:
    t, c = np.asarray(result["temperature"]), np.asarray(result["moisture"])
    temperature_increments = np.diff(t, axis=1)
    monotone_rows = np.all(temperature_increments >= -2e-7, axis=1)
    checks = {
        "finite": bool(np.all(np.isfinite(t)) and np.all(np.isfinite(c))),
        "temperature_min": float(t.min()), "temperature_max": float(t.max()),
        "moisture_min": float(c.min()), "moisture_max": float(c.max()),
        "radial_temperature_monotone": bool(np.all(monotone_rows)),
        "radial_temperature_monotone_time_fraction": float(np.mean(monotone_rows)),
        "maximum_temperature_reverse_increment_degC": float(max(0.0, -np.min(temperature_increments))),
        "radial_moisture_monotone": bool(np.all(np.diff(c, axis=1) <= 2e-7)),
        "max_nonlinear_iterations": result["max_nonlinear_iterations"],
        "mean_nonlinear_iterations": result["mean_nonlinear_iterations"],
        "max_algebraic_residual": result["max_algebraic_residual"],
        "heat_step_balance_max_degC_equivalent": result["heat_step_balance_max_degC_equivalent"],
        "mass_step_balance_max": result["mass_step_balance_max"],
        "runtime_s": result["runtime_s"],
        "radial_intervals": result["config"].intervals, "dt_s": result["config"].dt,
        "theta": result["config"].theta, "tolerance": result["config"].tolerance,
        "end_time_s": result["config"].end_time,
    }
    checks["passed"] = bool(checks["finite"] and checks["temperature_min"] >= T_INITIAL - 1e-7 and
        checks["temperature_max"] <= float(np.max(env_t)) + 1e-6 and checks["moisture_min"] >= 0 and
        checks["moisture_max"] <= C_INITIAL + 1e-7 and checks["radial_moisture_monotone"] and
        checks["max_algebraic_residual"] < 1e-7 and
        checks["heat_step_balance_max_degC_equivalent"] < 1e-8 and checks["mass_step_balance_max"] < 1e-8)
    if not checks["passed"]:
        raise AssertionError(f"物理或数值检查失败: {checks}")
    return checks


# 将时间、半径场变量矩阵写入 CSV 文件。
def write_csv_matrix(path: Path, times: np.ndarray, radii_cm: np.ndarray, values: np.ndarray,
                     time_label: str = "时间/s", decimals: int = 6) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.writer(stream); writer.writerow([time_label, *[f"{x:g} cm" for x in radii_cm]])
        for tm, row in zip(times, values):
            writer.writerow([f"{tm:g}", *[f"{v:.{decimals}f}" for v in row]])


# 根据结果矩阵成 Excel 工作表的 XML 内容。
def build_sheet_xml(times: np.ndarray, values: np.ndarray, numeric_style: int) -> bytes:
    rows = []
    cells = ['<c r="A1" s="2" t="inlineStr"><is><t>时间/s\\到药材中心的距离/cm</t></is></c>']
    for col, radius in enumerate(OUTPUT_RADII_CM, 2):
        cells.append(f'<c r="{xlsx_utils.excel_column(col)}1" s="2"><v>{radius:g}</v></c>')
    rows.append(f'<row r="1" spans="1:22">{"".join(cells)}</row>')
    for row_num, (tm, row) in enumerate(zip(times, values), 2):
        cells = [f'<c r="A{row_num}" s="1"><v>{tm:g}</v></c>']
        for col, value in enumerate(row, 2):
            cells.append(f'<c r="{xlsx_utils.excel_column(col)}{row_num}" s="{numeric_style}"><v>{value:.6f}</v></c>')
        rows.append(f'<row r="{row_num}" spans="1:22">{"".join(cells)}</row>')
    last = len(times) + 1
    xml = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        f'<dimension ref="A1:V{last}"/><sheetViews><sheetView workbookViewId="0"/></sheetViews>'
        '<sheetFormatPr defaultColWidth="8.75" defaultRowHeight="14.1"/>'
        '<cols><col min="1" max="1" width="24" customWidth="1"/><col min="2" max="22" width="10" customWidth="1"/></cols>'
        f'<sheetData>{"".join(rows)}</sheetData><pageMargins left="0.75" right="0.75" top="1" bottom="1" header="0.5" footer="0.5"/></worksheet>')
    return xml.encode("utf-8")


# 将温湿度计算果写入 Excel 模板。
def write_xlsx(template: Path, output: Path, times: np.ndarray, temperature: np.ndarray, moisture: np.ndarray) -> None:
    output.parent.mkdir(parents=True, exist_ok=True); temporary = output.with_suffix(".tmp.xlsx")
    with zipfile.ZipFile(template) as source, zipfile.ZipFile(temporary, "w", zipfile.ZIP_DEFLATED) as target:
        styles, style = xlsx_utils.add_four_decimal_style(source.read("xl/styles.xml"))
        replacements = {"xl/styles.xml": styles, "xl/worksheets/sheet1.xml": build_sheet_xml(times, temperature, style),
                        "xl/worksheets/sheet2.xml": build_sheet_xml(times, moisture, style)}
        for item in source.infolist(): target.writestr(item, replacements.get(item.filename, source.read(item.filename)))
    temporary.replace(output)


# 将图形同时导为论文所需的位图量格式。
def export_plot(fig, name: str, size=(6.4, 4.0)) -> None:
    fig.set_size_inches(*size)
    fig.savefig(FIGURES / f"{name}.png", dpi=300, bbox_inches=None)
    fig.savefig(FIGURES / f"{name}.svg", bbox_inches=None)
    fig.savefig(FIGURES / f"{name}.pdf", bbox_inches=None)
    plt.close(fig)


# 根据正式计算结果生成论文所需图件。
def make_figures(env_time, env_t, env_c, result, convergence_rows, solver_rows) -> None:
    apply_publication_style()
    plt.rcParams.update({"font.family": "sans-serif",
                         "font.sans-serif": ["Microsoft YaHei", "SimHei", "DejaVu Sans"],
                         "font.size": 8, "svg.fonttype": "none", "pdf.fonttype": 42})
    FIGURES.mkdir(parents=True, exist_ok=True)
    hours = env_time / 3600
    fig, ax = plt.subplots(layout="constrained"); ax.plot(hours[env_time <= END_TIME], env_t[env_time <= END_TIME], color=PALETTE["contrast"])
    ax.set(xlabel="时间（h）", ylabel="烘房温度（℃）"); export_plot(fig, "raw_q2_boundary_temperature")
    fig, ax = plt.subplots(layout="constrained"); ax.plot(hours[env_time <= END_TIME], env_c[env_time <= END_TIME], color=PALETTE["primary"])
    ax.set(xlabel="时间（h）", ylabel="烘房水分浓度（kg/kg）"); export_plot(fig, "raw_q2_boundary_moisture")
    sample_c = np.linspace(0.2, 2.55, 200); sample_t = np.array([28., 40., 50.])
    fig, ax = plt.subplots(layout="constrained")
    for i, tv in enumerate(sample_t): ax.plot(sample_c, diffusivity(np.full_like(sample_c, tv), sample_c), label=f"{tv:g} ℃", color=COLOR_SEQUENCE[i])
    ax.set(xlabel="药材水分浓度（kg/kg）", ylabel="扩散系数（m²/s）"); ax.legend(title="温度"); ax.ticklabel_format(axis="y", style="sci", scilimits=(0,0)); export_plot(fig, "raw_q2_diffusivity")

    fig, ax = plt.subplots(layout="constrained")
    for tm, color in zip(KEY_TIMES, COLOR_SEQUENCE):
        idx = int(tm); ax.plot(result["radius_m"]*100, result["temperature"][idx], color=color, label=f"{tm/3600:g} h")
    ax.set(xlabel="到药材中心的距离（cm）", ylabel="温度（℃）"); ax.legend(ncol=2); export_plot(fig, "process_q2_temperature_profiles")
    fig, ax = plt.subplots(layout="constrained")
    for tm, color in zip(KEY_TIMES, COLOR_SEQUENCE):
        idx = int(tm); ax.plot(result["radius_m"]*100, result["moisture"][idx], color=color, label=f"{tm/3600:g} h")
    ax.set(xlabel="到药材中心的距离（cm）", ylabel="水分浓度（kg/kg）"); ax.legend(ncol=2); export_plot(fig, "process_q2_moisture_profiles")
    fig, axes = plt.subplots(2, 3, layout="constrained")
    for col, (study, xlabel) in enumerate((("网格", "径向区间数 N"), ("时间步", "时间步长（s）"))):
        rows = [row for row in convergence_rows if row[0] == study]
        x = np.asarray([row[1] for row in rows], float)
        axes[0,col].loglog(x, [row[2] for row in rows], "o-", color=PALETTE["primary"])
        axes[1,col].loglog(x, [row[3] for row in rows], "s-", color=PALETTE["secondary"])
        axes[0,col].set(xlabel=xlabel, ylabel="温度最大差（℃）")
        axes[1,col].set(xlabel=xlabel, ylabel="水分最大差（kg/kg）")
        if study == "时间步": axes[0,col].invert_xaxis(); axes[1,col].invert_xaxis()
    names=[x[0] for x in solver_rows]
    axes[0,2].bar(names,[x[1] for x in solver_rows],color=[PALETTE["primary"],PALETTE["secondary"]]); axes[0,2].set(ylabel="300 s 运行时间（s）")
    axes[1,2].bar(names,[x[2] for x in solver_rows],color=[PALETTE["primary"],PALETTE["secondary"]]); axes[1,2].set(ylabel="平均迭代/调用次数")
    export_plot(fig, "process_q2_numerical_comparison", (8.0,5.4))

    time_h = result["time"]/3600; radius_cm=result["radius_m"]*100
    fig, ax = plt.subplots(layout="constrained"); mesh=ax.pcolormesh(time_h,radius_cm,result["temperature"].T,shading="auto",cmap="coolwarm",rasterized=True); fig.colorbar(mesh,ax=ax,label="温度（℃）"); ax.set(xlabel="时间（h）",ylabel="到药材中心的距离（cm）"); export_plot(fig,"result_q2_temperature_field")
    fig, ax = plt.subplots(layout="constrained"); mesh=ax.pcolormesh(time_h,radius_cm,result["moisture"].T,shading="auto",cmap="YlGnBu",rasterized=True); fig.colorbar(mesh,ax=ax,label="水分浓度（kg/kg）"); ax.set(xlabel="时间（h）",ylabel="到药材中心的距离（cm）"); export_plot(fig,"result_q2_moisture_field")
    fig, axes=plt.subplots(2,1,sharex=True,layout="constrained")
    for radius,color in zip(KEY_RADII_CM,COLOR_SEQUENCE):
        idx=int(round(radius/2*config_intervals(result))); axes[0].plot(time_h,result["temperature"][:,idx],color=color,label=f"{radius:g} cm"); axes[1].plot(time_h,result["moisture"][:,idx],color=color)
    axes[0].set(ylabel="温度（℃）"); axes[0].legend(ncol=3); axes[1].set(xlabel="时间（h）",ylabel="水分浓度（kg/kg）"); export_plot(fig,"result_q2_selected_positions",(7.2,4.8))


# 从模拟结果中读取实际采用的径向网格数。
def config_intervals(result: dict) -> int:
    return len(result["radius_m"])-1


# 直接读取已保存结果重新生成图件。
def figures_from_saved_outputs() -> None:
    """用已完成的1 s×0.1 cm正式输出重绘图片，避免重复数值计算。"""
    env_time, env_t, env_c = load_environment()
    t_data = np.genfromtxt(RESULTS/"q2_full_temperature.csv", delimiter=",", skip_header=1, encoding="utf-8-sig")
    c_data = np.genfromtxt(RESULTS/"q2_full_moisture.csv", delimiter=",", skip_header=1, encoding="utf-8-sig")
    times = np.concatenate(([0.0], t_data[:, 0]))
    temperature = np.vstack((np.full(OUTPUT_RADII_CM.size, T_INITIAL), t_data[:, 1:]))
    moisture = np.vstack((np.full(OUTPUT_RADII_CM.size, C_INITIAL), c_data[:, 1:]))
    result = {"time": times, "radius_m": OUTPUT_RADII_CM/100.0,
              "temperature": temperature, "moisture": moisture}
    convergence=[]
    with (RESULTS/"q2_convergence.csv").open(encoding="utf-8-sig") as f:
        for row in list(csv.reader(f))[1:]: convergence.append([row[0],float(row[1]),float(row[2]),float(row[3])])
    solver_rows=[]
    with (RESULTS/"q2_nonlinear_solver_comparison.csv").open(encoding="utf-8-sig") as f:
        for row in list(csv.reader(f))[1:3]: solver_rows.append([row[0],float(row[1]),float(row[2])])
    make_figures(env_time,env_t,env_c,result,convergence,solver_rows)


# 执行问题二正式求解、导出、验证和绘图流程。
def run_full(no_figures: bool = False) -> dict:
    RESULTS.mkdir(parents=True, exist_ok=True); FIGURES.mkdir(parents=True, exist_ok=True)
    env_time, env_t, env_c = load_environment()
    formal_config = Config()
    formal = simulate(env_time, env_t, env_c, formal_config); checks = validate(formal, env_t[env_time <= END_TIME])
    output_t = extract(formal, np.arange(1.0, END_TIME+1.0), OUTPUT_RADII_CM, "temperature")
    output_c = extract(formal, np.arange(1.0, END_TIME+1.0), OUTPUT_RADII_CM, "moisture")
    key_t = extract(formal, KEY_TIMES, KEY_RADII_CM, "temperature"); key_c = extract(formal, KEY_TIMES, KEY_RADII_CM, "moisture")
    write_csv_matrix(RESULTS/"q2_table3_temperature.csv",KEY_TIMES/3600.0,KEY_RADII_CM,key_t,"时间/h",4)
    write_csv_matrix(RESULTS/"q2_table4_moisture.csv",KEY_TIMES/3600.0,KEY_RADII_CM,key_c,"时间/h",4)
    write_csv_matrix(RESULTS/"q2_full_temperature.csv",np.arange(1.0,END_TIME+1.0),OUTPUT_RADII_CM,output_t)
    write_csv_matrix(RESULTS/"q2_full_moisture.csv",np.arange(1.0,END_TIME+1.0),OUTPUT_RADII_CM,output_c)
    write_xlsx(DATA/"附件3"/"result2.xlsx",RESULTS/"result2.xlsx",np.arange(1.0,END_TIME+1.0),output_t,output_c)

    convergence=[]
    grid_values={}
    for intervals in (50,100,200,400):
        candidate=formal if intervals == 200 else simulate(env_time,env_t,env_c,Config(intervals=intervals,dt=0.5,store_interval=1.0))
        grid_values[intervals]=(extract(candidate,KEY_TIMES,KEY_RADII_CM,"temperature"),extract(candidate,KEY_TIMES,KEY_RADII_CM,"moisture"))
    for intervals in (50,100,200):
        convergence.append(["网格",intervals,float(np.max(np.abs(grid_values[intervals][0]-grid_values[400][0]))),float(np.max(np.abs(grid_values[intervals][1]-grid_values[400][1])))])
    time_values={}
    for dt in (2.0,1.0,0.5,0.25):
        candidate=formal if dt == 0.5 else simulate(
            env_time, env_t, env_c, Config(intervals=200, dt=dt, store_interval=max(1.0, dt)))
        time_values[dt]=(extract(candidate,KEY_TIMES,KEY_RADII_CM,"temperature"),extract(candidate,KEY_TIMES,KEY_RADII_CM,"moisture"))
    for dt in (2.0,1.0,0.5):
        convergence.append(["时间步",dt,float(np.max(np.abs(time_values[dt][0]-time_values[0.25][0]))),float(np.max(np.abs(time_values[dt][1]-time_values[0.25][1])))])
    with (RESULTS/"q2_convergence.csv").open("w",encoding="utf-8-sig",newline="") as f:
        w=csv.writer(f);w.writerow(["加密类型","N或dt","温度最大差/℃","水分浓度最大差/(kg/kg)"]);w.writerows(convergence)

    short_picard=simulate(env_time,env_t,env_c,Config(intervals=20,dt=1.0,end_time=300,solver="picard"))
    short_newton=simulate(env_time,env_t,env_c,Config(intervals=20,dt=1.0,end_time=300,solver="newton"))
    solver_difference=max(float(np.max(np.abs(short_picard["temperature"]-short_newton["temperature"]))),float(np.max(np.abs(short_picard["moisture"]-short_newton["moisture"]))))
    solver_rows=[["Picard（300 s）",short_picard["runtime_s"],short_picard["mean_nonlinear_iterations"]],["Newton（300 s）",short_newton["runtime_s"],short_newton["mean_nonlinear_iterations"]]]
    with (RESULTS/"q2_nonlinear_solver_comparison.csv").open("w",encoding="utf-8-sig",newline="") as f:
        w=csv.writer(f);w.writerow(["方法","300s运行时间/s","平均迭代次数或函数调用次数"]);w.writerows(solver_rows);w.writerow(["两法最大结果差",solver_difference,""])
    summary={"formal_config":vars(formal_config),"validation":checks,"key_temperature":key_t.tolist(),"key_moisture":key_c.tolist(),"nonlinear_solver_max_difference_300s":solver_difference,"convergence":convergence}
    (RESULTS/"q2_validation_summary.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding="utf-8")
    if not no_figures: make_figures(env_time,env_t,env_c,formal,convergence,solver_rows)
    return summary


# 组织当前脚本的完整执行流程并返回运行状态。
def main() -> int:
    parser=argparse.ArgumentParser(); parser.add_argument("--smoke",action="store_true");parser.add_argument("--no-figures",action="store_true");parser.add_argument("--figures-only",action="store_true");args=parser.parse_args()
    env_time,env_t,env_c=load_environment()
    if args.figures_only:
        figures_from_saved_outputs(); return 0
    if args.smoke:
        result=simulate(env_time,env_t,env_c,Config(intervals=20,dt=1.0,end_time=60.0))
        print(json.dumps(validate(result,env_t[env_time<=60]),ensure_ascii=False,indent=2));return 0
    summary=run_full(args.no_figures);print(json.dumps(summary["validation"],ensure_ascii=False,indent=2));return 0


if __name__ == "__main__":
    raise SystemExit(main())
