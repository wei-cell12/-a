#!/usr/bin/env python3
"""问题四：收缩圆柱参考域上的非线性热湿耦合求解器。

空间采用节点中心有限体积法，时间采用 Crank--Nicolson 格式，
每一步使用 Picard 迭代同步更新物性。半径由附件2线性插值，
停止条件检查全部径向节点；可切换为附录4固定半径控制模型。
"""

from __future__ import annotations

import argparse
import json
import math
import time
import zipfile
from dataclasses import asdict, dataclass
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np
from scipy.linalg import solve_banded


# 不使用 resolve()：Windows 的 subst 盘符可规避旧版 Python 对中文长路径的误解码。
ROOT = Path(__file__).absolute().parent.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
RADIUS_M = 0.02
HEAT_TRANSFER = 25.0
MASS_TRANSFER = 8.0e-7
INITIAL_T_C = 28.0
INITIAL_C = 2.55
THRESHOLD = 0.15
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


@dataclass(frozen=True)
class Config:
    intervals: int = 200
    mesh_power: float = 1.75
    dt_s: float = 0.5
    max_time_s: float = 180.0 * 3600.0
    output_interval_s: float = 60.0
    theta: float = 0.5
    tolerance: float = 1.0e-10
    max_iterations: int = 30
    event_tolerance_s: float = 1.0e-4
    post_temperature_c: float = 50.0
    post_moisture: float = 0.05
    thermal_lock_tolerance_c: float = 1.0e-8


# 从 XLSX 工作表中读取未经表头推原始数据行。
def _xlsx_rows(path: Path) -> list[list[object]]:
    """只读解析题目附件的首张工作表，避免额外 Excel 依赖。"""
    with zipfile.ZipFile(path) as archive:
        shared: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            shared = ["".join(x.itertext()) for x in ET.fromstring(
                archive.read("xl/sharedStrings.xml"))]
        rows: list[list[object]] = []
        root = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
        for row in root.findall(".//m:sheetData/m:row", NS):
            values: dict[int, object] = {}
            for cell in row.findall("m:c", NS):
                letters = "".join(ch for ch in cell.attrib["r"] if ch.isalpha())
                column = 0
                for letter in letters:
                    column = column * 26 + ord(letter) - 64
                node = cell.find("m:v", NS)
                value: object = None if node is None else node.text
                if cell.attrib.get("t") == "s" and value is not None:
                    value = shared[int(value)]
                elif cell.attrib.get("t") == "inlineStr":
                    inline = cell.find("m:is", NS)
                    value = "" if inline is None else "".join(inline.itertext())
                elif value is not None:
                    value = float(value)
                values[column - 1] = value
            if values:
                rows.append([values.get(i) for i in range(max(values) + 1)])
    return rows


# 读取附件中的烘房与水分边界数据。
def load_environment() -> np.ndarray:
    values = np.asarray(_xlsx_rows(DATA / "附件1.xlsx")[1:], dtype=float)
    if values.shape != (241, 3):
        raise ValueError(f"附件1应有241×3条数值，实际为{values.shape}")
    if values[0, 0] != 0.0 or values[-1, 0] != 14400.0:
        raise ValueError("附件1时间范围应为0--14400 s")
    return values


# 读取药材半径随时间变的实测数据。
def load_radius_data() -> np.ndarray:
    """读取附件2，返回时间(s)-半径(m)，并执行完整输入审计。"""
    values = np.asarray(_xlsx_rows(DATA / "附件2.xlsx")[1:], dtype=float)
    if values.shape != (145, 2):
        raise ValueError(f"附件2应有145×2条数值，实际为{values.shape}")
    if values[0, 0] != 0.0 or values[-1, 0] != 259200.0:
        raise ValueError("附件2时间范围应为0--259200 s（72 h）")
    if not np.allclose(np.diff(values[:, 0]), 1800.0):
        raise ValueError("附件2时间间隔应恒为1800 s")
    values[:, 1] /= 100.0
    if not np.isclose(values[0, 1], RADIUS_M) or np.any(values[:, 1] <= 0.0):
        raise ValueError("附件2初始半径或单位异常")
    if np.any(np.diff(values[:, 1]) > 1.0e-12):
        raise ValueError("附件2半径应整体非增")
    return values


# 由半径实测数插值得到指定时刻的药材半径。
def radius_at(t_s: float, radius_data: np.ndarray, fixed_radius: bool) -> float:
    if fixed_radius:
        return RADIUS_M
    if t_s < radius_data[0, 0] - 1.0e-9 or t_s > radius_data[-1, 0] + 1.0e-9:
        raise ValueError("收缩模型超出附件2的0--72 h观测范围，禁止静默外推")
    return float(np.interp(t_s, radius_data[:, 0], radius_data[:, 1]))


# 返回指定时刻的烘房温度与环境含水状态。
def environment(t_s: float, measured: np.ndarray, cfg: Config) -> tuple[float, float]:
    if t_s <= measured[-1, 0]:
        return (float(np.interp(t_s, measured[:, 0], measured[:, 1])),
                float(np.interp(t_s, measured[:, 0], measured[:, 2])))
    return cfg.post_temperature_c, cfg.post_moisture


# 根据含水率计算单体积热容。
def volumetric_heat_capacity(c: np.ndarray) -> np.ndarray:
    rho = 760.0 + 90.0 * c
    cp = 1850.0 + 2150.0 * c / (c + 1.0)
    return rho * cp


# 根据含水率计算药材导热系数。
def conductivity(c: np.ndarray) -> np.ndarray:
    return 0.12 + 0.20 * c / (c + 1.0)


# 根据温度和含水率计水分有效扩散系数。
def diffusivity(t_c: np.ndarray, c: np.ndarray) -> np.ndarray:
    if np.any(c <= 0.0):
        raise FloatingPointError("含水率必须为正，程序不使用截断掩盖数值错误")
    return 4.2e-4 * np.exp(-0.30 / c) * np.exp(-3850.0 / (t_c + 273.15))


# 在收缩圆柱参考坐标下构造向网格及面积权重。
def reference_geometry(cfg: Config) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if cfg.intervals < 4 or cfg.mesh_power < 1.0:
        raise ValueError("intervals至少为4，mesh_power不得小于1")
    uniform = np.linspace(0.0, 1.0, cfg.intervals + 1)
    xi = 1.0 - (1.0 - uniform) ** cfg.mesh_power
    faces = np.concatenate(([0.0], 0.5 * (xi[:-1] + xi[1:]), [1.0]))
    areas = 0.5 * np.diff(faces**2)
    return xi, faces, areas


# 组装径向扩散方程的三对离散算子及边界源项。
def spatial_operator(coefficient: np.ndarray, boundary_coefficient: float,
                     boundary_value: float, radius_m: float,
                     geom: tuple[np.ndarray, ...]
                     ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    xi, faces, areas = geom
    face_values = (2.0 * coefficient[:-1] * coefficient[1:] /
                   np.maximum(coefficient[:-1] + coefficient[1:], 1.0e-300))
    conductance = faces[1:-1] * face_values / np.diff(xi) / radius_m**2
    lower = np.zeros_like(xi)
    upper = np.zeros_like(xi)
    lower[1:] = conductance / areas[1:]
    upper[:-1] = conductance / areas[:-1]
    diagonal = -(lower + upper)
    diagonal[-1] -= boundary_coefficient / radius_m / areas[-1]
    source = np.zeros_like(xi)
    source[-1] = boundary_coefficient * boundary_value / radius_m / areas[-1]
    return lower, diagonal, upper, source


# 将已组装的空间离散作用于场变量。
def apply_operator(op: tuple[np.ndarray, ...], field: np.ndarray) -> np.ndarray:
    lower, diagonal, upper, source = op
    result = diagonal * field + source
    result[1:] += lower[1:] * field[:-1]
    result[:-1] += upper[:-1] * field[1:]
    return result


# 完成一次 Crank--Nicolson 线性隐式推进。
def cn_linear_solve(old: np.ndarray, old_op: tuple[np.ndarray, ...],
                    new_op: tuple[np.ndarray, ...], capacity: np.ndarray,
                    dt_s: float, theta: float) -> tuple[np.ndarray, float]:
    lower, diagonal, upper, source = new_op
    rhs = capacity * old + dt_s * ((1.0 - theta) * apply_operator(old_op, old)
                                   + theta * source)
    matrix = np.zeros((3, old.size))
    matrix[0, 1:] = -dt_s * theta * upper[:-1]
    matrix[1] = capacity - dt_s * theta * diagonal
    matrix[2, :-1] = -dt_s * theta * lower[1:]
    answer = solve_banded((1, 1), matrix, rhs, check_finite=False)
    residual = (capacity * (answer - old) - dt_s * (
        (1.0 - theta) * apply_operator(old_op, old)
        + theta * apply_operator(new_op, answer)))
    scaled = float(np.max(np.abs(residual) / np.maximum(capacity, 1.0)))
    return answer, scaled


# 将温度场和含水率场向前推进一个时间步。
def advance(t0: np.ndarray, c0: np.ndarray, time0_s: float, dt_s: float,
            measured: np.ndarray, radius_data: np.ndarray,
            geom: tuple[np.ndarray, ...], cfg: Config,
            fixed_radius: bool = False, thermal_locked: bool = False
            ) -> tuple[np.ndarray, np.ndarray, int, float, float]:
    env0 = environment(time0_s, measured, cfg)
    env1 = environment(time0_s + dt_s, measured, cfg)
    radius_mid = radius_at(time0_s + 0.5 * dt_s, radius_data, fixed_radius)
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


# 按圆柱截面控制体权重计算场变量的面积平均值。
def area_average(field: np.ndarray, areas: np.ndarray) -> float:
    """参考截面的面积加权平均；sum(areas)=1/2。"""
    return float(np.dot(field, areas) / 0.5)


# 在一个时间步内定位平均含水率达到阈值的准确时刻。
def locate_event(t0: np.ndarray, c0: np.ndarray, time0_s: float, dt_s: float,
                 measured: np.ndarray, radius_data: np.ndarray,
                 geom: tuple[np.ndarray, ...], cfg: Config,
                 fixed_radius: bool, thermal_locked: bool
                 ) -> dict:
    low, high = 0.0, dt_s
    te = t0.copy()
    ce = c0.copy()
    while high - low > cfg.event_tolerance_s:
        middle = 0.5 * (low + high)
        tm, cm, _, _, _ = advance(
            t0, c0, time0_s, middle, measured, radius_data, geom, cfg,
            fixed_radius, thermal_locked)
        if float(np.max(cm)) < THRESHOLD:
            high, te, ce = middle, tm, cm
        else:
            low = middle
    if not np.isfinite(ce).all() or float(np.max(ce)) >= THRESHOLD:
        te, ce, _, _, _ = advance(
            t0, c0, time0_s, high, measured, radius_data, geom, cfg,
            fixed_radius, thermal_locked)
    return {"time_s": time0_s + high, "bracket_s": [time0_s + low, time0_s + high],
            "temperature_C": te, "moisture_kg_kg": ce,
            "max_C": float(np.max(ce)), "argmax_node": int(np.argmax(ce)),
            "average_C": area_average(ce, geom[2]), "surface_C": float(ce[-1]),
            "radius_m": radius_at(time0_s + high, radius_data, fixed_radius)}


# 按给定配置执行完整时程的热湿耦合数值模拟。
def simulate(cfg: Config, stop_at_threshold: bool = True, verbose: bool = True,
             fixed_radius: bool = False) -> dict:
    if cfg.dt_s <= 0 or cfg.output_interval_s <= 0:
        raise ValueError("时间步和输出间隔必须为正")
    if cfg.theta != 0.5:
        raise ValueError("第四问正式算法固定使用Crank--Nicolson(theta=0.5)")
    measured = load_environment()
    radius_data = load_radius_data()
    if not fixed_radius and cfg.max_time_s > radius_data[-1, 0]:
        # 收缩模型不得超过附件2；固定半径控制模型不受此限制。
        max_time_s = float(radius_data[-1, 0])
    else:
        max_time_s = cfg.max_time_s
    geom = reference_geometry(cfg)
    xi, _, areas = geom
    t = np.full(xi.size, INITIAL_T_C)
    c = np.full(xi.size, INITIAL_C)
    now = 0.0
    next_output = cfg.output_interval_s
    times = [0.0]
    stored_radius = [radius_at(0.0, radius_data, fixed_radius)]
    stored_t = [t.copy()]
    stored_c = [c.copy()]
    iterations: list[int] = []
    max_residual_t = 0.0
    max_residual_c = 0.0
    max_mass_balance = 0.0
    cumulative_outflow = 0.0
    monotonic_violations = 0
    first_violation: dict | None = None
    max_radial_increase = 0.0
    max_radial_increase_time_s = None
    max_center_gap = 0.0
    negative_count = 0
    nonfinite_count = 0
    event = None
    mean_event = None
    thermal_locked = False
    thermal_lock_time_s = None
    thermal_lock_deviation_c = None
    previous_average = area_average(c, areas)
    started = time.perf_counter()
    progress_mark = 6.0 * 3600.0
    while now < max_time_s - 1.0e-12:
        step = min(cfg.dt_s, max_time_s - now)
        if now < 14400.0 < now + step:
            step = 14400.0 - now
        if now < next_output < now + step:
            step = next_output - now
        tn, cn, count, residual_t, residual_c = advance(
            t, c, now, step, measured, radius_data, geom, cfg,
            fixed_radius, thermal_locked)
        iterations.append(count)
        max_residual_t = max(max_residual_t, residual_t)
        max_residual_c = max(max_residual_c, residual_c)
        nonfinite_count += int(not np.isfinite(tn).all() or not np.isfinite(cn).all())
        negative_count += int(float(np.min(cn)) < 0.0)
        outward_increment = np.diff(cn)
        violation = float(np.max(outward_increment))
        if violation > max_radial_increase:
            max_radial_increase = violation
            max_radial_increase_time_s = now + step
        if violation > 1.0e-9:
            monotonic_violations += 1
            if first_violation is None:
                idx = int(np.argmax(outward_increment))
                first_violation = {"time_s": now + step, "left_node": idx,
                                   "right_node": idx + 1, "increase": violation}
        max_center_gap = max(max_center_gap, float(np.max(cn) - cn[0]))
        env0 = environment(now, measured, cfg)[1]
        env1 = environment(now + step, measured, cfg)[1]
        radius_mid = radius_at(now + 0.5 * step, radius_data, fixed_radius)
        outflow = step * MASS_TRANSFER / radius_mid * (
            c[-1] - env0 + cn[-1] - env1)
        mass_error = area_average(cn - c, areas) + outflow
        max_mass_balance = max(max_mass_balance, abs(mass_error))
        cumulative_outflow += outflow
        current_average = area_average(cn, areas)
        if mean_event is None and previous_average >= THRESHOLD > current_average:
            fraction = ((previous_average - THRESHOLD) /
                        (previous_average - current_average))
            mean_event = now + fraction * step
        if event is None and float(np.max(c)) >= THRESHOLD > float(np.max(cn)):
            event = locate_event(
                t, c, now, step, measured, radius_data, geom, cfg,
                fixed_radius, thermal_locked)
        if (not thermal_locked and cfg.thermal_lock_tolerance_c > 0.0
                and now + step >= 14400.0):
            deviation = float(np.max(np.abs(tn - cfg.post_temperature_c)))
            if deviation <= cfg.thermal_lock_tolerance_c:
                thermal_locked = True
                thermal_lock_time_s = now + step
                thermal_lock_deviation_c = deviation
        t, c, now = tn, cn, now + step
        previous_average = current_average
        if math.isclose(now, next_output, abs_tol=1.0e-9):
            times.append(now)
            stored_radius.append(radius_at(now, radius_data, fixed_radius))
            stored_t.append(t.copy())
            stored_c.append(c.copy())
            next_output += cfg.output_interval_s
        if verbose and now >= progress_mark:
            print(f"N={cfg.intervals}, dt={cfg.dt_s:g}s, t={now/3600:.1f}h, "
                  f"Cmax={np.max(c):.6f}", flush=True)
            progress_mark += 6.0 * 3600.0
        if event is not None and stop_at_threshold:
            break
    if stop_at_threshold and event is None:
        suffix = "（附件2覆盖上限）" if not fixed_radius else ""
        raise RuntimeError(f"{max_time_s/3600:g} h内未达到阈值{suffix}")
    diagnostics = {
        "steps": len(iterations),
        "runtime_s": time.perf_counter() - started,
        "mean_picard_iterations": float(np.mean(iterations)),
        "max_picard_iterations": int(np.max(iterations)),
        "unconverged_steps": 0,
        "max_temperature_algebraic_residual": max_residual_t,
        "max_moisture_algebraic_residual": max_residual_c,
        "max_step_mass_balance_abs": max_mass_balance,
        "cumulative_reference_balance_abs": abs(
            area_average(c, areas) - INITIAL_C + cumulative_outflow),
        "radial_monotonicity_violations": monotonic_violations,
        "first_radial_violation": first_violation,
        "max_outward_radial_increase": max_radial_increase,
        "max_outward_radial_increase_time_s": max_radial_increase_time_s,
        "max_center_vs_global_gap": max_center_gap,
        "negative_moisture_steps": negative_count,
        "nonfinite_steps": nonfinite_count,
        "minimum_moisture": float(np.min(stored_c)),
        "maximum_temperature_C": float(np.max(stored_t)),
        "minimum_dxi": float(np.min(np.diff(xi))),
        "radius_positive": bool(np.all(radius_data[:, 1] > 0.0)),
        "radius_nonincreasing": bool(np.all(np.diff(radius_data[:, 1]) <= 1.0e-12)),
        "radius_data_end_s": float(radius_data[-1, 0]),
        "fixed_radius_control": fixed_radius,
        "thermal_lock_enabled": bool(cfg.thermal_lock_tolerance_c > 0.0),
        "thermal_lock_time_s": thermal_lock_time_s,
        "thermal_lock_deviation_c": thermal_lock_deviation_c,
    }
    return {"config": asdict(cfg), "xi": xi, "area_weights": areas,
            "radius_time_m": np.asarray(stored_radius),
            "time_s": np.asarray(times), "temperature_C": np.asarray(stored_t),
            "moisture_kg_kg": np.asarray(stored_c), "event": event,
            "mean_threshold_time_s": mean_event, "diagnostics": diagnostics,
            "fixed_radius_control": fixed_radius}


# 提取可序列化的关键结果，形成计算摘要。
def serializable_summary(result: dict) -> dict:
    event = result["event"]
    event_summary = None if event is None else {
        key: value for key, value in event.items()
        if key not in ("temperature_C", "moisture_kg_kg")}
    return {"config": result["config"], "event": event_summary,
            "fixed_radius_control": result["fixed_radius_control"],
            "mean_threshold_time_s": result["mean_threshold_time_s"],
            "diagnostics": result["diagnostics"]}


# 保存指定工况的完整数组结果和摘要文件。
def save_run(result: dict, label: str) -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    summary_path = RESULTS / f"{label}_summary.json"
    summary_path.write_text(json.dumps(serializable_summary(result), ensure_ascii=False, indent=2),
                            encoding="utf-8")
    event = result["event"]
    payload = {key: result[key] for key in (
        "xi", "area_weights", "radius_time_m", "time_s",
        "temperature_C", "moisture_kg_kg")}
    if event is not None:
        payload.update(event_time_s=event["time_s"],
                       event_radius_m=event["radius_m"],
                       event_temperature_C=event["temperature_C"],
                       event_moisture_kg_kg=event["moisture_kg_kg"])
    np.savez_compressed(RESULTS / f"{label}_fields.npz", **payload)


# 组织当前脚本的完整执行流程并返回运行状态。
def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--intervals", type=int, default=200)
    parser.add_argument("--mesh-power", type=float, default=1.75)
    parser.add_argument("--dt", type=float, default=0.5)
    parser.add_argument("--max-hours", type=float, default=180.0)
    parser.add_argument("--output-seconds", type=float, default=60.0)
    parser.add_argument("--label", default="q4_run")
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    parser.add_argument("--no-thermal-lock", action="store_true")
    parser.add_argument("--fixed-radius", action="store_true",
                        help="使用附录4物性但令R恒为2 cm，作为收缩效应控制模型")
    args = parser.parse_args()
    cfg = Config(intervals=20 if args.smoke else args.intervals,
                 mesh_power=args.mesh_power,
                 dt_s=2.0 if args.smoke else args.dt,
                 max_time_s=600.0 if args.smoke else args.max_hours * 3600.0,
                 output_interval_s=60.0 if args.smoke else args.output_seconds,
                 thermal_lock_tolerance_c=0.0 if args.no_thermal_lock else 1.0e-8)
    result = simulate(cfg, stop_at_threshold=not args.smoke,
                      verbose=not args.quiet, fixed_radius=args.fixed_radius)
    save_run(result, "smoke" if args.smoke else args.label)
    print(json.dumps(serializable_summary(result), ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
