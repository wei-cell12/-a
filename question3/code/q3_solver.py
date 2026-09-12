#!/usr/bin/env python3
"""问题三：固定半径圆柱的非线性热湿耦合求解器。

空间采用节点中心有限体积法，时间采用 Crank--Nicolson 格式，
每一步使用 Picard 迭代同步更新物性。停止条件检查全部径向节点。
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


ROOT = Path(__file__).resolve().parent.parent
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
    mesh_power: float = 2.0
    dt_s: float = 0.5
    max_time_s: float = 80.0 * 3600.0
    output_interval_s: float = 60.0
    theta: float = 0.5
    tolerance: float = 1.0e-10
    max_iterations: int = 30
    event_tolerance_s: float = 1.0e-4
    post_temperature_c: float = 50.0
    post_moisture: float = 0.05
    thermal_lock_tolerance_c: float = 1.0e-8


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


def load_environment() -> np.ndarray:
    values = np.asarray(_xlsx_rows(DATA / "附件1.xlsx")[1:], dtype=float)
    if values.shape != (241, 3):
        raise ValueError(f"附件1应有241×3条数值，实际为{values.shape}")
    if values[0, 0] != 0.0 or values[-1, 0] != 14400.0:
        raise ValueError("附件1时间范围应为0--14400 s")
    return values


def environment(t_s: float, measured: np.ndarray, cfg: Config) -> tuple[float, float]:
    if t_s <= measured[-1, 0]:
        return (float(np.interp(t_s, measured[:, 0], measured[:, 1])),
                float(np.interp(t_s, measured[:, 0], measured[:, 2])))
    return cfg.post_temperature_c, cfg.post_moisture


def volumetric_heat_capacity(c: np.ndarray) -> np.ndarray:
    rho = 650.0 + 128.0 * c
    cp = 1450.0 + 2736.0 * c / (c + 1.0)
    return rho * cp


def conductivity(c: np.ndarray) -> np.ndarray:
    return 0.21 + 0.38 * c / (c + 1.0)


def diffusivity(t_c: np.ndarray, c: np.ndarray) -> np.ndarray:
    if np.any(c <= 0.0):
        raise FloatingPointError("含水率必须为正，程序不使用截断掩盖数值错误")
    return 2.4e-3 * np.exp(-0.45 / c) * np.exp(-3850.0 / (t_c + 273.15))


def radial_geometry(cfg: Config) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if cfg.intervals < 4 or cfg.mesh_power < 1.0:
        raise ValueError("intervals至少为4，mesh_power不得小于1")
    xi = np.linspace(0.0, 1.0, cfg.intervals + 1)
    r = RADIUS_M * (1.0 - (1.0 - xi) ** cfg.mesh_power)
    faces = np.concatenate(([0.0], 0.5 * (r[:-1] + r[1:]), [RADIUS_M]))
    volumes = 0.5 * np.diff(faces**2)
    return r, faces, volumes


def spatial_operator(coefficient: np.ndarray, boundary_coefficient: float,
                     boundary_value: float, geom: tuple[np.ndarray, ...]
                     ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    r, faces, volumes = geom
    face_values = (2.0 * coefficient[:-1] * coefficient[1:] /
                   np.maximum(coefficient[:-1] + coefficient[1:], 1.0e-300))
    conductance = faces[1:-1] * face_values / np.diff(r)
    lower = np.zeros_like(r)
    upper = np.zeros_like(r)
    lower[1:] = conductance / volumes[1:]
    upper[:-1] = conductance / volumes[:-1]
    diagonal = -(lower + upper)
    diagonal[-1] -= RADIUS_M * boundary_coefficient / volumes[-1]
    source = np.zeros_like(r)
    source[-1] = RADIUS_M * boundary_coefficient * boundary_value / volumes[-1]
    return lower, diagonal, upper, source


def apply_operator(op: tuple[np.ndarray, ...], field: np.ndarray) -> np.ndarray:
    lower, diagonal, upper, source = op
    result = diagonal * field + source
    result[1:] += lower[1:] * field[:-1]
    result[:-1] += upper[:-1] * field[1:]
    return result


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


def advance(t0: np.ndarray, c0: np.ndarray, time0_s: float, dt_s: float,
            measured: np.ndarray, geom: tuple[np.ndarray, ...], cfg: Config,
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


def area_average(field: np.ndarray, volumes: np.ndarray) -> float:
    return float(np.dot(field, volumes) / (0.5 * RADIUS_M**2))


def locate_event(t0: np.ndarray, c0: np.ndarray, time0_s: float, dt_s: float,
                 measured: np.ndarray, geom: tuple[np.ndarray, ...], cfg: Config,
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


def simulate(cfg: Config, stop_at_threshold: bool = True, verbose: bool = True) -> dict:
    if cfg.dt_s <= 0 or cfg.output_interval_s <= 0:
        raise ValueError("时间步和输出间隔必须为正")
    if cfg.theta != 0.5:
        raise ValueError("第三问正式算法固定使用Crank--Nicolson(theta=0.5)")
    measured = load_environment()
    geom = radial_geometry(cfg)
    r, _, volumes = geom
    t = np.full(r.size, INITIAL_T_C)
    c = np.full(r.size, INITIAL_C)
    now = 0.0
    next_output = cfg.output_interval_s
    times = [0.0]
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
    previous_average = area_average(c, volumes)
    started = time.perf_counter()
    progress_mark = 6.0 * 3600.0
    while now < cfg.max_time_s - 1.0e-12:
        step = min(cfg.dt_s, cfg.max_time_s - now)
        if now < 14400.0 < now + step:
            step = 14400.0 - now
        if now < next_output < now + step:
            step = next_output - now
        tn, cn, count, residual_t, residual_c = advance(
            t, c, now, step, measured, geom, cfg, thermal_locked)
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
        outflow = step * MASS_TRANSFER / RADIUS_M * (
            c[-1] - env0 + cn[-1] - env1)
        mass_error = area_average(cn - c, volumes) + outflow
        max_mass_balance = max(max_mass_balance, abs(mass_error))
        cumulative_outflow += outflow
        current_average = area_average(cn, volumes)
        if mean_event is None and previous_average >= THRESHOLD > current_average:
            fraction = ((previous_average - THRESHOLD) /
                        (previous_average - current_average))
            mean_event = now + fraction * step
        if event is None and float(np.max(c)) >= THRESHOLD > float(np.max(cn)):
            event = locate_event(t, c, now, step, measured, geom, cfg, thermal_locked)
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
        raise RuntimeError(f"{cfg.max_time_s/3600:g} h内未达到阈值")
    diagnostics = {
        "steps": len(iterations),
        "runtime_s": time.perf_counter() - started,
        "mean_picard_iterations": float(np.mean(iterations)),
        "max_picard_iterations": int(np.max(iterations)),
        "unconverged_steps": 0,
        "max_temperature_algebraic_residual": max_residual_t,
        "max_moisture_algebraic_residual": max_residual_c,
        "max_step_mass_balance_abs": max_mass_balance,
        "cumulative_mass_balance_abs": abs(area_average(c, volumes) - INITIAL_C + cumulative_outflow),
        "radial_monotonicity_violations": monotonic_violations,
        "first_radial_violation": first_violation,
        "max_outward_radial_increase": max_radial_increase,
        "max_outward_radial_increase_time_s": max_radial_increase_time_s,
        "max_center_vs_global_gap": max_center_gap,
        "negative_moisture_steps": negative_count,
        "nonfinite_steps": nonfinite_count,
        "minimum_moisture": float(np.min(stored_c)),
        "maximum_temperature_C": float(np.max(stored_t)),
        "minimum_dr_m": float(np.min(np.diff(r))),
        "thermal_lock_enabled": bool(cfg.thermal_lock_tolerance_c > 0.0),
        "thermal_lock_time_s": thermal_lock_time_s,
        "thermal_lock_deviation_c": thermal_lock_deviation_c,
    }
    return {"config": asdict(cfg), "radius_m": r, "volume_m2": volumes,
            "time_s": np.asarray(times), "temperature_C": np.asarray(stored_t),
            "moisture_kg_kg": np.asarray(stored_c), "event": event,
            "mean_threshold_time_s": mean_event, "diagnostics": diagnostics}


def serializable_summary(result: dict) -> dict:
    event = result["event"]
    event_summary = None if event is None else {
        key: value for key, value in event.items()
        if key not in ("temperature_C", "moisture_kg_kg")}
    return {"config": result["config"], "event": event_summary,
            "mean_threshold_time_s": result["mean_threshold_time_s"],
            "diagnostics": result["diagnostics"]}


def save_run(result: dict, label: str) -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    summary_path = RESULTS / f"{label}_summary.json"
    summary_path.write_text(json.dumps(serializable_summary(result), ensure_ascii=False, indent=2),
                            encoding="utf-8")
    event = result["event"]
    payload = {key: result[key] for key in (
        "radius_m", "volume_m2", "time_s", "temperature_C", "moisture_kg_kg")}
    if event is not None:
        payload.update(event_time_s=event["time_s"],
                       event_temperature_C=event["temperature_C"],
                       event_moisture_kg_kg=event["moisture_kg_kg"])
    np.savez_compressed(RESULTS / f"{label}_fields.npz", **payload)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--intervals", type=int, default=200)
    parser.add_argument("--mesh-power", type=float, default=2.0)
    parser.add_argument("--dt", type=float, default=0.5)
    parser.add_argument("--max-hours", type=float, default=80.0)
    parser.add_argument("--output-seconds", type=float, default=60.0)
    parser.add_argument("--label", default="q3_run")
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    parser.add_argument("--no-thermal-lock", action="store_true")
    args = parser.parse_args()
    cfg = Config(intervals=20 if args.smoke else args.intervals,
                 mesh_power=args.mesh_power,
                 dt_s=2.0 if args.smoke else args.dt,
                 max_time_s=600.0 if args.smoke else args.max_hours * 3600.0,
                 output_interval_s=60.0 if args.smoke else args.output_seconds,
                 thermal_lock_tolerance_c=0.0 if args.no_thermal_lock else 1.0e-8)
    result = simulate(cfg, stop_at_threshold=not args.smoke, verbose=not args.quiet)
    save_run(result, "smoke" if args.smoke else args.label)
    print(json.dumps(serializable_summary(result), ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
