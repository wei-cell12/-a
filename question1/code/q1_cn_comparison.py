#!/usr/bin/env python3
"""问题一：Crank–Nicolson 主求解 + 后向欧拉基线与数值验证。"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.interpolate import Akima1DInterpolator, CubicSpline, PchipInterpolator
from scipy.linalg import solve_banded
from scipy.optimize import brentq
from scipy.special import j0, j1, jn_zeros

import 问题1_求解 as legacy
from utils.plot_style import COLOR_SEQUENCE, PALETTE, apply_publication_style


ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
R = legacy.RADIUS
RHO, CP, K = legacy.RHO, legacy.CP, legacy.K
H, HM = legacy.H, legacy.HM
T0, C0 = legacy.T0, legacy.C0
KEY_TIMES = legacy.KEY_TIMES
KEY_RADII_CM = legacy.KEY_RADII_CM
OUTPUT_RADII_CM = legacy.OUTPUT_RADII_CM
PICARD_TOL = 1.0e-11


@dataclass(frozen=True)
class RunConfig:
    method: str
    intervals: int
    dt: float
    interpolation: str = "linear"
    face_mean: str = "harmonic"
    end_time: float = 1800.0
    store_interval: float = 1.0

    @property
    def theta(self) -> float:
        if self.method == "be":
            return 1.0
        if self.method == "cn":
            return 0.5
        raise ValueError(f"未知时间推进方法: {self.method}")


def make_interpolator(x: np.ndarray, y: np.ndarray, method: str) -> Callable[[float | np.ndarray], np.ndarray]:
    if method == "linear":
        return lambda z: np.interp(z, x, y)
    if method == "pchip":
        return PchipInterpolator(x, y, extrapolate=True)
    if method == "cubic":
        return CubicSpline(x, y, bc_type="natural", extrapolate=True)
    if method == "akima":
        return Akima1DInterpolator(x, y)
    raise ValueError(f"未知插值方法: {method}")


def interpolation_diagnostics(x: np.ndarray, y: np.ndarray) -> list[dict[str, float | str]]:
    """内部留一验证；端点不外推，并检查区间过冲。"""
    rows: list[dict[str, float | str]] = []
    for method in ("linear", "pchip", "cubic", "akima"):
        errors = []
        for i in range(1, len(x) - 1):
            keep = np.ones(len(x), dtype=bool)
            keep[i] = False
            predictor = make_interpolator(x[keep], y[keep], method)
            errors.append(float(predictor(x[i])) - float(y[i]))
        errors_array = np.asarray(errors)
        predictor = make_interpolator(x, y, method)
        overshoot = 0.0
        for i in range(len(x) - 1):
            dense = np.linspace(x[i], x[i + 1], 21)
            values = np.asarray(predictor(dense), dtype=float)
            low, high = sorted((float(y[i]), float(y[i + 1])))
            overshoot = max(overshoot, float(np.max(values) - high), float(low - np.min(values)))
        rows.append({
            "method": method,
            "mae": float(np.mean(np.abs(errors_array))),
            "rmse": float(np.sqrt(np.mean(errors_array ** 2))),
            "max_abs_error": float(np.max(np.abs(errors_array))),
            "max_interval_overshoot": max(0.0, overshoot),
        })
    return rows


def radial_geometry(intervals: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if intervals < 2:
        raise ValueError("径向区间数必须至少为2")
    radii = np.linspace(0.0, R, intervals + 1)
    dr = R / intervals
    west_faces = np.maximum(radii - 0.5 * dr, 0.0)
    east_faces = np.minimum(radii + 0.5 * dr, R)
    volumes = 0.5 * (east_faces**2 - west_faces**2)
    return radii, volumes, east_faces


def face_values(values: np.ndarray, mean: str) -> np.ndarray:
    if mean == "harmonic":
        return 2.0 * values[:-1] * values[1:] / np.maximum(values[:-1] + values[1:], 1e-300)
    if mean == "arithmetic":
        return 0.5 * (values[:-1] + values[1:])
    raise ValueError(f"未知界面平均方法: {mean}")


def spatial_operator(
    coefficient: np.ndarray,
    capacity: float,
    boundary_coefficient: float,
    boundary_value: float,
    radii: np.ndarray,
    volumes: np.ndarray,
    east_faces: np.ndarray,
    mean: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """返回 du/dt=L u+s 的三对角系数和源项。"""
    n = radii.size
    dr = radii[1] - radii[0]
    conductance = east_faces[:-1] * face_values(coefficient, mean) / dr
    west = np.zeros(n)
    east = np.zeros(n)
    west[1:] = conductance
    east[:-1] = conductance
    boundary = np.zeros(n)
    boundary[-1] = R * boundary_coefficient
    scale = 1.0 / (capacity * volumes)
    lower = scale * west
    diagonal = -scale * (west + east + boundary)
    upper = scale * east
    source = scale * boundary * boundary_value
    return lower, diagonal, upper, source


def apply_tridiagonal(lower: np.ndarray, diagonal: np.ndarray, upper: np.ndarray, x: np.ndarray) -> np.ndarray:
    result = diagonal * x
    result[1:] += lower[1:] * x[:-1]
    result[:-1] += upper[:-1] * x[1:]
    return result


def theta_step(
    old: np.ndarray,
    old_op: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    new_op: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    dt: float,
    theta: float,
) -> tuple[np.ndarray, float]:
    lo0, di0, up0, src0 = old_op
    lo1, di1, up1, src1 = new_op
    rhs = old + dt * (1.0 - theta) * (apply_tridiagonal(lo0, di0, up0, old) + src0)
    rhs += dt * theta * src1
    banded = np.zeros((3, old.size))
    banded[0, 1:] = -dt * theta * up1[:-1]
    banded[1] = 1.0 - dt * theta * di1
    banded[2, :-1] = -dt * theta * lo1[1:]
    solved = solve_banded((1, 1), banded, rhs, check_finite=False)
    residual = solved - old - dt * (
        (1.0 - theta) * (apply_tridiagonal(lo0, di0, up0, old) + src0)
        + theta * (apply_tridiagonal(lo1, di1, up1, solved) + src1)
    )
    return solved, float(np.max(np.abs(residual)))


def area_average(field: np.ndarray, volumes: np.ndarray) -> float:
    return float(np.dot(field, volumes) / (0.5 * R**2))


def diffusivity(c: np.ndarray) -> np.ndarray:
    return 7.0e-9 * np.exp(-0.89 / np.maximum(c, 1e-12))


def simulate(
    env_time: np.ndarray,
    env_temp: np.ndarray,
    env_moisture: np.ndarray,
    config: RunConfig,
) -> dict:
    steps = int(round(config.end_time / config.dt))
    stride = int(round(config.store_interval / config.dt))
    if not math.isclose(steps * config.dt, config.end_time, abs_tol=1e-10):
        raise ValueError("end_time 必须能被 dt 整除")
    if not math.isclose(stride * config.dt, config.store_interval, abs_tol=1e-10):
        raise ValueError("输出间隔必须能被 dt 整除")
    theta = config.theta
    t_boundary = make_interpolator(env_time, env_temp, config.interpolation)
    c_boundary = make_interpolator(env_time, env_moisture, config.interpolation)
    radii, volumes, east_faces = radial_geometry(config.intervals)
    n = radii.size
    temperature = np.full(n, T0)
    moisture = np.full(n, C0)
    out_times = np.arange(0.0, config.end_time + 0.5 * config.store_interval, config.store_interval)
    out_t = np.empty((out_times.size, n))
    out_c = np.empty_like(out_t)
    out_t[0], out_c[0] = temperature, moisture
    picard_counts = np.zeros(steps, dtype=np.int16)
    heat_residuals = np.empty(steps)
    mass_residuals = np.empty(steps)
    linear_residual_max = 0.0
    heat_flux_integral = 0.0
    mass_flux_integral = 0.0
    coefficient_t = np.full(n, K)
    store_index = 1

    for step in range(1, steps + 1):
        time_old = (step - 1) * config.dt
        time_new = step * config.dt
        tb0, tb1 = float(t_boundary(time_old)), float(t_boundary(time_new))
        cb0, cb1 = float(c_boundary(time_old)), float(c_boundary(time_new))
        previous_t_average = area_average(temperature, volumes)
        previous_c_average = area_average(moisture, volumes)

        top0 = spatial_operator(coefficient_t, RHO * CP, H, tb0, radii, volumes, east_faces, config.face_mean)
        top1 = spatial_operator(coefficient_t, RHO * CP, H, tb1, radii, volumes, east_faces, config.face_mean)
        temperature_new, residual = theta_step(temperature, top0, top1, config.dt, theta)
        linear_residual_max = max(linear_residual_max, residual)

        mop0 = spatial_operator(diffusivity(moisture), 1.0, HM, cb0, radii, volumes, east_faces, config.face_mean)
        guess = moisture.copy()
        moisture_new = None
        for iteration in range(1, 51):
            mop1 = spatial_operator(diffusivity(guess), 1.0, HM, cb1, radii, volumes, east_faces, config.face_mean)
            updated, residual = theta_step(moisture, mop0, mop1, config.dt, theta)
            linear_residual_max = max(linear_residual_max, residual)
            if np.max(np.abs(updated - guess)) < PICARD_TOL:
                moisture_new = updated
                picard_counts[step - 1] = iteration
                break
            guess = updated
        if moisture_new is None:
            raise RuntimeError(f"Picard 迭代在 t={time_new:g} s 未收敛")

        heat_increment = config.dt * 2.0 * H / (RHO * CP * R) * (
            (1.0 - theta) * (tb0 - temperature[-1]) + theta * (tb1 - temperature_new[-1])
        )
        mass_increment = config.dt * 2.0 * HM / R * (
            (1.0 - theta) * (cb0 - moisture[-1]) + theta * (cb1 - moisture_new[-1])
        )
        heat_residuals[step - 1] = area_average(temperature_new, volumes) - previous_t_average - heat_increment
        mass_residuals[step - 1] = area_average(moisture_new, volumes) - previous_c_average - mass_increment
        heat_flux_integral += heat_increment
        mass_flux_integral += mass_increment
        temperature, moisture = temperature_new, moisture_new

        if step % stride == 0:
            out_t[store_index], out_c[store_index] = temperature, moisture
            store_index += 1

    return {
        "config": config,
        "time": out_times,
        "radius_m": radii,
        "temperature": out_t,
        "moisture": out_c,
        "picard_counts": picard_counts,
        "max_picard_iterations": int(np.max(picard_counts)),
        "mean_picard_iterations": float(np.mean(picard_counts)),
        "linear_system_residual_max": linear_residual_max,
        "temperature_step_conservation_max": float(np.max(np.abs(heat_residuals))),
        "moisture_step_conservation_max": float(np.max(np.abs(mass_residuals))),
        "temperature_cumulative_residual": area_average(temperature, volumes) - T0 - heat_flux_integral,
        "moisture_cumulative_residual": area_average(moisture, volumes) - C0 - mass_flux_integral,
        "temperature_conservation_history": heat_residuals,
        "moisture_conservation_history": mass_residuals,
        "temperature_physical_lower": min(T0, float(np.min(t_boundary(out_times)))),
        "temperature_physical_upper": max(T0, float(np.max(t_boundary(out_times)))),
        "moisture_physical_upper": max(C0, float(np.max(c_boundary(out_times)))),
    }


def extract(result: dict, times: np.ndarray = KEY_TIMES, radii_cm: np.ndarray = KEY_RADII_CM, key: str = "temperature") -> np.ndarray:
    source_times = np.asarray(result["time"])
    source_radii = np.asarray(result["radius_m"]) * 100.0
    field = np.asarray(result[key])
    values = []
    for target in times:
        index = int(np.argmin(np.abs(source_times - target)))
        if not math.isclose(float(source_times[index]), float(target), abs_tol=1e-9):
            raise ValueError(f"未存储目标时刻 {target}")
        values.append(np.interp(radii_cm, source_radii, field[index]))
    return np.asarray(values)


def difference_metrics(a: dict, b: dict) -> dict[str, float]:
    dt = extract(a, key="temperature") - extract(b, key="temperature")
    dc = extract(a, key="moisture") - extract(b, key="moisture")
    return {
        "temperature_max_abs": float(np.max(np.abs(dt))),
        "temperature_rmse": float(np.sqrt(np.mean(dt**2))),
        "moisture_max_abs": float(np.max(np.abs(dc))),
        "moisture_rmse": float(np.sqrt(np.mean(dc**2))),
    }


def validate_result(result: dict) -> dict[str, float | int | bool]:
    t = np.asarray(result["temperature"])
    c = np.asarray(result["moisture"])
    finite = bool(np.all(np.isfinite(t)) and np.all(np.isfinite(c)))
    monotone_t = bool(np.all(np.diff(t, axis=1) >= -2e-8))
    monotone_c = bool(np.all(np.diff(c, axis=1) <= 2e-8))
    passed = finite and np.min(t) >= float(result["temperature_physical_lower"]) - 1e-8
    passed = passed and np.max(t) <= float(result["temperature_physical_upper"]) + 1e-8
    passed = passed and np.min(c) >= -1e-10 and np.max(c) <= float(result["moisture_physical_upper"]) + 1e-8
    passed = passed and monotone_t and monotone_c
    if not passed:
        raise AssertionError("物理范围、有限性或径向单调性检查失败")
    return {
        "passed": passed,
        "temperature_min": float(np.min(t)),
        "temperature_max": float(np.max(t)),
        "moisture_min": float(np.min(c)),
        "moisture_max": float(np.max(c)),
        "radially_monotone_temperature": monotone_t,
        "radially_monotone_moisture": monotone_c,
        "max_picard_iterations": result["max_picard_iterations"],
        "mean_picard_iterations": result["mean_picard_iterations"],
        "linear_system_residual_max": result["linear_system_residual_max"],
        "temperature_cumulative_residual": result["temperature_cumulative_residual"],
        "moisture_cumulative_residual": result["moisture_cumulative_residual"],
    }


def bessel_roots(bi: float, count: int = 80) -> np.ndarray:
    zeros0 = jn_zeros(0, count + 1)
    roots = []
    f = lambda z: z * j1(z) - bi * j0(z)
    left = 1e-12
    for right in zeros0:
        if f(left) * f(right) < 0:
            roots.append(brentq(f, left, right))
        left = right + 1e-10
        if len(roots) == count:
            break
    return np.asarray(roots)


def analytic_cylinder(r: np.ndarray, t: float, diffusivity_value: float, bi: float, initial: float, ambient: float) -> np.ndarray:
    roots = bessel_roots(bi)
    coeff = 2.0 * j1(roots) / (roots * (j0(roots) ** 2 + j1(roots) ** 2))
    ratio = np.sum(
        coeff[:, None] * j0(roots[:, None] * r[None, :] / R)
        * np.exp(-(roots[:, None] ** 2) * diffusivity_value * t / R**2), axis=0
    )
    return ambient + (initial - ambient) * ratio


def simulate_constant_linear(coefficient: float, capacity: float, boundary_coefficient: float,
                             initial: float, ambient: float, intervals: int, dt: float) -> dict:
    """供解析解核验使用的常系数 CN 求解器。"""
    radii, volumes, east_faces = radial_geometry(intervals)
    state = np.full(intervals + 1, initial)
    coefficient_field = np.full(intervals + 1, coefficient)
    operator = spatial_operator(coefficient_field, capacity, boundary_coefficient, ambient,
                                radii, volumes, east_faces, "harmonic")
    times = np.arange(0.0, 1800.0 + 300.0, 300.0)
    stored = [state.copy()]
    target_index = 1
    for step in range(1, int(round(1800.0 / dt)) + 1):
        state, _ = theta_step(state, operator, operator, dt, 0.5)
        if math.isclose(step * dt, times[target_index], abs_tol=1e-10):
            stored.append(state.copy())
            target_index += 1
            if target_index == len(times):
                break
    return {"time": times, "radius_m": radii, "field": np.asarray(stored)}


def constant_boundary_validation(intervals: int = 200, dt: float = 0.5) -> dict[str, float]:
    thermal = simulate_constant_linear(K, RHO * CP, H, T0, 50.0, intervals, dt)
    d0 = float(diffusivity(np.array([C0]))[0])
    moisture = simulate_constant_linear(d0, 1.0, HM, C0, 0.05, intervals, dt)
    radii = np.asarray(thermal["radius_m"])
    temperature_errors, moisture_errors = [], []
    for row, target_time in enumerate(thermal["time"][1:], start=1):
        exact_t = analytic_cylinder(radii, float(target_time), K / (RHO * CP), H * R / K, T0, 50.0)
        exact_c = analytic_cylinder(radii, float(target_time), d0, HM * R / d0, C0, 0.05)
        temperature_errors.append(np.max(np.abs(np.asarray(thermal["field"])[row] - exact_t)))
        moisture_errors.append(np.max(np.abs(np.asarray(moisture["field"])[row] - exact_c)))
    return {
        "temperature_bessel_max_abs_degC": float(np.max(temperature_errors)),
        "moisture_bessel_max_abs_kg_per_kg": float(np.max(moisture_errors)),
    }


def write_rows(path: Path, header: list[str], rows: list[list[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def export_plot(fig, name: str, size: tuple[float, float] = (6.3, 3.9)) -> None:
    fig.set_size_inches(*size)
    svg_path = FIGURES / f"{name}.svg"
    png_path = FIGURES / f"{name}.png"
    fig.savefig(svg_path)
    fig.savefig(png_path, dpi=300)
    from PIL import Image
    qa = FIGURES / "_qa"
    qa.mkdir(exist_ok=True)
    Image.open(png_path).convert("L").save(qa / f"{name}_grayscale.png", dpi=(300, 300))
    plt.close(fig)


def make_convergence_figure(convergence_rows: list[list[object]]) -> None:
    grid_rows = [r for r in convergence_rows if r[0] == "grid"][:-1]
    time_rows = [r for r in convergence_rows if r[0] == "time"][:-1]
    fig, axes = plt.subplots(2, 2, layout="constrained")
    panels = (
        (axes[0, 0], grid_rows, 1, 3, "Radial intervals N", "Temperature difference (degC)", "Temperature: grid"),
        (axes[1, 0], grid_rows, 1, 4, "Radial intervals N", "Moisture difference (kg/kg)", "Moisture: grid"),
        (axes[0, 1], time_rows, 2, 3, "Time step (s)", "Temperature difference (degC)", "Temperature: time step"),
        (axes[1, 1], time_rows, 2, 4, "Time step (s)", "Moisture difference (kg/kg)", "Moisture: time step"),
    )
    for ax, rows, xcol, ycol, xlabel, ylabel, title in panels:
        ax.loglog([r[xcol] for r in rows], [r[ycol] for r in rows], marker="o", color=PALETTE["primary"])
        ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
        if xcol == 2:
            ax.invert_xaxis()
    export_plot(fig, "process_q1_cn_convergence", size=(7.2, 5.4))


def make_new_figures(env_time: np.ndarray, env_temp: np.ndarray, env_c: np.ndarray,
                     cn: dict, be: dict, matched_cn: dict,
                     interpolation_rows: list[dict], convergence_rows: list[list[object]]) -> None:
    apply_publication_style()
    FIGURES.mkdir(exist_ok=True)
    for values, label, filename in (
        (env_temp, "Ambient temperature (degC)", "raw_q1_temperature_interpolation"),
        (env_c, "Ambient moisture (kg/kg)", "raw_q1_moisture_interpolation"),
    ):
        fig, ax = plt.subplots(layout="constrained")
        dense = np.linspace(env_time[0], env_time[-1], 1201)
        for idx, method in enumerate(("linear", "pchip", "cubic", "akima")):
            ax.plot(dense, make_interpolator(env_time, values, method)(dense),
                    color=COLOR_SEQUENCE[idx], linestyle=("-", "--", "-.", ":")[idx], label=method)
        ax.scatter(env_time, values, s=12, color="black", zorder=5, label="measured input")
        ax.set(xlabel="Time (s)", ylabel=label, title="Boundary interpolation comparison")
        ax.legend(ncol=3)
        export_plot(fig, filename)

    fig, axes = plt.subplots(1, 2, layout="constrained")
    times = KEY_TIMES
    for ax, key, ylabel in ((axes[0], "temperature", "Max abs difference (degC)"),
                            (axes[1], "moisture", "Max abs difference (kg/kg)")):
        d = np.max(np.abs(extract(be, key=key) - extract(matched_cn, key=key)), axis=1)
        ax.plot(times, d, marker="o", color=PALETTE["primary"])
        ax.set(xlabel="Time (s)", ylabel=ylabel, title=f"BE versus CN: {key}")
    export_plot(fig, "process_q1_be_cn_difference", size=(7.2, 3.6))

    make_convergence_figure(convergence_rows)

    fig, axes = plt.subplots(2, 1, sharex=True, layout="constrained")
    r = np.asarray(cn["radius_m"]) * 100
    for idx, target in enumerate(KEY_TIMES):
        row = int(target)
        axes[0].plot(r, cn["temperature"][row], color=COLOR_SEQUENCE[idx], label=f"{target:g} s")
        axes[1].plot(r, cn["moisture"][row], color=COLOR_SEQUENCE[idx])
    axes[0].set(ylabel="Temperature (degC)", title="CN radial profiles")
    axes[0].legend(ncol=4)
    axes[1].set(xlabel="Radius (cm)", ylabel="Moisture (kg/kg)")
    export_plot(fig, "result_q1_cn_profiles", size=(7.2, 4.8))


def run_studies(env_time: np.ndarray, env_temp: np.ndarray, env_c: np.ndarray) -> tuple[dict, list[list[object]], dict]:
    runs: dict[tuple[int, float], dict] = {}
    for intervals in (100, 200, 400, 800):
        cfg = RunConfig("cn", intervals, 0.25)
        runs[(intervals, 0.25)] = simulate(env_time, env_temp, env_c, cfg)
    for dt in (1.0, 0.5, 0.125):
        cfg = RunConfig("cn", 400, dt)
        runs[(400, dt)] = simulate(env_time, env_temp, env_c, cfg)
    grid_reference = runs[(800, 0.25)]
    time_reference = runs[(400, 0.125)]
    rows: list[list[object]] = []
    for intervals in (100, 200, 400, 800):
        metrics = difference_metrics(runs[(intervals, 0.25)], grid_reference)
        rows.append(["grid", intervals, 0.25, metrics["temperature_max_abs"], metrics["moisture_max_abs"]])
    for dt in (1.0, 0.5, 0.25, 0.125):
        metrics = difference_metrics(runs[(400, dt)], time_reference)
        rows.append(["time", 400, dt, metrics["temperature_max_abs"], metrics["moisture_max_abs"]])
    selected = simulate(env_time, env_temp, env_c, RunConfig("cn", 800, 0.125))
    return selected, rows, runs


def full_main(no_figures: bool) -> dict:
    started = time.perf_counter()
    environment_file, template_file = legacy.locate_inputs()
    env_time, env_temp, env_c = legacy.load_environment(environment_file)
    RESULTS.mkdir(exist_ok=True)
    temp_interp = interpolation_diagnostics(env_time, env_temp)
    moisture_interp = interpolation_diagnostics(env_time, env_c)
    interpolation_rows = []
    for variable, entries in (("temperature", temp_interp), ("moisture", moisture_interp)):
        for entry in entries:
            interpolation_rows.append([variable, entry["method"], entry["mae"], entry["rmse"],
                                       entry["max_abs_error"], entry["max_interval_overshoot"]])
    write_rows(RESULTS / "问题1_边界插值比较.csv",
               ["variable", "method", "loo_mae", "loo_rmse", "loo_max_abs", "max_interval_overshoot"],
               interpolation_rows)

    cn, convergence_rows, study_runs = run_studies(env_time, env_temp, env_c)
    be = simulate(env_time, env_temp, env_c, RunConfig("be", 200, 1.0))
    matched_cn = simulate(env_time, env_temp, env_c, RunConfig("cn", 200, 1.0))
    arithmetic = simulate(env_time, env_temp, env_c, RunConfig("cn", 400, 0.25, face_mean="arithmetic"))
    pchip = simulate(env_time, env_temp, env_c, RunConfig("cn", 400, 0.25, interpolation="pchip"))
    validate = {"cn": validate_result(cn), "be": validate_result(be), "matched_cn": validate_result(matched_cn)}
    bessel = constant_boundary_validation()

    write_rows(RESULTS / "问题1_CN收敛性.csv",
               ["study", "radial_intervals", "dt_s", "temperature_max_abs", "moisture_max_abs"],
               convergence_rows)
    method_metrics = difference_metrics(be, matched_cn)
    # 只改变径向网格、固定 dt=0.25 s，避免把时间步误差混入“分辨率影响”。
    resolution_metrics = difference_metrics(study_runs[(200, 0.25)], study_runs[(800, 0.25)])
    face_metrics = difference_metrics(arithmetic, study_runs[(400, 0.25)])
    interpolation_metrics = difference_metrics(pchip, study_runs[(400, 0.25)])
    comparison_rows = [
        ["时间方法（相同 N=200, dt=1 s）", *method_metrics.values()],
        ["空间网格（固定dt=0.25 s，N=200与800）", *resolution_metrics.values()],
        ["界面平均（算术与调和）", *face_metrics.values()],
        ["边界插值（PCHIP与线性）", *interpolation_metrics.values()],
    ]
    write_rows(RESULTS / "问题1_差异来源量化.csv",
               ["source", "temperature_max_abs", "temperature_rmse", "moisture_max_abs", "moisture_rmse"], comparison_rows)

    full_times = np.arange(1, 1801, dtype=float)
    cn_t_full = extract(cn, full_times, OUTPUT_RADII_CM, "temperature")
    cn_c_full = extract(cn, full_times, OUTPUT_RADII_CM, "moisture")
    be_t_full = extract(be, full_times, OUTPUT_RADII_CM, "temperature")
    be_c_full = extract(be, full_times, OUTPUT_RADII_CM, "moisture")
    for prefix, t_values, c_values in (("CN", cn_t_full, cn_c_full), ("BE", be_t_full, be_c_full)):
        legacy.write_csv_matrix(RESULTS / f"问题1_{prefix}_完整温度.csv", full_times, OUTPUT_RADII_CM, t_values)
        legacy.write_csv_matrix(RESULTS / f"问题1_{prefix}_完整水分浓度.csv", full_times, OUTPUT_RADII_CM, c_values)
        legacy.write_template_result(template_file, RESULTS / f"result1_{prefix}.xlsx", t_values, c_values)
        legacy.write_csv_matrix(RESULTS / f"问题1_{prefix}_指定时刻温度.csv", KEY_TIMES, KEY_RADII_CM,
                                extract(cn if prefix == "CN" else be, key="temperature"))
        legacy.write_csv_matrix(RESULTS / f"问题1_{prefix}_指定时刻水分浓度.csv", KEY_TIMES, KEY_RADII_CM,
                                extract(cn if prefix == "CN" else be, key="moisture"))
    # 主提交结果使用经收敛检验后的 CN；BE 文件始终单独保留。
    legacy.write_template_result(template_file, RESULTS / "result1.xlsx", cn_t_full, cn_c_full)
    legacy.write_csv_matrix(RESULTS / "问题1_完整温度.csv", full_times, OUTPUT_RADII_CM, cn_t_full)
    legacy.write_csv_matrix(RESULTS / "问题1_完整水分浓度.csv", full_times, OUTPUT_RADII_CM, cn_c_full)
    legacy.write_csv_matrix(RESULTS / "问题1_指定时刻温度.csv", KEY_TIMES, KEY_RADII_CM,
                            extract(cn, key="temperature"))
    legacy.write_csv_matrix(RESULTS / "问题1_指定时刻水分浓度.csv", KEY_TIMES, KEY_RADII_CM,
                            extract(cn, key="moisture"))

    # 用户提供截图中的公开结果，仅用于逐点复核；个别截图版本末位存在差别。
    paper_t = np.array([
        [28.0001, 28.0003, 28.0040, 28.0327, 28.1801],
        [28.0408, 28.0635, 28.1514, 28.3681, 28.8489],
        [28.4533, 28.5360, 28.8040, 29.3159, 30.1652],
        [29.3243, 29.4583, 29.8755, 30.6161, 31.7304],
        [30.5427, 30.7098, 31.2223, 32.1126, 33.4276],
        [31.9957, 32.1867, 32.7660, 33.7463, 35.1203],
        [33.5753, 33.7720, 34.3642, 35.3621, 36.7856],
    ])
    paper_c = np.array([
        [2.5500, 2.5500, 2.5500, 2.5500, 2.2470],
        [2.5500, 2.5500, 2.5500, 2.5492, 2.0508],
        [2.5500, 2.5500, 2.5500, 2.5353, 1.8770],
        [2.5500, 2.5500, 2.5497, 2.5045, 1.7547],
        [2.5500, 2.5500, 2.5482, 2.4646, 1.6586],
        [2.5500, 2.5499, 2.5445, 2.4206, 1.5787],
        [2.5500, 2.5497, 2.5383, 2.3755, 1.5103],
    ])
    cn_key_t, cn_key_c = extract(cn, key="temperature"), extract(cn, key="moisture")
    paper_rows = []
    for i, target_time in enumerate(KEY_TIMES):
        for j, target_radius in enumerate(KEY_RADII_CM):
            paper_rows.append([target_time, target_radius, "temperature", cn_key_t[i, j], paper_t[i, j], cn_key_t[i, j] - paper_t[i, j]])
            paper_rows.append([target_time, target_radius, "moisture", cn_key_c[i, j], paper_c[i, j], cn_key_c[i, j] - paper_c[i, j]])
    write_rows(RESULTS / "问题1_论文截图逐点对比.csv",
               ["time_s", "radius_cm", "variable", "our_CN", "screenshot_value", "difference"], paper_rows)

    if not no_figures:
        make_new_figures(env_time, env_temp, env_c, cn, be, matched_cn, interpolation_rows, convergence_rows)

    summary = {
        "formal_CN": {"method": "Crank-Nicolson", "radial_intervals": 800, "node_count": 801,
                      "dr_m": R / 800, "dt_s": 0.125, "interpolation": "linear", "face_mean": "harmonic"},
        "BE_baseline": {"method": "backward Euler", "radial_intervals": 200, "node_count": 201,
                        "dr_m": R / 200, "dt_s": 1.0},
        "validation": validate,
        "bessel_validation": bessel,
        "difference_sources": {
            "time_method_same_grid": method_metrics,
            "spatial_grid_fixed_dt_0.25s_N200_vs_N800": resolution_metrics,
            "face_mean": face_metrics,
            "boundary_interpolation": interpolation_metrics,
            "paper_screenshot": {
                "temperature_max_abs": float(np.max(np.abs(cn_key_t - paper_t))),
                "moisture_max_abs": float(np.max(np.abs(cn_key_c - paper_c))),
                "note": "values transcribed from user screenshots; different screenshots differ in the last digit",
            },
        },
        "CN_key_temperature": extract(cn, key="temperature").tolist(),
        "CN_key_moisture": extract(cn, key="moisture").tolist(),
        "BE_key_temperature": extract(be, key="temperature").tolist(),
        "BE_key_moisture": extract(be, key="moisture").tolist(),
        "runtime_seconds": time.perf_counter() - started,
    }
    (RESULTS / "问题1_CN_BE校验摘要.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    (RESULTS / "问题1_校验摘要.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    manifest = {
        "command": "python code/q1_solver.py",
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "inputs": {str(environment_file.relative_to(ROOT)): sha256(environment_file),
                   str(template_file.relative_to(ROOT)): sha256(template_file)},
        "code": {"code/q1_solver.py": sha256(ROOT / "code" / "q1_solver.py"),
                 "code/q1_cn_comparison.py": sha256(ROOT / "code" / "q1_cn_comparison.py"),
                 "code/问题1_求解.py": sha256(ROOT / "code" / "问题1_求解.py")},
        "transcribed_reference": {
            "results/问题1_论文截图逐点对比.csv": sha256(RESULTS / "问题1_论文截图逐点对比.csv"),
            "note": "截图数值由用户提供图片人工转录，仅作比较，不作为模型输入",
        },
        "outputs": ["results/result1.xlsx", "results/result1_CN.xlsx", "results/result1_BE.xlsx"],
        "random_seed": None,
        "deterministic": True,
    }
    (RESULTS / "复现清单_CN_BE.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true", help="运行 CN/BE 60 s 最小测试")
    parser.add_argument("--no-figures", action="store_true")
    args = parser.parse_args()
    environment_file, _ = legacy.locate_inputs()
    env_time, env_temp, env_c = legacy.load_environment(environment_file)
    if args.smoke:
        output = {}
        for method in ("be", "cn"):
            cfg = RunConfig(method, 40, 0.5, end_time=60.0, store_interval=1.0)
            result = simulate(env_time, env_temp, env_c, cfg)
            output[method] = {"configuration": {"radial_intervals": cfg.intervals, "dt_s": cfg.dt,
                                                  "theta": cfg.theta, "picard_tolerance": PICARD_TOL},
                              "validation": validate_result(result),
                              "T_center_60s": float(result["temperature"][-1, 0]),
                              "T_surface_60s": float(result["temperature"][-1, -1]),
                              "C_center_60s": float(result["moisture"][-1, 0]),
                              "C_surface_60s": float(result["moisture"][-1, -1])}
        print(json.dumps(output, ensure_ascii=False, indent=2))
        return 0
    print(json.dumps(full_main(args.no_figures), ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
