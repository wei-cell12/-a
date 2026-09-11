#!/usr/bin/env python3
"""2026 A 题问题一：固定半径圆柱径向导热与非线性水分扩散。"""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import sys
import time
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import solve_banded

from utils.plot_style import COLOR_SEQUENCE, PALETTE, apply_publication_style


ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
RADIUS = 0.02
RHO = 820.0
CP = 2600.0
K = 0.36
H = 25.0
HM = 8.0e-7
T0 = 28.0
C0 = 2.55
KEY_TIMES = np.array([100, 300, 600, 900, 1200, 1500, 1800], dtype=float)
KEY_RADII_CM = np.array([0, 0.5, 1.0, 1.5, 2.0], dtype=float)
OUTPUT_RADII_CM = np.round(np.arange(0.0, 2.0001, 0.1), 10)

NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
      "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}


def _xlsx_rows(path: Path, sheet_index: int = 0) -> list[list[object]]:
    """Read raw rows from XLSX with the standard library; no header inference."""
    with zipfile.ZipFile(path) as archive:
        shared: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for item in root.findall("m:si", NS):
                shared.append("".join(node.text or "" for node in item.iter(f"{{{NS['m']}}}t")))
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        relation_map = {item.attrib["Id"]: item.attrib["Target"] for item in relationships}
        sheets = list(workbook.find("m:sheets", NS))
        relation_id = sheets[sheet_index].attrib[f"{{{NS['r']}}}id"]
        target = relation_map[relation_id].lstrip("/")
        if not target.startswith("xl/"):
            target = "xl/" + target
        worksheet = ET.fromstring(archive.read(target))
        rows: list[list[object]] = []
        for row in worksheet.findall(".//m:sheetData/m:row", NS):
            values: dict[int, object] = {}
            for cell in row.findall("m:c", NS):
                reference = cell.attrib["r"]
                letters = "".join(ch for ch in reference if ch.isalpha())
                column = 0
                for letter in letters:
                    column = column * 26 + ord(letter.upper()) - 64
                value_node = cell.find("m:v", NS)
                cell_type = cell.attrib.get("t")
                value: object = None if value_node is None else value_node.text
                if cell_type == "s" and value is not None:
                    value = shared[int(value)]
                elif cell_type == "inlineStr":
                    value = "".join(node.text or "" for node in cell.iter(f"{{{NS['m']}}}t"))
                elif value is not None:
                    value = float(value)
                values[column - 1] = value
            if values:
                width = max(values) + 1
                rows.append([values.get(index) for index in range(width)])
        return rows


def locate_inputs() -> tuple[Path, Path]:
    environment = None
    template = None
    for path in ROOT.glob("**/*.xlsx"):
        try:
            rows = _xlsx_rows(path)
        except Exception:
            continue
        if len(rows) == 242 and len(rows[0]) >= 3:
            environment = path
        if path.name.lower() == "result1.xlsx":
            template = path
    if environment is None or template is None:
        raise FileNotFoundError("无法定位附件1或 result1.xlsx 模板")
    return environment, template


def load_environment(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows = _xlsx_rows(path)
    if len(rows) != 242:
        raise ValueError(f"附件1有效行数应为242（含表头），实际为{len(rows)}")
    numeric = np.asarray([[float(row[0]), float(row[1]), float(row[2])] for row in rows[1:]], dtype=float)
    if numeric[0, 0] != 0 or numeric[-1, 0] != 14400:
        raise ValueError("附件1首末时间不符合题目说明")
    if not np.all(np.diff(numeric[:, 0]) == 60):
        raise ValueError("附件1时间间隔不是统一的60 s")
    selected = numeric[numeric[:, 0] <= 1800]
    if selected.shape != (31, 3):
        raise ValueError(f"0~1800 s 应有31条数据，实际为{selected.shape[0]}")
    return selected[:, 0], selected[:, 1], selected[:, 2]


def diffusivity(moisture: np.ndarray) -> np.ndarray:
    safe = np.maximum(moisture, 1.0e-12)
    return 7.0e-9 * np.exp(-0.89 / safe)


def radial_geometry(node_count: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    radii = np.linspace(0.0, RADIUS, node_count)
    dr = radii[1] - radii[0]
    west_faces = np.maximum(radii - 0.5 * dr, 0.0)
    east_faces = np.minimum(radii + 0.5 * dr, RADIUS)
    volumes = 0.5 * (east_faces ** 2 - west_faces ** 2)
    return radii, volumes, east_faces


def build_implicit_system(
    coefficient: np.ndarray,
    capacity: float,
    dt: float,
    boundary_coefficient: float,
    boundary_value: float,
    radii: np.ndarray,
    volumes: np.ndarray,
    east_faces: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    n = radii.size
    dr = radii[1] - radii[0]
    faces = east_faces[:-1]
    harmonic = 2.0 * coefficient[:-1] * coefficient[1:] / np.maximum(
        coefficient[:-1] + coefficient[1:], 1.0e-300
    )
    conductance = faces * harmonic / dr
    diagonal = np.ones(n)
    lower = np.zeros(n)
    upper = np.zeros(n)
    source = np.zeros(n)
    for i in range(n):
        west = conductance[i - 1] if i > 0 else 0.0
        east = conductance[i] if i < n - 1 else 0.0
        boundary = RADIUS * boundary_coefficient if i == n - 1 else 0.0
        scale = dt / (capacity * volumes[i])
        diagonal[i] += scale * (west + east + boundary)
        if i > 0:
            lower[i] = -scale * west
        if i < n - 1:
            upper[i] = -scale * east
        if i == n - 1:
            source[i] = scale * boundary * boundary_value
    banded = np.zeros((3, n))
    banded[0, 1:] = upper[:-1]
    banded[1, :] = diagonal
    banded[2, :-1] = lower[1:]
    return banded, source


def area_average(field: np.ndarray, volumes: np.ndarray) -> float:
    return float(np.dot(field, volumes) / (0.5 * RADIUS ** 2))


def simulate(
    env_time: np.ndarray,
    env_temp: np.ndarray,
    env_moisture: np.ndarray,
    *,
    node_count: int = 201,
    dt: float = 1.0,
    end_time: float = 1800.0,
    store_interval: float = 1.0,
) -> dict[str, np.ndarray | float | int]:
    steps = int(round(end_time / dt))
    if not math.isclose(steps * dt, end_time, abs_tol=1e-12):
        raise ValueError("end_time 必须为 dt 的整数倍")
    store_stride = int(round(store_interval / dt))
    if not math.isclose(store_stride * dt, store_interval, abs_tol=1e-12):
        raise ValueError("store_interval 必须为 dt 的整数倍")
    radii, volumes, east_faces = radial_geometry(node_count)
    temperature = np.full(node_count, T0)
    moisture = np.full(node_count, C0)
    out_times = np.arange(0.0, end_time + 0.5 * store_interval, store_interval)
    out_t = np.empty((out_times.size, node_count))
    out_c = np.empty_like(out_t)
    out_t[0] = temperature
    out_c[0] = moisture
    temperature_flux_integral = 0.0
    moisture_flux_integral = 0.0
    temperature_residuals = []
    moisture_residuals = []
    maximum_picard_iterations = 0

    thermal_matrix, _ = build_implicit_system(
        np.full(node_count, K), RHO * CP, dt, H, 0.0,
        radii, volumes, east_faces,
    )
    store_index = 1
    for step in range(1, steps + 1):
        current_time = step * dt
        ambient_t = float(np.interp(current_time, env_time, env_temp))
        ambient_c = float(np.interp(current_time, env_time, env_moisture))
        _, thermal_source = build_implicit_system(
            np.full(node_count, K), RHO * CP, dt, H, ambient_t,
            radii, volumes, east_faces,
        )
        previous_average_t = area_average(temperature, volumes)
        previous_average_c = area_average(moisture, volumes)
        temperature = solve_banded((1, 1), thermal_matrix, temperature + thermal_source)

        guess = moisture.copy()
        converged = False
        for iteration in range(1, 51):
            matrix, source = build_implicit_system(
                diffusivity(guess), 1.0, dt, HM, ambient_c,
                radii, volumes, east_faces,
            )
            updated = solve_banded((1, 1), matrix, moisture + source)
            if np.max(np.abs(updated - guess)) < 1.0e-11:
                converged = True
                guess = updated
                maximum_picard_iterations = max(maximum_picard_iterations, iteration)
                break
            guess = updated
        if not converged:
            raise RuntimeError(f"水分 Picard 迭代在 t={current_time} s 未收敛")
        moisture = guess

        current_average_t = area_average(temperature, volumes)
        current_average_c = area_average(moisture, volumes)
        delta_t_flux = dt * (2.0 * H / (RHO * CP * RADIUS)) * (ambient_t - temperature[-1])
        delta_c_flux = dt * (2.0 * HM / RADIUS) * (ambient_c - moisture[-1])
        temperature_flux_integral += delta_t_flux
        moisture_flux_integral += delta_c_flux
        temperature_residuals.append((current_average_t - previous_average_t) - delta_t_flux)
        moisture_residuals.append((current_average_c - previous_average_c) - delta_c_flux)

        if step % store_stride == 0:
            out_t[store_index] = temperature
            out_c[store_index] = moisture
            store_index += 1

    return {
        "time": out_times,
        "radius_m": radii,
        "temperature": out_t,
        "moisture": out_c,
        "temperature_conservation_abs": float(np.max(np.abs(temperature_residuals))),
        "moisture_conservation_abs": float(np.max(np.abs(moisture_residuals))),
        "temperature_cumulative_residual": float(
            area_average(temperature, volumes) - T0 - temperature_flux_integral
        ),
        "moisture_cumulative_residual": float(
            area_average(moisture, volumes) - C0 - moisture_flux_integral
        ),
        "max_picard_iterations": maximum_picard_iterations,
    }


def interpolate_result(result: dict, times: np.ndarray, radii_cm: np.ndarray, key: str) -> np.ndarray:
    source_times = np.asarray(result["time"])
    source_radii = np.asarray(result["radius_m"]) * 100.0
    field = np.asarray(result[key])
    rows = []
    for target_time in times:
        time_index = int(np.argmin(np.abs(source_times - target_time)))
        if not math.isclose(source_times[time_index], target_time, abs_tol=1e-9):
            raise ValueError("目标时刻不在存储网格上")
        rows.append(np.interp(radii_cm, source_radii, field[time_index]))
    return np.asarray(rows)


def write_csv_matrix(path: Path, times: np.ndarray, radii_cm: np.ndarray, values: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["time_s", *[f"r_{value:g}_cm" for value in radii_cm]])
        for time_value, row in zip(times, values):
            writer.writerow([f"{time_value:g}", *[f"{value:.4f}" for value in row]])


def excel_column(index: int) -> str:
    letters = ""
    while index:
        index, remainder = divmod(index - 1, 26)
        letters = chr(65 + remainder) + letters
    return letters


def build_sheet_xml(values: np.ndarray, numeric_style: int) -> bytes:
    rows = []
    header_cells = ['<c r="A1" s="2" t="inlineStr"><is><t>time / radius (cm)</t></is></c>']
    for column, radius in enumerate(OUTPUT_RADII_CM, start=2):
        header_cells.append(f'<c r="{excel_column(column)}1" s="2"><v>{radius:g}</v></c>')
    rows.append(f'<row r="1" spans="1:22">{"".join(header_cells)}</row>')
    for row_number, (time_value, row) in enumerate(zip(range(1, 1801), values), start=2):
        cells = [f'<c r="A{row_number}" s="1"><v>{time_value}</v></c>']
        for column, value in enumerate(row, start=2):
            cells.append(
                f'<c r="{excel_column(column)}{row_number}" s="{numeric_style}"><v>{value:.4f}</v></c>'
            )
        rows.append(f'<row r="{row_number}" spans="1:22">{"".join(cells)}</row>')
    xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        '<dimension ref="A1:V1801"/><sheetViews><sheetView workbookViewId="0"/></sheetViews>'
        '<sheetFormatPr defaultColWidth="8.75" defaultRowHeight="14.1"/>'
        '<cols><col min="1" max="1" width="19.625" customWidth="1"/>'
        '<col min="2" max="22" width="9.625" customWidth="1"/></cols>'
        f'<sheetData>{"".join(rows)}</sheetData>'
        '<pageMargins left="0.75" right="0.75" top="1" bottom="1" header="0.5" footer="0.5"/>'
        '</worksheet>'
    )
    return xml.encode("utf-8")


def add_four_decimal_style(styles_xml: bytes) -> tuple[bytes, int]:
    namespace = NS["m"]
    ET.register_namespace("", namespace)
    root = ET.fromstring(styles_xml)
    num_formats = root.find(f"{{{namespace}}}numFmts")
    if num_formats is None:
        num_formats = ET.Element(f"{{{namespace}}}numFmts", {"count": "1"})
        root.insert(0, num_formats)
    else:
        num_formats.attrib["count"] = str(int(num_formats.attrib.get("count", "0")) + 1)
    ET.SubElement(num_formats, f"{{{namespace}}}numFmt", {"numFmtId": "164", "formatCode": "0.0000"})
    cell_xfs = root.find(f"{{{namespace}}}cellXfs")
    style_index = len(list(cell_xfs))
    ET.SubElement(cell_xfs, f"{{{namespace}}}xf", {
        "numFmtId": "164", "fontId": "0", "fillId": "0", "borderId": "0",
        "xfId": "0", "applyNumberFormat": "1", "applyAlignment": "1",
    }).append(ET.Element(f"{{{namespace}}}alignment", {"horizontal": "center", "vertical": "center"}))
    cell_xfs.attrib["count"] = str(style_index + 1)
    return ET.tostring(root, encoding="utf-8", xml_declaration=True), style_index


def write_template_result(template: Path, output: Path, temperature: np.ndarray, moisture: np.ndarray) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(".tmp.xlsx")
    with zipfile.ZipFile(template) as source, zipfile.ZipFile(temporary, "w", zipfile.ZIP_DEFLATED) as target:
        styles, numeric_style = add_four_decimal_style(source.read("xl/styles.xml"))
        replacements = {
            "xl/styles.xml": styles,
            "xl/worksheets/sheet1.xml": build_sheet_xml(temperature, numeric_style),
            "xl/worksheets/sheet2.xml": build_sheet_xml(moisture, numeric_style),
        }
        for item in source.infolist():
            target.writestr(item, replacements.get(item.filename, source.read(item.filename)))
    temporary.replace(output)


def export_plot(fig, basename: str, size: tuple[float, float] = (6.3, 3.9)) -> None:
    skill_scripts = Path(r"C:\Users\17216\.codex\skills\math-modeling\tools\figure\scripts")
    sys.path.insert(0, str(skill_scripts))
    from export_figure import export_figure
    paths = export_figure(fig, str(FIGURES / basename), formats=["svg", "png"], dpi=300,
                          size_inches=size, grayscale_preview=False, tight=False)
    qa_directory = FIGURES / "_qa"
    qa_directory.mkdir(exist_ok=True)
    from PIL import Image
    png_path = next(Path(exported) for exported in paths if str(exported).endswith(".png"))
    grayscale_path = qa_directory / f"{basename}_grayscale.png"
    Image.open(png_path).convert("L").save(grayscale_path, dpi=(300, 300))
    plt.close(fig)


def make_figures(
    env_time: np.ndarray,
    env_temp: np.ndarray,
    env_moisture: np.ndarray,
    result: dict,
    convergence: dict,
) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    apply_publication_style()
    r_cm = np.asarray(result["radius_m"]) * 100.0
    times = np.asarray(result["time"])
    temperature = np.asarray(result["temperature"])
    moisture = np.asarray(result["moisture"])

    fig, ax = plt.subplots(layout="constrained")
    ax.plot(env_time, env_temp, color=PALETTE["primary"], marker="o", markevery=5)
    ax.set(xlabel="Time (s)", ylabel="Ambient temperature (degC)", title="Boundary temperature")
    export_plot(fig, "raw_q1_boundary_temperature")

    fig, axes = plt.subplots(1, 2, layout="constrained")
    axes[0].hist(env_moisture, bins=8, color=PALETTE["sky"], edgecolor="white")
    axes[0].set(xlabel="Moisture input (kg/kg)", ylabel="Count", title="Distribution")
    axes[1].boxplot(env_moisture, vert=True, widths=0.4)
    axes[1].set(ylabel="Moisture input (kg/kg)", xticks=[], title="Range")
    export_plot(fig, "raw_q1_boundary_moisture_distribution")

    fig, ax = plt.subplots(layout="constrained")
    scatter = ax.scatter(env_temp, env_moisture, c=env_time, cmap="viridis", s=24)
    colorbar = fig.colorbar(scatter, ax=ax)
    colorbar.set_label("Time (s)")
    ax.set(xlabel="Ambient temperature (degC)", ylabel="Moisture input (kg/kg)", title="Boundary-input relation")
    export_plot(fig, "raw_q1_boundary_relation")

    profile_times = [300, 900, 1800]
    fig, ax = plt.subplots(layout="constrained")
    for index, target in enumerate(profile_times):
        ax.plot(r_cm, temperature[target], color=COLOR_SEQUENCE[index], label=f"{target} s",
                linestyle=["-", "--", "-."][index])
    ax.set(xlabel="Radius (cm)", ylabel="Temperature (degC)", title="Radial temperature profiles")
    ax.legend()
    export_plot(fig, "process_q1_temperature_profiles")

    fig, ax = plt.subplots(layout="constrained")
    for index, target in enumerate(profile_times):
        ax.plot(r_cm, moisture[target], color=COLOR_SEQUENCE[index], label=f"{target} s",
                linestyle=["-", "--", "-."][index])
    ax.set(xlabel="Radius (cm)", ylabel="Dry-basis moisture (kg/kg)", title="Radial moisture profiles")
    ax.legend()
    export_plot(fig, "process_q1_moisture_profiles")

    fig, axes = plt.subplots(2, 2, layout="constrained")
    categories = ["Grid", "Time step"]
    axes[0, 0].bar(categories,
                   [convergence["grid_temperature"], convergence["time_temperature"]],
                   color=[PALETTE["primary"], PALETTE["secondary"]])
    axes[0, 0].set(ylabel="Maximum difference (degC)", title="Temperature discretization")
    axes[0, 1].bar(categories,
                   [convergence["grid_moisture"], convergence["time_moisture"]],
                   color=[PALETTE["primary"], PALETTE["secondary"]])
    axes[0, 1].set(ylabel="Maximum difference (kg/kg)", title="Moisture discretization")
    axes[1, 0].bar(["Heat balance"], [abs(result["temperature_cumulative_residual"])],
                   color=PALETTE["positive"], width=0.55)
    axes[1, 0].set_yscale("log")
    axes[1, 0].set(ylabel="Cumulative residual (degC)", title="Heat conservation")
    axes[1, 1].bar(["Mass balance"], [abs(result["moisture_cumulative_residual"])],
                   color=PALETTE["accent"], width=0.55)
    axes[1, 1].set_yscale("log")
    axes[1, 1].set(ylabel="Cumulative residual (kg/kg)", title="Moisture conservation")
    export_plot(fig, "process_q1_numerical_verification", size=(7.2, 5.4))

    fig, ax = plt.subplots(layout="constrained")
    image = ax.imshow(temperature, origin="lower", aspect="auto", cmap="inferno",
                      extent=[0, 2, 0, 1800])
    colorbar = fig.colorbar(image, ax=ax)
    colorbar.set_label("Temperature (degC)")
    ax.set(xlabel="Radius (cm)", ylabel="Time (s)", title="Temperature field")
    export_plot(fig, "result_q1_temperature_field")

    fig, ax = plt.subplots(layout="constrained")
    image = ax.imshow(moisture, origin="lower", aspect="auto", cmap="viridis",
                      extent=[0, 2, 0, 1800])
    colorbar = fig.colorbar(image, ax=ax)
    colorbar.set_label("Dry-basis moisture (kg/kg)")
    ax.set(xlabel="Radius (cm)", ylabel="Time (s)", title="Moisture field")
    export_plot(fig, "result_q1_moisture_field")

    fig, axes = plt.subplots(2, 1, sharex=True, layout="constrained")
    selected_radii = [0.0, 0.5, 1.0, 1.5, 2.0]
    for index, radius in enumerate(selected_radii):
        node = int(round(radius / (r_cm[1] - r_cm[0])))
        line_style = ["-", "--", "-.", ":", (0, (5, 2))][index]
        axes[0].plot(times, temperature[:, node], label=f"r={radius:g} cm",
                     color=COLOR_SEQUENCE[index], linestyle=line_style)
        axes[1].plot(times, moisture[:, node], color=COLOR_SEQUENCE[index], linestyle=line_style)
    axes[0].set(ylabel="Temperature (degC)", title="Responses at selected radii")
    axes[0].legend(ncol=5, loc="upper center")
    axes[1].set(xlabel="Time (s)", ylabel="Moisture (kg/kg)")
    export_plot(fig, "result_q1_selected_positions", size=(7.2, 4.8))


def write_contract() -> None:
    rows = [
        ("raw_q1_boundary_temperature", "raw", "Environment temperature rises nonlinearly", "line", "6.3x3.9 in"),
        ("raw_q1_boundary_moisture_distribution", "raw", "Boundary moisture input has a narrow changing range", "histogram+box", "6.3x3.9 in"),
        ("raw_q1_boundary_relation", "raw", "Temperature and moisture inputs co-vary over time", "scatter", "6.3x3.9 in"),
        ("process_q1_temperature_profiles", "process", "Surface heating precedes the cylinder center", "line profiles", "6.3x3.9 in"),
        ("process_q1_moisture_profiles", "process", "Surface drying precedes the cylinder center", "line profiles", "6.3x3.9 in"),
        ("process_q1_numerical_verification", "process", "Temperature and moisture errors are quantified without cross-unit comparison", "four unit-specific bar panels", "7.2x5.4 in"),
        ("result_q1_temperature_field", "result", "The radial temperature field shows inward propagation", "heatmap", "6.3x3.9 in"),
        ("result_q1_moisture_field", "result", "The moisture field retains a wet core", "heatmap", "6.3x3.9 in"),
        ("result_q1_selected_positions", "result", "Responses are ordered by radius throughout preheating", "stacked lines", "7.2x4.8 in"),
    ]
    with (RESULTS / "图表契约.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["figure", "category", "core_claim", "chart_type", "final_size"])
        writer.writerows(rows)


def validate_result(result: dict) -> dict[str, float]:
    temperature = np.asarray(result["temperature"])
    moisture = np.asarray(result["moisture"])
    metrics = {
        "temperature_min": float(np.min(temperature)),
        "temperature_max": float(np.max(temperature)),
        "moisture_min": float(np.min(moisture)),
        "moisture_max": float(np.max(moisture)),
        "temperature_cumulative_residual": float(result["temperature_cumulative_residual"]),
        "moisture_cumulative_residual": float(result["moisture_cumulative_residual"]),
        "max_picard_iterations": int(result["max_picard_iterations"]),
    }
    if metrics["temperature_min"] < T0 - 1e-9 or metrics["moisture_min"] < 0:
        raise AssertionError("物理范围检查失败")
    if metrics["temperature_max"] > 41.513 + 1e-9:
        raise AssertionError("温度最大值原理检查失败")
    if abs(metrics["temperature_cumulative_residual"]) > 1e-9:
        raise AssertionError("热量离散守恒检查失败")
    if abs(metrics["moisture_cumulative_residual"]) > 1e-9:
        raise AssertionError("水分离散守恒检查失败")
    return metrics


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true", help="运行60 s、41节点的最小纵向切片")
    parser.add_argument("--no-figures", action="store_true")
    args = parser.parse_args()
    start = time.perf_counter()
    environment_file, template_file = locate_inputs()
    env_time, env_temp, env_moisture = load_environment(environment_file)
    if args.smoke:
        result = simulate(env_time, env_temp, env_moisture, node_count=41, dt=1.0, end_time=60.0)
        metrics = validate_result(result)
        print(json.dumps({"mode": "smoke", "T_center_60s": float(result["temperature"][-1, 0]),
                          "T_surface_60s": float(result["temperature"][-1, -1]),
                          "C_center_60s": float(result["moisture"][-1, 0]),
                          "C_surface_60s": float(result["moisture"][-1, -1]), **metrics},
                         ensure_ascii=False, indent=2))
        return 0

    RESULTS.mkdir(exist_ok=True)
    result = simulate(env_time, env_temp, env_moisture, node_count=201, dt=1.0)
    metrics = validate_result(result)
    coarse = simulate(env_time, env_temp, env_moisture, node_count=101, dt=1.0)
    half_step = simulate(env_time, env_temp, env_moisture, node_count=201, dt=0.5)
    fine_key_t = interpolate_result(result, KEY_TIMES, KEY_RADII_CM, "temperature")
    fine_key_c = interpolate_result(result, KEY_TIMES, KEY_RADII_CM, "moisture")
    coarse_key_t = interpolate_result(coarse, KEY_TIMES, KEY_RADII_CM, "temperature")
    coarse_key_c = interpolate_result(coarse, KEY_TIMES, KEY_RADII_CM, "moisture")
    half_key_t = interpolate_result(half_step, KEY_TIMES, KEY_RADII_CM, "temperature")
    half_key_c = interpolate_result(half_step, KEY_TIMES, KEY_RADII_CM, "moisture")
    convergence = {
        "grid_temperature": float(np.max(np.abs(fine_key_t - coarse_key_t))),
        "grid_moisture": float(np.max(np.abs(fine_key_c - coarse_key_c))),
        "time_temperature": float(np.max(np.abs(fine_key_t - half_key_t))),
        "time_moisture": float(np.max(np.abs(fine_key_c - half_key_c))),
    }

    write_csv_matrix(RESULTS / "问题1_指定时刻温度.csv", KEY_TIMES, KEY_RADII_CM, fine_key_t)
    write_csv_matrix(RESULTS / "问题1_指定时刻水分浓度.csv", KEY_TIMES, KEY_RADII_CM, fine_key_c)
    full_times = np.arange(1, 1801, dtype=float)
    full_t = interpolate_result(result, full_times, OUTPUT_RADII_CM, "temperature")
    full_c = interpolate_result(result, full_times, OUTPUT_RADII_CM, "moisture")
    write_csv_matrix(RESULTS / "问题1_完整温度.csv", full_times, OUTPUT_RADII_CM, full_t)
    write_csv_matrix(RESULTS / "问题1_完整水分浓度.csv", full_times, OUTPUT_RADII_CM, full_c)
    write_template_result(template_file, RESULTS / "result1.xlsx", full_t, full_c)
    with (RESULTS / "environment_q1.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["time_s", "ambient_temperature_degC", "boundary_moisture_kg_per_kg"])
        writer.writerows(zip(env_time, env_temp, env_moisture))
    write_contract()
    if not args.no_figures:
        make_figures(env_time, env_temp, env_moisture, result, convergence)
    summary = {
        "parameters": {"node_count": 201, "dr_m": 0.0001, "dt_s": 1.0, "picard_tolerance": 1e-11},
        "validation": metrics,
        "convergence": convergence,
        "key_temperature": fine_key_t.tolist(),
        "key_moisture": fine_key_c.tolist(),
        "runtime_seconds": time.perf_counter() - start,
    }
    (RESULTS / "问题1_校验摘要.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
