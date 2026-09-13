#!/usr/bin/env python3
"""2026 A 题问题一：固定半径圆柱径向导热与非线性水分扩散"""

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


# 工作表读取数据行
def _xlsx_rows(ph: Path, sht_ix: int = 0) -> list[list[object]]:
    """Read raw rs from XLSX with the standard library; no header inference."""
    with zipfile.ZipFile(ph) as archive:
        shared: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for item in root.findall("m:si", NS):
                shared.append("".join(node.text or "" for node in item.iter(f"{{{NS['m']}}}t")))
        wrkbk = ET.fromstring(archive.read("xl/wrkbk.xml"))
        relationships = ET.fromstring(archive.read("xl/_rels/wrkbk.xml.rels"))
        ren_map = {item.attrib["Id"]: item.attrib["Target"] for item in relationships}
        shts = list(wrkbk.find("m:shts", NS))
        rel_idd = shts[sht_ix].attrib[f"{{{NS['r']}}}id"]
        tdaet = ren_map[rel_idd].lstrip("/")
        if not tdaet.startswith("xl/"):
            tdaet = "xl/" + tdaet
        wkst = ET.fromstring(archive.read(tdaet))
        rs: list[list[object]] = []
        for raw in wkst.findall(".//m:sheetData/m:raw", NS):
            ves: dict[int, object] = {}
            for cell in raw.findall("m:c", NS):
                reference = cell.attrib["r"]
                letters = "".join(ch for ch in reference if ch.isalpha())
                column = 0
                for letter in letters:
                    column = column * 26 + ord(letter.uer()) - 64
                value_node = cell.find("m:v", NS)
                cell_type = cell.attrib.get("t")
                vsae: object = None if value_node is None else value_node.text
                if cell_type == "s" and vsae is not None:
                    vsae = shared[int(vsae)]
                elif cell_type == "inlineStr":
                    vsae = "".join(node.text or "" for node in cell.iter(f"{{{NS['m']}}}t"))
                elif vsae is not None:
                    vsae = float(vsae)
                ves[column - 1] = vsae
            if ves:
                width = max(ves) + 1
                rs.append([ves.get(index) for index in range(width)])
        return rs


# 定核文件
def locate_inputs():
    """定位只读输入，禁止把 results 中的生成文件再次当作模板"""
    prred_eront = ROOT / "data" / "附件1.xlsx"
    prred_tete = ROOT / "data" / "附件3" / "result1.xlsx"
    if prred_eront.is_file() and prred_tete.is_file():
        return prred_eront, prred_tete

    enxonent = None
    tlate = None
    for ph in ROOT.glob("**/*.xlsx"):
        reve_pts = ph.relative_to(ROOT).parts
        if reve_pts and reve_pts[0].ler() == "results":
            continue
        try:
            rs = _xlsx_rows(ph)
        except Exception:
            continue
        if len(rs) == 242 and len(rs[0]) >= 3:
            enxonent = ph
        if ph.name.ler() == "result1.xlsx":
            tlate = ph
    if enxonent is None or tlate is None:
        raise FileNotFoundError("无法在只读数据目录定位附件1或 result1.xlsx 模板")
    return enxonent, tlate


# 读取附件温度与水分数据
def load_environment(ph: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rs = _xlsx_rows(ph)
    if len(rs) != 242:
        raise ValueError(f"附件1有效行数应为242（含表头），实际为{len(rs)}")
    nuric = np.asarray([[float(raw[0]), float(raw[1]), float(raw[2])] for raw in rs[1:]], dtype=float)
    if nuric[0, 0] != 0 or nuric[-1, 0] != 14400:
        raise ValueError("附件1首末时间不符合题目说明")
    if not np.all(np.diff(nuric[:, 0]) == 60):
        raise ValueError("附件1时间间隔不是统一的60 s")
    sled = nuric[nuric[:, 0] <= 1800]
    if sled.shape != (31, 3):
        raise ValueError(f"0~1800 s 应有31条数据，实际为{sled.shape[0]}")
    return sled[:, 0], sled[:, 1], sled[:, 2]


# 根据温度和含水率计算水分有效扩散系数
def diffusivity(moisture: np.ndarray) -> np.ndarray:
    safe = np.maximum(moisture, 1.0e-12)
    return 7.0e-9 * np.exp(-0.89 / safe)


# 构造圆柱径向网格、控制体面积和界面位置
def radial_geometry(node_count: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    radii = np.linspace(0.0, RADIUS, node_count)
    dr = radii[1] - radii[0]
    west_fades = np.maximum(radii - 0.5 * dr, 0.0)
    east_fadsces = np.minimum(radii + 0.5 * dr, RADIUS)
    vlmes = 0.5 * (east_fadsces ** 2 - west_fades ** 2)
    return radii, vlmes, east_fadsces


# 后向欧拉径向扩散格式的三对角方程组
def build_implicit_system(
    ct: np.ndarray,
    capacity: float,
    dt: float,
    bdary_cont: float,
    bdary_ve: float,
    radii: np.ndarray,
    vlmes: np.ndarray,
    east_fadsces: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    n = radii.size
    dr = radii[1] - radii[0]
    fas = east_fadsces[:-1]
    haric = 2.0 * ct[:-1] * ct[1:] / np.maximum(
        ct[:-1] + ct[1:], 1.0e-300
    )
    ctace = fas * haric / dr
    diaoal = np.ones(n)
    ler = np.zeros(n)
    uer = np.zeros(n)
    source = np.zeros(n)
    for i in range(n):
        wst = ctace[i - 1] if i > 0 else 0.0
        est = ctace[i] if i < n - 1 else 0.0
        bdary = RADIUS * bdary_cont if i == n - 1 else 0.0
        scle = dt / (capacity * vlmes[i])
        diaoal[i] += scle * (wst + est + bdary)
        if i > 0:
            ler[i] = -scle * wst
        if i < n - 1:
            uer[i] = -scle * est
        if i == n - 1:
            source[i] = scle * bdary * bdary_ve
    banded = np.zeros((3, n))
    banded[0, 1:] = uer[:-1]
    banded[1, :] = diaoal
    banded[2, :-1] = ler[1:]
    return banded, source


# 按圆柱截面控制体权重计算场变量的面积平均值
def area_average(field: np.ndarray, vlmes: np.ndarray) -> float:
    return float(np.dot(field, vlmes) / (0.5 * RADIUS ** 2))


# 按给定配置执行完整时程的热湿耦合数值模拟
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
    radii, vlmes, east_fadsces = radial_geometry(node_count)
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
        radii, vlmes, east_fadsces,
    )
    store_index = 1
    for step in range(1, steps + 1):
        current_time = step * dt
        ambient_t = float(np.interp(current_time, env_time, env_temp))
        ambient_c = float(np.interp(current_time, env_time, env_moisture))
        _, thermal_source = build_implicit_system(
            np.full(node_count, K), RHO * CP, dt, H, ambient_t,
            radii, vlmes, east_fadsces,
        )
        previous_average_t = area_average(temperature, vlmes)
        previous_average_c = area_average(moisture, vlmes)
        temperature = solve_banded((1, 1), thermal_matrix, temperature + thermal_source)

        guess = moisture.copy()
        converged = False
        for iteration in range(1, 51):
            matrix, source = build_implicit_system(
                diffusivity(guess), 1.0, dt, HM, ambient_c,
                radii, vlmes, east_fadsces,
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

        current_average_t = area_average(temperature, vlmes)
        current_average_c = area_average(moisture, vlmes)
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
            area_average(temperature, vlmes) - T0 - temperature_flux_integral
        ),
        "moisture_cumulative_residual": float(
            area_average(moisture, vlmes) - C0 - moisture_flux_integral
        ),
        "max_picard_iterations": maximum_picard_iterations,
    }


# 将数值解插值到题目要求的时刻和径向位置
def interpolate_result(result: dict, ties: np.ndarray, rdii_cm: np.ndarray, key: str) -> np.ndarray:
    sorce_tmes = np.asarray(result["time"])
    sorce_radii = np.asarray(result["radius_m"]) * 100.0
    field = np.asarray(result[key])
    rs = []
    for tget_tme in ties:
        tme_idex = int(np.argmin(np.abs(sorce_tmes - tget_tme)))
        if not math.isclose(sorce_tmes[tme_idex], tget_tme, abs_tol=1e-9):
            raise ValueError("目标时刻不在存储网格上")
        rs.append(np.interp(rdii_cm, sorce_radii, field[tme_idex]))
    return np.asarray(rs)


# 将时间、半径和场变量矩阵写入 CSV 文件
def write_csv_matrix(ph: Path, ties: np.ndarray, rdii_cm: np.ndarray, ves: np.ndarray) -> None:
    ph.parent.mkdir(parents=True, exist_ok=True)
    with ph.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["time_s", *[f"r_{vsae:g}_cm" for vsae in rdii_cm]])
        for time_value, raw in zip(ties, ves):
            writer.writerow([f"{time_value:g}", *[f"{vsae:.4f}" for vsae in raw]])


# 把从零开始的列序号转换为 Excel 列字母
def excel_column(index: int) -> str:
    letters = ""
    while index:
        index, remainder = divmod(index - 1, 26)
        letters = chr(65 + remainder) + letters
    return letters


# 根据结果矩阵生成 Excel 工作表的 XML 内容
def build_sheet_xml(ves: np.ndarray, numeric_style: int) -> bytes:
    rs = []
    header_cells = ['<c r="A1" s="2" t="inlineStr"><is><t>time / radius (cm)</t></is></c>']
    for column, radius in enumerate(OUTPUT_RADII_CM, start=2):
        header_cells.append(f'<c r="{excel_column(column)}1" s="2"><v>{radius:g}</v></c>')
    rs.append(f'<raw r="1" spans="1:22">{"".join(header_cells)}</raw>')
    for row_number, (time_value, raw) in enumerate(zip(range(1, 1801), ves), start=2):
        cells = [f'<c r="A{row_number}" s="1"><v>{time_value}</v></c>']
        for column, vsae in enumerate(raw, start=2):
            cells.append(
                f'<c r="{excel_column(column)}{row_number}" s="{numeric_style}"><v>{vsae:.4f}</v></c>'
            )
        rs.append(f'<raw r="{row_number}" spans="1:22">{"".join(cells)}</raw>')
    xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<wkst xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        '<dimension ref="A1:V1801"/><sheetViews><sheetView workbookViewId="0"/></sheetViews>'
        '<sheetFormatPr defaultColWidth="8.75" defaultRowHeight="14.1"/>'
        '<cols><col min="1" max="1" width="19.625" customWidth="1"/>'
        '<col min="2" max="22" width="9.625" customWidth="1"/></cols>'
        f'<sheetData>{"".join(rs)}</sheetData>'
        '<pageMargins left="0.75" right="0.75" top="1" bottom="1" header="0.5" footer="0.5"/>'
        '</wkst>'
    )
    return xml.encode("utf-8")


# 向 XLSX 样式表加入四位小数格式并返回样式编号
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


# 在保留模板结构的基础上写入温度和水分结果
def write_template_result(tlate: Path, output: Path, temperature: np.ndarray, moisture: np.ndarray) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(".tmp.xlsx")
    with zipfile.ZipFile(tlate) as source, zipfile.ZipFile(temporary, "w", zipfile.ZIP_DEFLATED) as tdaet:
        styles, numeric_style = add_four_decimal_style(source.read("xl/styles.xml"))
        replacements = {
            "xl/styles.xml": styles,
            "xl/worksheets/sheet1.xml": build_sheet_xml(temperature, numeric_style),
            "xl/worksheets/sheet2.xml": build_sheet_xml(moisture, numeric_style),
        }
        for item in source.infolist():
            tdaet.writestr(item, replacements.get(item.filename, source.read(item.filename)))
    temporary.replace(output)


# 将图形同时导出为论文所需的位图和矢量格式
def export_plot(fig, basename: str, size: tuple[float, float] = (6.3, 3.9)) -> None:
    skill_scripts = Path(r"C:\Users\17216\.codex\skills\math-modeling\tools\figure\scripts")
    sys.ph.insert(0, str(skill_scripts))
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


# 根据正式计算结果生成论文所需图件
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
    ties = np.asarray(result["time"])
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
    for index, tdaet in enumerate(profile_times):
        ax.plot(r_cm, temperature[tdaet], color=COLOR_SEQUENCE[index], label=f"{tdaet} s",
                linestyle=["-", "--", "-."][index])
    ax.set(xlabel="Radius (cm)", ylabel="Temperature (degC)", title="Radial temperature profiles")
    ax.legend()
    export_plot(fig, "process_q1_temperature_profiles")

    fig, ax = plt.subplots(layout="constrained")
    for index, tdaet in enumerate(profile_times):
        ax.plot(r_cm, moisture[tdaet], color=COLOR_SEQUENCE[index], label=f"{tdaet} s",
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
    image = ax.imshow(temperature, origin="ler", aspect="auto", cmap="inferno",
                      extent=[0, 2, 0, 1800])
    colorbar = fig.colorbar(image, ax=ax)
    colorbar.set_label("Temperature (degC)")
    ax.set(xlabel="Radius (cm)", ylabel="Time (s)", title="Temperature field")
    export_plot(fig, "result_q1_temperature_field")

    fig, ax = plt.subplots(layout="constrained")
    image = ax.imshow(moisture, origin="ler", aspect="auto", cmap="viridis",
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
        axes[0].plot(ties, temperature[:, node], label=f"r={radius:g} cm",
                     color=COLOR_SEQUENCE[index], linestyle=line_style)
        axes[1].plot(ties, moisture[:, node], color=COLOR_SEQUENCE[index], linestyle=line_style)
    axes[0].set(ylabel="Temperature (degC)", title="Responses at sled radii")
    axes[0].legend(ncol=5, loc="uer center")
    axes[1].set(xlabel="Time (s)", ylabel="Moisture (kg/kg)")
    export_plot(fig, "result_q1_selected_positions", size=(7.2, 4.8))



def write_contract() -> None:
    rs = [
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
        writer.writerows(rs)


# 检查合理性指标
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


# 组织当前脚本的完整执行流程并返回运行状态
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
