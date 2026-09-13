#!/usr/bin/env python3
"""将问题四正式计算导出为题目表格、诊断、MAT和复现清单。"""

from __future__ import annotations

import csv
import hashlib
import json
import platform
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import matplotlib
import numpy as np
import scipy
from scipy.io import savemat

import xlsx_helpers
from q4_solver import load_environment, load_radius_data


# 保留 subst 盘符，兼容旧版 Windows Python 的中文路径。
ROOT = Path(__file__).absolute().parent.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
FORMAL = "shrink200_dt05"
CONTROL = "fixed200_dt05"
FULL_RADII_CM = np.arange(0.0, 2.0, 0.1)  # 0--1.9 cm，另设移动表面列
TABLE_RADII_CM = np.arange(0.0, 2.0, 0.5)  # 0、0.5、1.0、1.5 cm


# 计算文件的 SHA-256 哈希值结果复现校验。
def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


# 读取指定标签对应的模拟摘要和结果。
def load_run(label: str) -> tuple[dict, dict[str, np.ndarray]]:
    summary = json.loads((RESULTS / f"{label}_summary.json").read_text(encoding="utf-8"))
    return summary, dict(np.load(RESULTS / f"{label}_fields.npz"))


# 将参考坐标场转换并插值到物理半径位置。
def interpolate_physical(xi: np.ndarray, field: np.ndarray, radius_m: float,
                         positions_cm: np.ndarray) -> list[float | None]:
    values: list[float | None] = []
    for position_cm in positions_cm:
        position_m = float(position_cm) / 100.0
        if position_m > radius_m + 1.0e-12:
            values.append(None)
        else:
            values.append(float(np.interp(position_m / radius_m, xi, field)))
    return values


# 将结构化数据行写入指定 CSV 文件。
def write_csv(path: Path, header: list[object], rows: list[list[object]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)


# 把表格数据编码为可写入 XLSX 的工作表 XML。
def _sheet_xml(rows: list[list[object]], numeric_style: int) -> bytes:
    xml_rows: list[str] = []
    width = max(len(row) for row in rows)
    for row_number, row in enumerate(rows, 1):
        cells: list[str] = []
        for column, value in enumerate(row, 1):
            if value is None:
                continue
            reference = f"{xlsx_helpers.excel_column(column)}{row_number}"
            if isinstance(value, (int, float, np.integer, np.floating)):
                cells.append(f'<c r="{reference}" s="{numeric_style}"><v>{float(value):.12g}</v></c>')
            else:
                escaped = (str(value).replace("&", "&amp;").replace("<", "&lt;")
                           .replace(">", "&gt;"))
                cells.append(f'<c r="{reference}" t="inlineStr"><is><t>{escaped}</t></is></c>')
        xml_rows.append(f'<row r="{row_number}">{"".join(cells)}</row>')
    last_col = xlsx_helpers.excel_column(width)
    xml = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
           '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
           f'<dimension ref="A1:{last_col}{len(rows)}"/>'
           '<sheetViews><sheetView workbookViewId="0"/></sheetViews>'
           '<sheetFormatPr defaultRowHeight="15"/><cols>'
           '<col min="1" max="1" width="24" customWidth="1"/>'
           f'<col min="2" max="{width}" width="13" customWidth="1"/></cols>'
           f'<sheetData>{"".join(xml_rows)}</sheetData>'
           '<pageMargins left="0.7" right="0.7" top="0.75" bottom="0.75" '
           'header="0.3" footer="0.3"/></worksheet>')
    return xml.encode("utf-8")


# 按原题模板样式生成结果 XLSX 工作簿。
def write_template_xlsx(output: Path, rows: list[list[object]], sheet_name: str) -> None:
    template = DATA / "result4_template.xlsx"
    temporary = output.with_suffix(".tmp.xlsx")
    with zipfile.ZipFile(template) as source, zipfile.ZipFile(
            temporary, "w", zipfile.ZIP_DEFLATED) as target:
        styles, numeric_style = xlsx_helpers.add_four_decimal_style(
            source.read("xl/styles.xml"))
        workbook = ET.fromstring(source.read("xl/workbook.xml"))
        sheets = workbook.find(f"{{{xlsx_helpers.NS['m']}}}sheets")
        if sheets is not None and len(sheets):
            sheets[0].attrib["name"] = sheet_name
        replacements = {
            "xl/styles.xml": styles,
            "xl/workbook.xml": ET.tostring(workbook, encoding="utf-8", xml_declaration=True),
            "xl/worksheets/sheet1.xml": _sheet_xml(rows, numeric_style),
        }
        for item in source.infolist():
            target.writestr(item, replacements.get(item.filename, source.read(item.filename)))
    temporary.replace(output)


# 整理问题四完整时程结果的输出表格行。
def full_rows(z: dict[str, np.ndarray]) -> list[list[object]]:
    event_time = float(z["event_time_s"])
    rows: list[list[object]] = [[
        "时间/s\\到药材中心的距离/cm", *FULL_RADII_CM.tolist(), "药材表面"]]
    for time_s, radius_m, field in zip(z["time_s"], z["radius_time_m"],
                                       z["moisture_kg_kg"]):
        if time_s <= 0.0 or time_s >= event_time:
            continue
        values = interpolate_physical(z["xi"], field, float(radius_m), FULL_RADII_CM)
        rows.append([float(time_s), *values, float(field[-1])])
    event_radius = float(z["event_radius_m"])
    event_field = z["event_moisture_kg_kg"]
    rows.append([event_time,
                 *interpolate_physical(z["xi"], event_field, event_radius, FULL_RADII_CM),
                 float(event_field[-1])])
    return rows


# 整理问题四论文正文所需的精简结果表格行。
def paper_rows(z: dict[str, np.ndarray]) -> list[list[object]]:
    event_time = float(z["event_time_s"])
    rows: list[list[object]] = [[
        "时间/h", *[f"{x:g} cm" for x in TABLE_RADII_CM], "药材表面"]]
    for hour in np.arange(6.0, event_time / 3600.0, 6.0):
        idx = int(np.argmin(np.abs(z["time_s"] - hour * 3600.0)))
        radius_m = float(z["radius_time_m"][idx])
        field = z["moisture_kg_kg"][idx]
        values = interpolate_physical(z["xi"], field, radius_m, TABLE_RADII_CM)
        rows.append([float(hour), *values, float(field[-1])])
    event_radius = float(z["event_radius_m"])
    event_field = z["event_moisture_kg_kg"]
    rows.append([event_time / 3600.0,
                 *interpolate_physical(z["xi"], event_field, event_radius, TABLE_RADII_CM),
                 float(event_field[-1])])
    return rows


# 组织当前脚本的完整执行流程并返回运行状态。
def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    formal_summary, formal = load_run(FORMAL)
    control_summary, control = load_run(CONTROL)
    fine_space_summary, _ = load_run("shrink400_dt05")
    fine_time_summary, _ = load_run("shrink400_dt025")
    formal_time = float(formal["event_time_s"])
    control_time = float(control["event_time_s"])
    spatial_difference = abs(fine_space_summary["event"]["time_s"] - formal_time)
    temporal_difference = abs(
        fine_time_summary["event"]["time_s"] - fine_space_summary["event"]["time_s"])

    complete = full_rows(formal)
    write_template_xlsx(RESULTS / "result4.xlsx", complete, "水分浓度")
    table = paper_rows(formal)
    write_template_xlsx(RESULTS / "q4_result_table.xlsx", table, "表6")
    formatted = [[f"{float(row[0]):.4f}",
                  *[("" if value is None else f"{float(value):.4f}") for value in row[1:]]]
                 for row in table[1:]]
    write_csv(RESULTS / "q4_result_table.csv", table[0], formatted)

    convergence = [
        ["模型", "N", "加密指数", "时间步/s", "烘干时间/h", "与正式结果差/s", "径向反向增量步数"],
        ["收缩正式", 200, 1.75, 0.5, formal_time / 3600.0, 0.0,
         formal_summary["diagnostics"]["radial_monotonicity_violations"]],
        ["收缩空间细化", 400, 1.75, 0.5,
         fine_space_summary["event"]["time_s"] / 3600.0,
         fine_space_summary["event"]["time_s"] - formal_time,
         fine_space_summary["diagnostics"]["radial_monotonicity_violations"]],
        ["收缩时间细化", 400, 1.75, 0.25,
         fine_time_summary["event"]["time_s"] / 3600.0,
         fine_time_summary["event"]["time_s"] - formal_time,
         fine_time_summary["diagnostics"]["radial_monotonicity_violations"]],
        ["固定半径控制", 200, 1.75, 0.5, control_time / 3600.0,
         control_time - formal_time,
         control_summary["diagnostics"]["radial_monotonicity_violations"]],
    ]
    write_template_xlsx(RESULTS / "q4_convergence.xlsx", convergence, "收敛与对照")
    write_csv(RESULTS / "q4_convergence.csv", convergence[0], convergence[1:])

    seam_index = int(np.flatnonzero(np.isclose(formal["time_s"], 14400.0))[0])
    seam = {
        "max_C_change_14340_to_14400": float(np.max(np.abs(
            formal["moisture_kg_kg"][seam_index] - formal["moisture_kg_kg"][seam_index - 1]))),
        "max_C_change_14400_to_14460": float(np.max(np.abs(
            formal["moisture_kg_kg"][seam_index + 1] - formal["moisture_kg_kg"][seam_index]))),
        "max_T_change_14340_to_14400_C": float(np.max(np.abs(
            formal["temperature_C"][seam_index] - formal["temperature_C"][seam_index - 1]))),
        "max_T_change_14400_to_14460_C": float(np.max(np.abs(
            formal["temperature_C"][seam_index + 1] - formal["temperature_C"][seam_index]))),
    }
    summary = {
        "formal_grid_N": 200,
        "formal_mesh_power": 1.75,
        "formal_dt_s": 0.5,
        "t4_shrink_s": formal_time,
        "t4_shrink_h": formal_time / 3600.0,
        "strict_integer_second_s": int(np.floor(formal_time)) + 1,
        "t4_fixed_s": control_time,
        "t4_fixed_h": control_time / 3600.0,
        "shrink_time_reduction_h": (control_time - formal_time) / 3600.0,
        "shrink_relative_reduction_percent": (control_time - formal_time) / control_time * 100.0,
        "event_radius_cm": float(formal["event_radius_m"]) * 100.0,
        "event_max_C": formal_summary["event"]["max_C"],
        "event_surface_C": formal_summary["event"]["surface_C"],
        "event_average_C": formal_summary["event"]["average_C"],
        "argmax_always_center": formal_summary["diagnostics"]["max_center_vs_global_gap"] < 1.0e-10,
        "spatial_N200_to_N400_difference_s": spatial_difference,
        "spatial_relative_difference": spatial_difference / formal_time,
        "temporal_dt05_to_dt025_difference_s": temporal_difference,
        "temporal_relative_difference": temporal_difference / fine_space_summary["event"]["time_s"],
        "four_hour_continuity": seam,
        "formal_diagnostics": formal_summary["diagnostics"],
        "control_diagnostics": control_summary["diagnostics"],
    }
    (RESULTS / "q4_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    summary_rows = [["项目", "数值", "单位"],
                    ["正式网格N", 200, ""], ["正式时间步", 0.5, "s"],
                    ["收缩模型烘干时间", formal_time / 3600.0, "h"],
                    ["固定半径控制烘干时间", control_time / 3600.0, "h"],
                    ["收缩减少时间", (control_time - formal_time) / 3600.0, "h"],
                    ["相对减少率", (control_time - formal_time) / control_time * 100.0, "%"],
                    ["结束半径", float(formal["event_radius_m"]) * 100.0, "cm"],
                    ["空间细化时刻差", spatial_difference, "s"],
                    ["时间步细化时刻差", temporal_difference, "s"],
                    ["Picard平均迭代", formal_summary["diagnostics"]["mean_picard_iterations"], "次/步"],
                    ["Picard最大迭代", formal_summary["diagnostics"]["max_picard_iterations"], "次/步"],
                    ["未收敛时间步", formal_summary["diagnostics"]["unconverged_steps"], "步"]]
    write_csv(RESULTS / "q4_summary.csv", summary_rows[0], summary_rows[1:])
    lines = [f"{row[0]}: {row[1]} {row[2]}".rstrip() for row in summary_rows[1:]]
    lines.extend([
        f"全程最大含水率位于中心: {summary['argmax_always_center']}",
        f"负含水率步数: {formal_summary['diagnostics']['negative_moisture_steps']}",
        f"NaN/Inf步数: {formal_summary['diagnostics']['nonfinite_steps']}",
        "说明: reference balance仅为当前简化参考域方程的离散闭合检查，不宣称完整物理总质量守恒。",
    ])
    (RESULTS / "diagnostics.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    radius_data = load_radius_data()
    environment = load_environment()
    write_csv(RESULTS / "q4_radius_data.csv", ["time_s", "radius_cm"],
              [[float(t), float(r * 100)] for t, r in radius_data])
    write_csv(RESULTS / "q4_environment.csv",
              ["time_s", "temperature_C", "moisture_kg_kg"], environment.tolist())
    savemat(RESULTS / "problem4_raw.mat", {
        "time_s": formal["time_s"], "xi": formal["xi"],
        "radius_time_m": formal["radius_time_m"],
        "temperature_C": formal["temperature_C"],
        "moisture_kg_kg": formal["moisture_kg_kg"],
        "event_time_s": formal["event_time_s"],
        "event_radius_m": formal["event_radius_m"],
        "event_temperature_C": formal["event_temperature_C"],
        "event_moisture_kg_kg": formal["event_moisture_kg_kg"],
    })

    manifest = {
        "formal_run": FORMAL,
        "random_seed": None,
        "parameters": {"N": 200, "mesh_power": 1.75, "dt_s": 0.5,
                       "theta": 0.5, "threshold_kg_kg": 0.15},
        "unique_command": "python code/run_all_q4.py",
        "commands": [
            "python code/q4_solver.py --intervals 200 --mesh-power 1.75 --dt 0.5 --max-hours 72 --label shrink200_dt05",
            "python code/q4_solver.py --intervals 400 --mesh-power 1.75 --dt 0.5 --max-hours 72 --label shrink400_dt05",
            "python code/q4_solver.py --intervals 400 --mesh-power 1.75 --dt 0.25 --max-hours 72 --label shrink400_dt025",
            "python code/q4_solver.py --intervals 200 --mesh-power 1.75 --dt 0.5 --max-hours 180 --fixed-radius --label fixed200_dt05",
            "python code/export_q4_results.py", "python code/make_q4_figures.py",
            "python code/export_q4_results.py", "python code/verify_q4.py"],
        "inputs": {path.name: sha256(path) for path in (
            DATA / "附件1.xlsx", DATA / "附件2.xlsx", DATA / "result4_template.xlsx")},
        "software": {"python": platform.python_version(), "numpy": np.__version__,
                     "scipy": scipy.__version__, "matplotlib": matplotlib.__version__},
        "outputs": {},
    }
    for directory in (RESULTS, ROOT / "figures"):
        if not directory.exists():
            continue
        for path in directory.iterdir():
            auxiliary = path.name.startswith(("rough_shrink", "smoke"))
            if path.is_file() and path.name != "复现清单.json" and not auxiliary:
                manifest["outputs"][str(path.relative_to(ROOT)).replace("\\", "/")] = sha256(path)
    (RESULTS / "复现清单.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
