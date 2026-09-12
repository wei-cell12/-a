#!/usr/bin/env python3
"""从正式 NPZ 生成第三问题目表格、XLSX、MAT 和诊断汇总。"""

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


ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
DATA = ROOT / "data"
FORMAL_LABEL = "optimized400_dt05"
RADII_CM = np.arange(0.0, 2.0001, 0.1)
TABLE_RADII_CM = np.arange(0.0, 2.0001, 0.5)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_run(label: str, with_arrays: bool = True) -> tuple[dict, dict | None]:
    summary = json.loads((RESULTS / f"{label}_summary.json").read_text(encoding="utf-8"))
    arrays = dict(np.load(RESULTS / f"{label}_fields.npz")) if with_arrays else None
    return summary, arrays


def interpolate_rows(source_r_m: np.ndarray, field: np.ndarray,
                     radii_cm: np.ndarray) -> np.ndarray:
    target_m = radii_cm / 100.0
    return np.asarray([np.interp(target_m, source_r_m, row) for row in field])


def write_csv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)


def _sheet_xml(rows: list[list[object]], numeric_style: int) -> bytes:
    xml_rows: list[str] = []
    width = max(len(row) for row in rows)
    for row_number, row in enumerate(rows, 1):
        cells: list[str] = []
        for column, value in enumerate(row, 1):
            reference = f"{xlsx_helpers.excel_column(column)}{row_number}"
            if isinstance(value, (int, float, np.integer, np.floating)):
                cells.append(f'<c r="{reference}" s="{numeric_style}"><v>{float(value):.12g}</v></c>')
            else:
                escaped = (str(value).replace("&", "&amp;").replace("<", "&lt;")
                           .replace(">", "&gt;"))
                cells.append(f'<c r="{reference}" t="inlineStr"><is><t>{escaped}</t></is></c>')
        xml_rows.append(f'<row r="{row_number}">{"".join(cells)}</row>')
    last = len(rows)
    last_col = xlsx_helpers.excel_column(width)
    xml = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
           '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
           f'<dimension ref="A1:{last_col}{last}"/><sheetViews><sheetView workbookViewId="0"/></sheetViews>'
           '<sheetFormatPr defaultRowHeight="15"/><cols><col min="1" max="1" width="22" customWidth="1"/>'
           f'<col min="2" max="{width}" width="13" customWidth="1"/></cols>'
           f'<sheetData>{"".join(xml_rows)}</sheetData><pageMargins left="0.7" right="0.7" '
           'top="0.75" bottom="0.75" header="0.3" footer="0.3"/></worksheet>')
    return xml.encode("utf-8")


def write_template_xlsx(output: Path, rows: list[list[object]], sheet_name: str) -> None:
    template = DATA / "result3_template.xlsx"
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


def convergence_rows() -> list[dict]:
    labels = [
        ("表面过密", "graded200_dt05"), ("表面过密", "graded400_dt05"),
        ("表面过密", "graded400_dt025"), ("优化加密", "optimized200_dt05"),
        ("优化加密", "optimized400_dt05"), ("优化加密", "optimized400_dt025"),
    ]
    records = []
    for mesh, label in labels:
        summary, _ = load_run(label, with_arrays=False)
        cfg, event, diagnostic = summary["config"], summary["event"], summary["diagnostics"]
        records.append({"mesh": mesh, "label": label, "N": cfg["intervals"],
                        "power": cfg["mesh_power"], "dt": cfg["dt_s"],
                        "time_s": event["time_s"], "time_h": event["time_s"] / 3600.0,
                        "violations": diagnostic["radial_monotonicity_violations"]})
    return records


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    summary, arrays = load_run(FORMAL_LABEL)
    assert arrays is not None
    event_time = float(arrays["event_time_s"])
    regular_mask = (arrays["time_s"] > 0.0) & (arrays["time_s"] < event_time)
    regular_times = arrays["time_s"][regular_mask]
    regular_field = arrays["moisture_kg_kg"][regular_mask]
    full_values = interpolate_rows(arrays["radius_m"], regular_field, RADII_CM)
    event_values = np.interp(RADII_CM / 100.0, arrays["radius_m"],
                             arrays["event_moisture_kg_kg"])
    xlsx_rows = [["时间/s\\到药材中心的距离/cm", *RADII_CM.tolist()]]
    xlsx_rows.extend([[float(tm), *row.tolist()] for tm, row in zip(regular_times, full_values)])
    xlsx_rows.append([event_time, *event_values.tolist()])
    write_template_xlsx(RESULTS / "result3.xlsx", xlsx_rows, "水分浓度")
    write_template_xlsx(RESULTS / "q3_full_result.xlsx", xlsx_rows, "水分浓度")

    key_times = np.arange(6.0, 60.0, 6.0) * 3600.0
    key_times = key_times[key_times < event_time]
    source_indices = [int(np.argmin(np.abs(arrays["time_s"] - tm))) for tm in key_times]
    table_values = interpolate_rows(arrays["radius_m"],
                                    arrays["moisture_kg_kg"][source_indices], TABLE_RADII_CM)
    event_table = np.interp(TABLE_RADII_CM / 100.0, arrays["radius_m"],
                            arrays["event_moisture_kg_kg"])
    table_rows = [["时间/h", *[f"{r:g} cm" for r in TABLE_RADII_CM]]]
    table_rows.extend([[float(tm / 3600.0), *row.tolist()]
                       for tm, row in zip(key_times, table_values)])
    table_rows.append([float(event_time / 3600.0), *event_table.tolist()])
    write_template_xlsx(RESULTS / "q3_result_table.xlsx", table_rows, "表5")
    paper_csv_rows = [[f"{float(row[0]):.4f}", *[f"{float(x):.4f}" for x in row[1:]]]
                      for row in table_rows[1:]]
    write_csv(RESULTS / "q3_result_table.csv", table_rows[0], paper_csv_rows)

    records = convergence_rows()
    reference = next(x for x in records if x["label"] == FORMAL_LABEL)
    convergence_table = [["网格类型", "N", "加密指数", "时间步/s", "烘干时间/h",
                          "与正式结果差/s", "相对差", "单调性破坏步数"]]
    for row in records:
        relative = abs(row["time_s"] - reference["time_s"]) / reference["time_s"]
        convergence_table.append([row["mesh"], row["N"], row["power"], row["dt"],
                                  row["time_h"], row["time_s"] - reference["time_s"],
                                  f"{relative:.6e}", row["violations"]])
    write_template_xlsx(RESULTS / "q3_final_summary.xlsx", convergence_table, "收敛性")
    write_csv(RESULTS / "q3_convergence.csv", convergence_table[0], convergence_table[1:])

    formal_diag = summary["diagnostics"]
    diagnostics = {
        "formal_label": FORMAL_LABEL,
        "critical_time_s": event_time,
        "critical_time_h": event_time / 3600.0,
        "strict_integer_second_s": int(np.floor(event_time)) + 1,
        "mean_threshold_time_s": summary["mean_threshold_time_s"],
        "mean_threshold_time_h": summary["mean_threshold_time_s"] / 3600.0,
        "premature_stop_gap_h": (event_time - summary["mean_threshold_time_s"]) / 3600.0,
        "event_max_C": summary["event"]["max_C"],
        "event_argmax_node": summary["event"]["argmax_node"],
        "event_surface_C": summary["event"]["surface_C"],
        "event_average_C": summary["event"]["average_C"],
        "spatial_N200_to_N400_difference_s": abs(
            next(x["time_s"] for x in records if x["label"] == "optimized200_dt05") - event_time),
        "temporal_dt05_to_dt025_difference_s": abs(
            next(x["time_s"] for x in records if x["label"] == "optimized400_dt025") - event_time),
        "formal_run": formal_diag,
    }
    (RESULTS / "q3_diagnostics.json").write_text(
        json.dumps(diagnostics, ensure_ascii=False, indent=2), encoding="utf-8")

    savemat(RESULTS / "problem3_raw.mat", {
        "time_s": arrays["time_s"], "radius_m": arrays["radius_m"],
        "temperature_C": arrays["temperature_C"],
        "moisture_kg_kg": arrays["moisture_kg_kg"],
        "event_time_s": event_time,
        "event_temperature_C": arrays["event_temperature_C"],
        "event_moisture_kg_kg": arrays["event_moisture_kg_kg"],
    })
    environment_rows = _environment_rows()
    write_csv(RESULTS / "q3_environment.csv", ["time_s", "temperature_C", "moisture_kg_kg"],
              environment_rows)
    manifest = {
        "formal_run": FORMAL_LABEL,
        "commands": [
            "python code/run_all_q3.py",
            "python code/q3_solver.py --intervals 200 --mesh-power 2 --dt 0.5 --label graded200_dt05",
            "python code/q3_solver.py --intervals 400 --mesh-power 2 --dt 0.5 --label graded400_dt05",
            "python code/q3_solver.py --intervals 400 --mesh-power 2 --dt 0.25 --label graded400_dt025",
            "python code/q3_solver.py --intervals 200 --mesh-power 1.75 --dt 0.5 --label optimized200_dt05",
            "python code/q3_solver.py --intervals 400 --mesh-power 1.75 --dt 0.5 --label optimized400_dt05",
            "python code/q3_solver.py --intervals 400 --mesh-power 1.75 --dt 0.25 --label optimized400_dt025",
            "python code/export_q3_results.py",
            "python code/make_q3_figures.py",
            "python code/export_q3_results.py",
            "python code/verify_q3.py",
        ],
        "inputs": {"附件1.xlsx": sha256(DATA / "附件1.xlsx"),
                   "result3_template.xlsx": sha256(DATA / "result3_template.xlsx")},
        "software": {"python": platform.python_version(), "numpy": np.__version__,
                     "scipy": scipy.__version__, "matplotlib": matplotlib.__version__},
        "outputs": {},
    }
    for path in RESULTS.iterdir():
        if path.is_file() and path.name != "q3_repro_manifest.json":
            manifest["outputs"][path.name] = sha256(path)
    figures = ROOT / "figures"
    if figures.exists():
        for path in figures.iterdir():
            if path.is_file() and path.suffix.lower() in {".png", ".svg", ".pdf"}:
                manifest["outputs"][f"figures/{path.name}"] = sha256(path)
    (RESULTS / "q3_repro_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(diagnostics, ensure_ascii=False, indent=2))
    return 0


def _environment_rows() -> list[list[float]]:
    from q3_solver import load_environment
    return load_environment().tolist()


if __name__ == "__main__":
    raise SystemExit(main())
