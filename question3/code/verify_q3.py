#!/usr/bin/env python3
"""第三问产物的一键一致性验证。"""

from __future__ import annotations

import json
import zipfile
from pathlib import Path

import numpy as np

from xlsx_helpers import xlsx_rows


ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"


def main() -> int:
    formal = json.loads((RESULTS / "optimized400_dt05_summary.json").read_text(encoding="utf-8"))
    refined = json.loads((RESULTS / "optimized400_dt025_summary.json").read_text(encoding="utf-8"))
    coarse = json.loads((RESULTS / "optimized200_dt05_summary.json").read_text(encoding="utf-8"))
    z = np.load(RESULTS / "optimized400_dt05_fields.npz")
    event = z["event_moisture_kg_kg"]
    event_time = float(z["event_time_s"])
    regular_time = z["time_s"]
    checks = {
        "finite_fields": bool(np.isfinite(z["moisture_kg_kg"]).all()
                              and np.isfinite(z["temperature_C"]).all()),
        "event_all_points_below_threshold": bool(np.max(event) < 0.15),
        "event_argmax_is_center": bool(int(np.argmax(event)) == 0),
        "stored_radial_monotone": bool(np.max(np.diff(z["moisture_kg_kg"], axis=1)) < 1e-9),
        "regular_output_interval_60s": bool(np.allclose(np.diff(regular_time), 60.0)),
        "time_refinement_difference_below_1s": bool(abs(
            formal["event"]["time_s"] - refined["event"]["time_s"]) < 1.0),
        "space_refinement_relative_below_1e_3": bool(abs(
            coarse["event"]["time_s"] - event_time) / event_time < 1e-3),
        "no_unconverged_steps": bool(formal["diagnostics"]["unconverged_steps"] == 0),
        "no_negative_or_nonfinite_steps": bool(
            formal["diagnostics"]["negative_moisture_steps"] == 0
            and formal["diagnostics"]["nonfinite_steps"] == 0),
    }
    expected = {"result3.xlsx": (3450, 22), "q3_full_result.xlsx": (3450, 22),
                "q3_result_table.xlsx": (11, 6), "q3_final_summary.xlsx": (7, 8)}
    xlsx_details = {}
    for name, shape in expected.items():
        path = RESULTS / name
        with zipfile.ZipFile(path) as archive:
            zip_ok = archive.testzip() is None
        rows = xlsx_rows(path)
        actual = (len(rows), len(rows[0]), len(rows[-1]))
        checks[f"xlsx_{name}_valid"] = bool(zip_ok and actual == (shape[0], shape[1], shape[1]))
        xlsx_details[name] = actual
    for category in ("raw", "process", "result"):
        png = list(FIGURES.glob(f"{category}_q3_*.png"))
        svg = list(FIGURES.glob(f"{category}_q3_*.svg"))
        checks[f"figures_{category}_three_png_svg"] = len(png) >= 3 and len(svg) >= 3
    report = {"passed": bool(all(checks.values())), "checks": checks,
              "event_time_s": event_time, "event_time_h": event_time / 3600.0,
              "xlsx_shapes": xlsx_details}
    (RESULTS / "q3_verification_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
