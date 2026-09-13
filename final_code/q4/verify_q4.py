#!/usr/bin/env python3
"""问题四正式产物的自动验收。任一关键检查失败时返回非零状态。"""

from __future__ import annotations

import json
import zipfile
from pathlib import Path

import numpy as np

from xlsx_helpers import xlsx_rows


# 保留 subst 盘符，兼容旧版 Windows Python 的中文路径。
ROOT = Path(__file__).absolute().parent.parent
RESULTS = ROOT / "results"


# 把单项验收条件转换为结构化检查结果。
def check(condition: bool, detail: str) -> dict[str, object]:
    return {"passed": bool(condition), "detail": detail}


# 组织当前脚本的完整执行流程并返回运行状态。
def main() -> int:
    summary = json.loads((RESULTS / "q4_summary.json").read_text(encoding="utf-8"))
    run_summary = json.loads(
        (RESULTS / "shrink200_dt05_summary.json").read_text(encoding="utf-8"))
    z = dict(np.load(RESULTS / "shrink200_dt05_fields.npz"))
    checks: dict[str, dict[str, object]] = {}

    fields_finite = all(np.isfinite(z[name]).all() for name in (
        "time_s", "xi", "radius_time_m", "temperature_C", "moisture_kg_kg",
        "event_temperature_C", "event_moisture_kg_kg"))
    checks["数值场无NaN或Inf"] = check(fields_finite, "温度、水分、半径与时空坐标均为有限数")
    checks["水分非负"] = check(float(np.min(z["moisture_kg_kg"])) >= 0.0,
                            f"最小值={np.min(z['moisture_kg_kg']):.8g} kg/kg")
    checks["半径始终为正且不增加"] = check(
        bool(np.all(z["radius_time_m"] > 0.0) and
             np.all(np.diff(z["radius_time_m"]) <= 1e-12)),
        f"R: {z['radius_time_m'][0]*100:.4f} -> {float(z['event_radius_m'])*100:.4f} cm")

    event_field = z["event_moisture_kg_kg"]
    event_time = float(z["event_time_s"])
    checks["全域阈值达标"] = check(
        float(np.max(event_field)) <= 0.15 + 2e-10,
        f"t*={event_time:.6f} s，max(C)={np.max(event_field):.12f} kg/kg")
    checks["最湿点位于中心"] = check(
        int(np.argmax(event_field)) == 0 and
        run_summary["diagnostics"]["max_center_vs_global_gap"] < 1e-10,
        f"事件时刻argmax节点={int(np.argmax(event_field))}")
    checks["径向剖面单调"] = check(
        run_summary["diagnostics"]["radial_monotonicity_violations"] == 0,
        "全部存储时刻均未出现由中心向表面的反向增湿")
    checks["非线性迭代全部收敛"] = check(
        run_summary["diagnostics"]["unconverged_steps"] == 0,
        f"Picard最大迭代次数={run_summary['diagnostics']['max_picard_iterations']}")
    checks["离散方程残差"] = check(
        run_summary["diagnostics"]["max_temperature_algebraic_residual"] < 1e-8 and
        run_summary["diagnostics"]["max_moisture_algebraic_residual"] < 1e-10,
        "温度最大残差={:.3e}，水分最大残差={:.3e}".format(
            run_summary["diagnostics"]["max_temperature_algebraic_residual"],
            run_summary["diagnostics"]["max_moisture_algebraic_residual"]))
    checks["空间网格收敛"] = check(
        summary["spatial_relative_difference"] < 1e-3,
        "N=200与N=400的达标时刻相对差={:.3e}".format(
            summary["spatial_relative_difference"]))
    checks["时间步长收敛"] = check(
        summary["temporal_dt05_to_dt025_difference_s"] < 1.0,
        "dt=0.5 s与0.25 s的达标时刻差={:.6f} s".format(
            summary["temporal_dt05_to_dt025_difference_s"]))
    seam = summary["four_hour_continuity"]
    checks["4小时状态连续"] = check(
        seam["max_C_change_14340_to_14400"] < 0.01 and
        seam["max_C_change_14400_to_14460"] < 0.01 and
        seam["max_T_change_14340_to_14400_C"] < 0.1 and
        seam["max_T_change_14400_to_14460_C"] < 0.1,
        "边界条件切换前后未重置内部温湿度场")
    checks["收缩效应方向合理"] = check(
        summary["t4_shrink_s"] < summary["t4_fixed_s"],
        "收缩模型{:.4f} h，固定半径对照{:.4f} h".format(
            summary["t4_shrink_h"], summary["t4_fixed_h"]))

    result_rows = xlsx_rows(RESULTS / "result4.xlsx")
    table_rows = xlsx_rows(RESULTS / "q4_result_table.xlsx")
    checks["result4表格尺寸"] = check(
        len(result_rows) == 3067 and len(result_rows[0]) == 22,
        f"实际{len(result_rows)}行×{len(result_rows[0])}列（含表头）")
    checks["论文表格尺寸"] = check(
        len(table_rows) == 10 and len(table_rows[0]) == 6,
        f"实际{len(table_rows)}行×{len(table_rows[0])}列（含表头）")
    result_times = np.asarray([row[0] for row in result_rows[1:]], dtype=float)
    checks["常规输出间隔60秒"] = check(
        bool(np.allclose(np.diff(result_times[:-1]), 60.0)),
        f"常规行数={len(result_times)-1}，末行另列临界时刻{result_times[-1]:.6f} s")

    headers = result_rows[0]
    outside_blank = True
    surface_consistent = True
    for row in result_rows[1:]:
        t = float(row[0])
        radius_cm = float(np.interp(t, z["time_s"], z["radius_time_m"])) * 100.0
        for j, position in enumerate(np.arange(0.0, 2.0, 0.1), start=1):
            value = row[j] if j < len(row) else None
            if position > radius_cm + 1e-8 and value is not None:
                outside_blank = False
        if len(row) < 22 or row[21] is None:
            surface_consistent = False
    checks["移动表面外固定点留空"] = check(outside_blank, "所有超出当前半径的固定物理位置均为空")
    checks["独立药材表面列完整"] = check(surface_consistent, "每个输出时刻均含移动表面水分值")
    checks["表头固定位置正确"] = check(
        len(headers) == 22 and np.allclose(np.asarray(headers[1:21], float),
                                           np.arange(0.0, 2.0, 0.1)),
        "固定距离为0.0--1.9 cm，步长0.1 cm，另设药材表面列")

    xlsx_ok = True
    xlsx_detail: list[str] = []
    for name in ("result4.xlsx", "q4_result_table.xlsx", "q4_convergence.xlsx"):
        with zipfile.ZipFile(RESULTS / name) as archive:
            bad = archive.testzip()
        xlsx_ok &= bad is None
        xlsx_detail.append(f"{name}:{'OK' if bad is None else bad}")
    checks["XLSX文件结构完整"] = check(xlsx_ok, "；".join(xlsx_detail))

    passed = all(item["passed"] for item in checks.values())
    report = {
        "status": "PASS" if passed else "FAIL",
        "checks_passed": sum(bool(item["passed"]) for item in checks.values()),
        "checks_total": len(checks),
        "checks": checks,
    }
    (RESULTS / "q4_verification_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
