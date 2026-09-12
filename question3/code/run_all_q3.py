#!/usr/bin/env python3
"""从原始附件开始，一键复现问题三的求解、表格、图件与验证。"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
CODE = ROOT / "code"

RUNS = [
    (200, 2.00, 0.50, "graded200_dt05"),
    (400, 2.00, 0.50, "graded400_dt05"),
    (400, 2.00, 0.25, "graded400_dt025"),
    (200, 1.75, 0.50, "optimized200_dt05"),
    (400, 1.75, 0.50, "optimized400_dt05"),
    (400, 1.75, 0.25, "optimized400_dt025"),
]


def run(*args: str) -> None:
    command = [sys.executable, *args]
    print("\n>", " ".join(command), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


def main() -> int:
    print("问题三全流程复现开始；六组长时程求解预计耗时较长。", flush=True)
    for intervals, power, dt_s, label in RUNS:
        run(str(CODE / "q3_solver.py"), "--intervals", str(intervals),
            "--mesh-power", str(power), "--dt", str(dt_s), "--label", label)

    # 收敛性组仅需摘要；保留正式组全场数据供表格和绘图使用。
    for _, _, _, label in RUNS:
        if label != "optimized400_dt05":
            path = ROOT / "results" / f"{label}_fields.npz"
            if path.exists():
                path.unlink()

    run(str(CODE / "export_q3_results.py"))
    run(str(CODE / "make_q3_figures.py"))
    # 图件生成后再次导出，使复现清单记录最终图件哈希。
    run(str(CODE / "export_q3_results.py"))
    run(str(CODE / "verify_q3.py"))
    print("问题三全流程复现及验证完成。", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
