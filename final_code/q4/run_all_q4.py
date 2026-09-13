#!/usr/bin/env python3
"""从原始附件开始一键复现问题四的求解、表格、图片与自动验收。"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path


# 保留 subst 盘符，兼容旧版 Windows Python 的中文路径。
ROOT = Path(__file__).absolute().parent.parent
CODE = ROOT / "code"


# 调用指定 Python 子程序并在失败时终止总流程。
def run(*arguments: str) -> None:
    command = [sys.executable, *arguments]
    print("\n[运行]", " ".join(command), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


# 必要时在纯 ASCII 临时路径下重新启动计算流程。
def relaunch_on_ascii_drive() -> int | None:
    """旧版 Windows Python 遇到中文工程路径时，自动转到临时盘符。"""
    if os.name != "nt" or os.environ.get("Q4_ASCII_PATH") == "1":
        return None
    if all(ord(character) < 128 for character in str(ROOT)):
        return None
    drive = next((f"{letter}:" for letter in "ZYXWVUTSRQPONMLKJIHGFED"
                  if not Path(f"{letter}:\\").exists()), None)
    if drive is None:
        raise RuntimeError("没有可用盘符用于规避中文路径，请手工将 question4 复制到纯英文路径。")
    subprocess.run(["subst", drive, str(ROOT)], check=True)
    mapped_script = f"{drive}\\code\\run_all_q4.py"
    environment = os.environ.copy()
    environment["Q4_ASCII_PATH"] = "1"
    try:
        completed = subprocess.run([sys.executable, mapped_script, *sys.argv[1:]],
                                   cwd=f"{drive}\\", env=environment)
        return completed.returncode
    finally:
        subprocess.run(["subst", drive, "/d"], check=False)


# 组织当前脚本的完整执行流程并返回运行状态。
def main() -> int:
    relaunched = relaunch_on_ascii_drive()
    if relaunched is not None:
        return relaunched
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reuse-results", action="store_true",
                        help="跳过耗时求解，仅用已有正式场重新生成表、图和验收报告")
    args = parser.parse_args()
    solver = str(CODE / "q4_solver.py")
    if not args.reuse_results:
        run(solver, "--intervals", "200", "--mesh-power", "1.75", "--dt", "0.5",
            "--max-hours", "72", "--label", "shrink200_dt05", "--quiet")
        run(solver, "--intervals", "400", "--mesh-power", "1.75", "--dt", "0.5",
            "--max-hours", "72", "--label", "shrink400_dt05", "--quiet")
        run(solver, "--intervals", "400", "--mesh-power", "1.75", "--dt", "0.25",
            "--max-hours", "72", "--label", "shrink400_dt025", "--quiet")
        run(solver, "--intervals", "200", "--mesh-power", "1.75", "--dt", "0.5",
            "--max-hours", "180", "--fixed-radius", "--label", "fixed200_dt05", "--quiet")
    run(str(CODE / "export_q4_results.py"))
    run(str(CODE / "make_q4_figures.py"))
    run(str(CODE / "export_q4_results.py"))
    run(str(CODE / "verify_q4.py"))
    # 最后更新一次哈希清单，使验收报告也进入复现记录。
    run(str(CODE / "export_q4_results.py"))
    print("\n问题四全部产物已生成并通过自动验收。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
