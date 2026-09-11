#!/usr/bin/env python3
"""ASCII-path entry point for the reproducibility command."""

from pathlib import Path
import runpy


if __name__ == "__main__":
    project_root = Path(__file__).resolve().parent
    runpy.run_path(str(project_root / "q1_cn_comparison.py"), run_name="__main__")
