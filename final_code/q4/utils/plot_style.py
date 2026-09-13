"""第三问中文论文图的统一样式。"""

from __future__ import annotations

import matplotlib as mpl
from cycler import cycler


COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7",
          "#E69F00", "#56B4E9", "#4D4D4D"]

# 应用本问题论文图表的统一绘图样式。
def apply_style() -> None:
    mpl.rcParams.update({
        "figure.facecolor": "white", "axes.facecolor": "white",
        "font.family": "sans-serif",
        "font.sans-serif": ["Microsoft YaHei", "SimHei", "Arial", "DejaVu Sans"],
        "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 9,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.unicode_minus": False, "legend.frameon": False,
        "axes.prop_cycle": cycler(color=COLORS),
        "svg.fonttype": "none", "pdf.fonttype": 42,
        "savefig.facecolor": "white", "savefig.dpi": 400,
    })
