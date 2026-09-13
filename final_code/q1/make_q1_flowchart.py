#!/usr/bin/env python3
"""从问题一正式代码与结果重新生成论文流程图"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Polygon
from PIL import Image


ROOT = Path(__file__).absolute().parent.parent
FIGURES = ROOT / "figures"
QA = FIGURES / "_qa"

NAVY = "#18324B"
BLUE = "#3977A6"
ORANGE = "#C66B3D"
TEAL = "#2A7F78"
GRAY = "#687684"
PALE_BLUE = "#EFF5F9"
PALE_ORANGE = "#FBF2EC"
PALE_TEAL = "#EDF7F5"
PALE_GRAY = "#F5F7F8"


# 配置中文字体、配色和论文图形的统一样式
def setup_style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Microsoft YaHei", "SimHei", "Arial", "DejaVu Sans"],
        "axes.unicode_minus": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    })


# 在流程图中绘制一个包含标题和说明的步骤卡片
def card(ax, center, width, height, title, lines, *, edge=NAVY,
         face="white", number=None, title_color=None, fontsize=7.5):
    x, y = center
    box = FancyBboxPatch(
        (x - width / 2, y - height / 2), width, height,
        boxstyle="round,pad=0.03,rounding_size=0.08",
        linewidth=1.05, edgecolor=edge, facecolor=face, zorder=2)
    ax.add_patch(box)
    ax.text(x, y + height * 0.25, title, ha="center", va="center",
            fontsize=9, fontweight="bold", color=title_color or edge, zorder=3)
    ax.text(x, y - height * 0.12, "\n".join(lines), ha="center", va="center",
            fontsize=fontsize, color="#17212B", linespacing=1.35, zorder=3)
    if number is not None:
        badge = Circle((x - width / 2 + 0.24, y + height / 2 - 0.23), 0.18,
                       facecolor=edge, edgecolor="none", zorder=4)
        ax.add_patch(badge)
        ax.text(badge.center[0], badge.center[1], str(number), color="white",
                ha="center", va="center", fontsize=6.7, fontweight="bold", zorder=5)
    return box


# 在流程图两个步骤之间绘制有向连接箭头
def arrow(ax, start, end, *, color=NAVY, connectionstyle="arc3", lw=1.05):
    item = FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=9,
                           linewidth=lw, color=color,
                           connectionstyle=connectionstyle, zorder=1)
    ax.add_patch(item)
    return item


# 组织当前脚本的完整执行流程并返回运行状态
def main() -> int:
    setup_style()
    FIGURES.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 6.0), layout="constrained")
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 12.2)
    ax.axis("off")

    # 第一层：真实输入、插值选择与物理抽象
    card(ax, (2.25, 10.65), 3.55, 1.65, "附件 1：环境边界",
         ["0–1800 s，每 60 s 观测", "烘房温度 Ta(t)、水分浓度 Cb(t)"],
         edge=BLUE, face=PALE_BLUE, number=1, fontsize=7.0)
    card(ax, (8.0, 10.65), 4.45, 1.65, "边界曲线重构",
         ["Linear / PCHIP / Spline / Akima", "留一检验最优且无过冲 → 分段线性插值"],
         edge=BLUE, face="white", number=2, fontsize=7.0)
    card(ax, (13.75, 10.65), 3.55, 1.65, "圆柱径向抽象",
         ["长圆柱、轴对称、固定 R=2 cm", "忽略端面效应，仅保留径向传递"],
         edge=NAVY, face=PALE_GRAY, number=3, fontsize=7.0)
    arrow(ax, (4.03, 10.65), (5.75, 10.65), color=BLUE)
    arrow(ax, (10.23, 10.65), (11.97, 10.65), color=NAVY)

    # 第二层：两个物理场并行建模，强调方向与边界
    card(ax, (4.55, 8.10), 5.45, 1.95, "温度场 T(r,t)",
         ["Fourier 径向导热；ρ、cp、k 取常数", "中心零通量；表面对流换热 h=25 W/(m²·K)"],
         edge=ORANGE, face=PALE_ORANGE, number="4T", title_color=ORANGE, fontsize=7.1)
    card(ax, (11.45, 8.10), 5.45, 1.95, "水分场 C(r,t)",
         [r"Fick 径向扩散；$D(C)=7\times10^{-9}\exp(-0.89/C)$",
          r"中心零通量；表面对流传质 $h_m=8\times10^{-7}\ \mathrm{m/s}$"],
         edge=TEAL, face=PALE_TEAL, number="4C", title_color=TEAL, fontsize=7.1)
    arrow(ax, (13.75, 9.82), (6.3, 9.13), color=ORANGE,
          connectionstyle="arc3,rad=0.08")
    arrow(ax, (13.75, 9.82), (11.45, 9.08), color=TEAL,
          connectionstyle="arc3,rad=-0.12")

    # 第三层：统一离散与水分非线性闭环
    card(ax, (8.0, 5.35), 9.0, 2.25, "守恒型数值求解",
         ["节点中心有限体积：中心/表面通量直接入方程",
          "Crank–Nicolson：N=800，Δt=0.125 s",
          r"温度三对角直接求解；水分 Picard 更新 $D(C)$，容差 $10^{-11}\ \mathrm{kg/kg}$"],
         edge=NAVY, face=PALE_BLUE, number=5, fontsize=7.45)
    arrow(ax, (4.55, 7.12), (6.1, 6.48), color=ORANGE,
          connectionstyle="angle3,angleA=270,angleB=180")
    arrow(ax, (11.45, 7.12), (9.9, 6.48), color=TEAL,
          connectionstyle="angle3,angleA=270,angleB=0")
    # 第四层：五条证据先检查，再进入带返工分支的质量门
    checks = [
        ("网格/步长", "接近二阶趋势"),
        ("代数/通量", "残差闭合"),
        ("解析退化", "Bessel 对照"),
        ("物理约束", "范围与单调性"),
        ("跨实现", "70 点复核"),
    ]
    x_positions = [2.05, 5.0, 8.0, 11.0, 13.95]
    for x, (top, bottom) in zip(x_positions, checks):
        card(ax, (x, 3.15), 2.55, 1.10, top, [bottom],
             edge=GRAY, face=PALE_GRAY, fontsize=6.7)
    ax.plot([2.05, 13.95], [3.92, 3.92], color=GRAY, lw=0.72, zorder=1)
    arrow(ax, (8.0, 4.22), (8.0, 3.92), color=NAVY)
    for x in x_positions:
        arrow(ax, (x, 3.92), (x, 3.70), color=GRAY, lw=0.72)
    ax.plot([2.05, 13.95], [2.45, 2.45], color=GRAY, lw=0.72, zorder=1)
    for x in x_positions:
        ax.plot([x, x], [2.60, 2.45], color=GRAY, lw=0.72, zorder=1)

    diamond = Polygon([(8.0, 2.35), (9.45, 1.82), (8.0, 1.29), (6.55, 1.82)],
                      closed=True, facecolor="white", edgecolor=NAVY,
                      linewidth=1.05, zorder=2)
    ax.add_patch(diamond)
    ax.text(8.0, 1.82, "五层检验\n是否全部通过？", ha="center", va="center",
            fontsize=7.7, fontweight="bold", color=NAVY, linespacing=1.2, zorder=3)
    arrow(ax, (8.0, 2.45), (8.0, 2.35), color=NAVY)

    # “否”分支回到求解器，形成完整的检查—返工闭环
    ax.plot([9.45, 15.35, 15.35], [1.82, 1.82, 5.35],
            color=ORANGE, lw=0.9, zorder=1)
    arrow(ax, (15.35, 5.35), (12.52, 5.35), color=ORANGE, lw=0.9)
    ax.text(14.25, 3.85, "否：调整网格/步长\n或检查离散实现", ha="center", va="center",
            fontsize=6.7, color=ORANGE, linespacing=1.25,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.2})

    # 输出采用底部横条，展示题目规定的最终交付
    output = FancyBboxPatch((2.15, 0.10), 11.7, 0.76,
                            boxstyle="round,pad=0.03,rounding_size=0.08",
                            linewidth=1.05, edgecolor=BLUE,
                            facecolor=PALE_BLUE, zorder=2)
    ax.add_patch(output)
    ax.text(8.0, 0.48,
            "输出 T(r,t)、C(r,t)：7 个指定时刻 × 5 个径向位置 + 每 1 s、每 0.1 cm 完整场 + result1.xlsx",
            ha="center", va="center", fontsize=7.25, color=NAVY,
            fontweight="bold", zorder=3)
    arrow(ax, (8.0, 1.29), (8.0, 0.86), color=BLUE)
    ax.text(8.32, 1.08, "是", ha="left", va="center", fontsize=6.8, color=BLUE)

    basename = FIGURES / "process_q1_analysis_framework_award"
    fig.savefig(basename.with_suffix(".png"), dpi=400)
    fig.savefig(basename.with_suffix(".svg"))
    fig.savefig(basename.with_suffix(".pdf"))
    with Image.open(basename.with_suffix(".png")) as source:
        source.convert("L").save(QA / "process_q1_analysis_framework_award_grayscale.png",
                                 dpi=(400, 400))
    plt.close(fig)
    print(f"已重新生成第一问流程图：{basename}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
