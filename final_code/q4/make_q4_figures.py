#!/usr/bin/env python3
"""生成问题四 raw/process/result 三类中文论文候选图。"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
import numpy as np
from PIL import Image

from q4_solver import diffusivity, load_environment, load_radius_data
from utils.plot_style import COLORS, apply_style


# 保留 subst 盘符，兼容旧版 Windows Python 的中文路径。
ROOT = Path(__file__).absolute().parent.parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
QA = FIGURES / "_qa"
FORMAL = "shrink200_dt05"
CONTROL = "fixed200_dt05"
MOISTURE_CMAP = LinearSegmentedColormap.from_list(
    "q4_blue_yellow_red", ["#1646A0", "#35A7C8", "#F2D45C", "#C7392F"], N=256)


# 按论文尺寸和分辨率导出当前图形。
def export(fig: plt.Figure, name: str, size: tuple[float, float] = (6.8, 4.2)) -> None:
    fig.set_size_inches(*size)
    fig.savefig(FIGURES / f"{name}.png", dpi=400)
    fig.savefig(FIGURES / f"{name}.svg")
    fig.savefig(FIGURES / f"{name}.pdf")
    with Image.open(FIGURES / f"{name}.png") as source:
        source.convert("L").save(QA / f"{name}_grayscale.png")
    plt.close(fig)


# 读取指定工况的模拟摘要与场变量数组。
def load(label: str) -> tuple[dict, dict[str, np.ndarray]]:
    summary = json.loads((RESULTS / f"{label}_summary.json").read_text(encoding="utf-8"))
    return summary, dict(np.load(RESULTS / f"{label}_fields.npz"))


# 组织当前脚本的完整执行流程并返回运行状态。
def main() -> int:
    FIGURES.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    apply_style()
    summary = json.loads((RESULTS / "q4_summary.json").read_text(encoding="utf-8"))
    formal_summary, z = load(FORMAL)
    _, fixed = load(CONTROL)
    radius_data = load_radius_data()
    environment = load_environment()
    hours = z["time_s"] / 3600.0
    fixed_hours = fixed["time_s"] / 3600.0
    cmax = np.max(z["moisture_kg_kg"], axis=1)
    cavg = z["moisture_kg_kg"] @ (z["area_weights"] / 0.5)
    fixed_cmax = np.max(fixed["moisture_kg_kg"], axis=1)
    event_h = float(z["event_time_s"]) / 3600.0
    fixed_event_h = float(fixed["event_time_s"]) / 3600.0

    # raw 1/4：附件2原始点与分段线性插值。
    dense_t = np.linspace(radius_data[0, 0], radius_data[-1, 0], 2000)
    dense_r = np.interp(dense_t, radius_data[:, 0], radius_data[:, 1]) * 100.0
    fig, ax = plt.subplots(layout="constrained")
    ax.plot(dense_t / 3600.0, dense_r, color=COLORS[0], lw=1.2, label="分段线性插值")
    ax.scatter(radius_data[:, 0] / 3600.0, radius_data[:, 1] * 100.0,
               s=9, facecolor="white", edgecolor=COLORS[0], linewidth=0.55,
               label="附件2观测点", zorder=3)
    ax.axvline(event_h, color=COLORS[1], ls="--", lw=0.9,
               label=f"烘干结束 {event_h:.4f} h")
    ax.set(xlabel="时间（h）", ylabel="药材半径（cm）", xlim=(0, 72), ylim=(1.15, 2.05))
    ax.legend(loc="upper right")
    export(fig, "raw_q4_radius_interpolation")

    # raw 2/4：前4 h环境边界，使用双面板避免双Y轴。
    fig, axes = plt.subplots(2, 1, sharex=True, layout="constrained")
    axes[0].plot(environment[:, 0] / 3600.0, environment[:, 1], color=COLORS[1], lw=1.1)
    axes[0].set(ylabel="烘房温度（℃）", xlim=(0, 4))
    axes[1].plot(environment[:, 0] / 3600.0, environment[:, 2], color=COLORS[0], lw=1.1)
    axes[1].set(xlabel="时间（h）", ylabel="环境水分浓度（kg/kg）", xlim=(0, 4))
    export(fig, "raw_q4_environment", (6.8, 5.0))

    # raw 3/4：每半小时收缩量的时序与分布。
    shrink_mm = -np.diff(radius_data[:, 1]) * 1000.0
    midpoint_h = 0.5 * (radius_data[1:, 0] + radius_data[:-1, 0]) / 3600.0
    fig, axes = plt.subplots(1, 2, layout="constrained")
    axes[0].plot(midpoint_h, shrink_mm, color=COLORS[2], lw=0.9)
    axes[0].scatter(midpoint_h[::5], shrink_mm[::5], s=8, color=COLORS[2])
    axes[0].set(xlabel="时间（h）", ylabel="每0.5 h半径减少量（mm）", xlim=(0, 72))
    axes[1].hist(shrink_mm, bins=18, color=COLORS[2], alpha=0.75,
                 edgecolor="white", linewidth=0.5)
    axes[1].set(xlabel="每0.5 h半径减少量（mm）", ylabel="区间数")
    export(fig, "raw_q4_shrinkage_distribution", (7.2, 3.2))

    # raw 4/4：附录4扩散系数的温湿关系。
    c_axis = np.linspace(0.05, 2.55, 220)
    t_axis = np.linspace(28.0, 50.0, 150)
    cc, tt = np.meshgrid(c_axis, t_axis)
    log_d = np.log10(diffusivity(tt, cc))
    fig, ax = plt.subplots(layout="constrained")
    levels = np.linspace(log_d.min(), log_d.max(), 65)
    mesh = ax.contourf(c_axis, t_axis, log_d, levels=levels, cmap="cividis")
    cb = fig.colorbar(mesh, ax=ax)
    if cb.solids is not None:
        cb.solids.set_rasterized(False)
    cb.set_label("扩散系数常用对数 log10(D / (m²/s))")
    ax.set(xlabel="药材水分浓度（kg/kg）", ylabel="药材温度（℃）")
    export(fig, "raw_q4_diffusivity_relation")

    # process 1/4：参考网格映射到实际物理空间后的间距演化。
    fig, ax = plt.subplots(layout="constrained")
    for target_h, color, ls in zip((0.0, 24.0, event_h), COLORS[:3], (":", "--", "-")):
        idx = int(np.argmin(np.abs(hours - target_h)))
        radius_m = float(z["radius_time_m"][idx]) if target_h < event_h else float(z["event_radius_m"])
        physical = z["xi"] * radius_m
        ax.semilogy(0.5 * (physical[:-1] + physical[1:]) * 100.0,
                    np.diff(physical) * 1e6, color=color, ls=ls, lw=1.1,
                    label=(f"{target_h:g} h" if target_h < event_h else f"结束 {event_h:.2f} h"))
    ax.set(xlabel="到药材中心的距离（cm）", ylabel="相邻节点物理间距（μm）")
    ax.legend(ncol=3)
    export(fig, "process_q4_mesh_evolution")

    # process 2/4：空间与时间收敛。
    s200 = formal_summary["event"]["time_s"]
    s400 = json.loads((RESULTS / "shrink400_dt05_summary.json").read_text(encoding="utf-8"))["event"]["time_s"]
    s400fine = json.loads((RESULTS / "shrink400_dt025_summary.json").read_text(encoding="utf-8"))["event"]["time_s"]
    fig, axes = plt.subplots(1, 2, layout="constrained")
    axes[0].plot([200, 400], [s200 / 3600.0, s400 / 3600.0], "o-",
                 color=COLORS[0], lw=1.1)
    axes[0].set(xlabel="径向区间数 N", ylabel="临界烘干时间（h）")
    axes[0].ticklabel_format(axis="y", style="plain", useOffset=False)
    axes[1].plot([0.5, 0.25], [abs(s400 - s400fine) * 1000.0, 0.0], "s--",
                 color=COLORS[1], lw=1.1)
    axes[1].set(xlabel="时间步长（s）", ylabel="相对细步长的时刻差（ms）")
    axes[1].set_ylim(bottom=-0.02)
    export(fig, "process_q4_convergence", (7.2, 3.2))

    # process 3/4：4 h边界切换前后内部状态连续。
    mask = (hours >= 3.7) & (hours <= 4.3)
    fig, axes = plt.subplots(2, 1, sharex=True, layout="constrained")
    axes[0].plot(hours[mask], z["temperature_C"][mask, 0], color=COLORS[1], lw=1.1, label="中心")
    axes[0].plot(hours[mask], z["temperature_C"][mask, -1], color=COLORS[0], ls="--", lw=1.1, label="表面")
    axes[0].axvline(4.0, color="#666666", ls=":", lw=0.8)
    axes[0].set(ylabel="温度（℃）")
    axes[0].legend(ncol=2)
    axes[1].plot(hours[mask], z["moisture_kg_kg"][mask, 0], color=COLORS[1], lw=1.1)
    axes[1].plot(hours[mask], z["moisture_kg_kg"][mask, -1], color=COLORS[0], ls="--", lw=1.1)
    axes[1].axvline(4.0, color="#666666", ls=":", lw=0.8, label="边界函数切换")
    axes[1].set(xlabel="时间（h）", ylabel="水分浓度（kg/kg）")
    axes[1].legend(loc="upper right")
    export(fig, "process_q4_four_hour_continuity", (6.8, 5.0))

    # process 4/4：平均判据和全域判据的差别。
    fig, ax = plt.subplots(layout="constrained")
    ax.plot(hours, cmax, color=COLORS[0], lw=1.2, label="全域最大含水率")
    ax.plot(hours, cavg, color=COLORS[1], ls="--", lw=1.2, label="截面平均含水率")
    ax.axhline(0.15, color="#333333", ls=":", lw=0.9, label="达标阈值")
    mean_h = formal_summary["mean_threshold_time_s"] / 3600.0
    ax.axvspan(mean_h, event_h, color=COLORS[1], alpha=0.12,
               label=f"误停风险区间 {event_h-mean_h:.2f} h")
    ax.set(xlabel="时间（h）", ylabel="水分浓度（kg/kg）", xlim=(0, 53), ylim=(0, 2.65))
    ax.legend(loc="upper right")
    export(fig, "process_q4_average_vs_max")

    # result 1/4：临界事件局部图。
    fig, ax = plt.subplots(layout="constrained")
    mask = hours >= 30.0
    ax.plot(hours[mask], cmax[mask], color=COLORS[0], lw=1.4, label="全域最大含水率")
    ax.axhline(0.15, color=COLORS[1], ls="--", lw=1.0, label="阈值 0.15 kg/kg")
    ax.axvline(event_h, color=COLORS[2], ls=":", lw=1.0,
               label=f"临界时刻 {event_h:.4f} h")
    ax.scatter([event_h], [0.15], color=COLORS[2], s=24, zorder=4)
    ax.set(xlabel="时间（h）", ylabel="最大含水率（kg/kg）",
           xlim=(30, 52), ylim=(0.145, 0.23))
    ax.legend(loc="upper right")
    export(fig, "result_q4_Cmax_time")

    # result 2/4：实际物理坐标下的径向剖面，端点随R(t)移动。
    fig, axes = plt.subplots(2, 1, layout="constrained")
    for target_h, color, ls in zip((6.0, 18.0), COLORS[:2], ("-", "--")):
        idx = int(np.argmin(np.abs(hours - target_h)))
        axes[0].plot(z["xi"] * z["radius_time_m"][idx] * 100.0,
                     z["moisture_kg_kg"][idx], color=color, ls=ls, lw=1.15,
                     label=f"{target_h:g} h（R={z['radius_time_m'][idx]*100:.3f} cm）")
    axes[0].set(ylabel="水分浓度（kg/kg）", xlim=(0, 2.02))
    axes[0].legend()
    for target_h, color, ls in zip((30.0, 42.0), COLORS[2:4], ("-.", ":")):
        idx = int(np.argmin(np.abs(hours - target_h)))
        axes[1].plot(z["xi"] * z["radius_time_m"][idx] * 100.0,
                     z["moisture_kg_kg"][idx], color=color, ls=ls, lw=1.15,
                     label=f"{target_h:g} h（R={z['radius_time_m'][idx]*100:.3f} cm）")
    axes[1].plot(z["xi"] * float(z["event_radius_m"]) * 100.0,
                 z["event_moisture_kg_kg"], color="#222222", lw=1.35,
                 label=f"结束 {event_h:.2f} h（R={float(z['event_radius_m'])*100:.3f} cm）")
    axes[1].axhline(0.15, color="#777777", ls="--", lw=0.8)
    axes[1].set(xlabel="到药材中心的实际距离（cm）", ylabel="水分浓度（kg/kg）",
                xlim=(0, 2.02), ylim=(0.045, 0.24))
    axes[1].legend(ncol=2)
    export(fig, "result_q4_radial_profiles", (6.8, 5.6))

    # result 3/4：相同附录4物性下，仅改变半径的受控比较。
    fig, ax = plt.subplots(layout="constrained")
    ax.plot(hours, cmax, color=COLORS[0], lw=1.25,
            label=f"收缩半径：{event_h:.2f} h")
    ax.plot(fixed_hours, fixed_cmax, color=COLORS[1], ls="--", lw=1.25,
            label=f"固定2 cm：{fixed_event_h:.2f} h")
    ax.axhline(0.15, color="#333333", ls=":", lw=0.9)
    ax.axvspan(event_h, fixed_event_h, color=COLORS[2], alpha=0.10,
               label=f"缩短 {summary['shrink_time_reduction_h']:.2f} h（{summary['shrink_relative_reduction_percent']:.1f}%）")
    ax.set(xlabel="时间（h）", ylabel="最大含水率（kg/kg）", xlim=(0, 135), ylim=(0, 2.65))
    ax.legend(loc="upper right")
    export(fig, "result_q4_shrink_vs_fixed")

    # result 4/4：移动物理域中的时空水分场。
    take_t = slice(None, None, 5)
    take_x = slice(None, None, 2)
    x = np.tile(hours[take_t], (z["xi"][take_x].size, 1))
    y = z["xi"][take_x, None] * z["radius_time_m"][take_t][None, :] * 100.0
    field = z["moisture_kg_kg"][take_t, take_x].T
    fig, ax = plt.subplots(layout="constrained")
    levels = np.linspace(0.05, 2.55, 65)
    mesh = ax.contourf(x, y, field, levels=levels, cmap=MOISTURE_CMAP)
    cb = fig.colorbar(mesh, ax=ax)
    if cb.solids is not None:
        cb.solids.set_rasterized(False)
    cb.set_label("水分浓度（kg/kg）")
    ax.contour(x, y, field, levels=[0.15], colors="white", linewidths=0.9,
               linestyles="--")
    boundary_line, = ax.plot(hours, z["radius_time_m"] * 100.0,
                             color="#222222", lw=1.0)
    ax.set(xlabel="时间（h）", ylabel="到药材中心的实际距离（cm）",
           xlim=(0, event_h), ylim=(0, 2.02))
    threshold_proxy = Line2D([0], [0], color="#777777", lw=0.9, ls="--")
    ax.legend([boundary_line, threshold_proxy],
              ["移动表面 R(t)", "0.15 kg/kg 等值线"], loc="upper right")
    export(fig, "result_q4_moisture_field")

    print(f"已生成12张问题四中文候选图：{FIGURES}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
