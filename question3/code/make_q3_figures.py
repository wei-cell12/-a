#!/usr/bin/env python3
"""生成第三问 raw/process/result 三类中文论文候选图。"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from PIL import Image

from q3_solver import Config, diffusivity, radial_geometry
from utils.plot_style import COLORS, apply_style


ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
QA = FIGURES / "_qa"
FORMAL = "optimized400_dt05"
MOISTURE_CMAP = LinearSegmentedColormap.from_list(
    "blue_yellow_red", ["#1646A0", "#3AAED8", "#F2D45C", "#C7392F"], N=256)


def export(fig: plt.Figure, name: str, size: tuple[float, float] = (6.8, 4.2)) -> None:
    fig.set_size_inches(*size)
    fig.savefig(FIGURES / f"{name}.png", dpi=400)
    fig.savefig(FIGURES / f"{name}.svg")
    fig.savefig(FIGURES / f"{name}.pdf")
    # 灰度预览只用于可辨识度质检，不作为论文图。
    with Image.open(FIGURES / f"{name}.png") as source:
        source.convert("L").save(QA / f"{name}_grayscale.png")
    plt.close(fig)


def load_formal() -> tuple[dict, dict]:
    summary = json.loads((RESULTS / f"{FORMAL}_summary.json").read_text(encoding="utf-8"))
    arrays = dict(np.load(RESULTS / f"{FORMAL}_fields.npz"))
    return summary, arrays


def main() -> int:
    FIGURES.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    apply_style()
    summary, z = load_formal()
    env = np.genfromtxt(RESULTS / "q3_environment.csv", delimiter=",", skip_header=1,
                        encoding="utf-8-sig")
    hours = z["time_s"] / 3600.0
    radius_cm = z["radius_m"] * 100.0
    weights = z["volume_m2"] / (0.5 * 0.02**2)
    cmax = np.max(z["moisture_kg_kg"], axis=1)
    cavg = z["moisture_kg_kg"] @ weights
    event_h = float(z["event_time_s"]) / 3600.0

    # raw 1/3：实测温度时间序列。
    fig, ax = plt.subplots(layout="constrained")
    ax.plot(env[:, 0] / 3600.0, env[:, 1], color=COLORS[0], lw=1.2)
    ax.scatter(env[::10, 0] / 3600.0, env[::10, 1], s=8, color=COLORS[0])
    ax.axhline(50.0, color=COLORS[1], ls="--", lw=0.9, label="后期稳定近似 50 ℃")
    ax.set(xlabel="时间（h）", ylabel="烘房温度（℃）", xlim=(0, 4))
    ax.legend(loc="lower right")
    export(fig, "raw_q3_boundary_temperature")

    # raw 2/3：实测环境水分时间序列。
    fig, ax = plt.subplots(layout="constrained")
    ax.plot(env[:, 0] / 3600.0, env[:, 2], color=COLORS[2], lw=1.2)
    ax.scatter(env[::10, 0] / 3600.0, env[::10, 2], s=8, color=COLORS[2])
    ax.axhline(0.05, color=COLORS[1], ls="--", lw=0.9, label="后期稳定近似 0.05 kg/kg")
    ax.set(xlabel="时间（h）", ylabel="烘房水分浓度（kg/kg）", xlim=(0, 4))
    ax.legend(loc="upper right")
    export(fig, "raw_q3_boundary_moisture")

    # raw 3/3：经验扩散系数的二维关系，说明末期扩散变慢。
    c_axis = np.linspace(0.05, 2.55, 220)
    t_axis = np.linspace(28.0, 50.0, 150)
    cc, tt = np.meshgrid(c_axis, t_axis)
    dd = diffusivity(tt, cc)
    fig, ax = plt.subplots(layout="constrained")
    log_dd = np.log10(dd)
    diffusion_levels = np.linspace(log_dd.min(), log_dd.max(), 65)
    mesh = ax.contourf(c_axis, t_axis, log_dd, levels=diffusion_levels, cmap="cividis")
    cb = fig.colorbar(mesh, ax=ax)
    if cb.solids is not None:
        cb.solids.set_rasterized(False)
    cb.set_label("扩散系数常用对数 log10(D / (m²/s))")
    ax.set(xlabel="药材水分浓度（kg/kg）", ylabel="药材温度（℃）")
    export(fig, "raw_q3_diffusivity_relation")

    # process 1/3：网格间距比较，直接展示 p=2 过密问题。
    fig, ax = plt.subplots(layout="constrained")
    for power, color, ls in zip((1.0, 1.75, 2.0), COLORS[:3], (":", "-", "--")):
        r, _, _ = radial_geometry(Config(intervals=400, mesh_power=power))
        ax.semilogy(0.5 * (r[:-1] + r[1:]) * 100.0, np.diff(r) * 1e6,
                    color=color, ls=ls, lw=1.1, label=f"p={power:g}")
    ax.set(xlabel="径向位置（cm）", ylabel="相邻节点间距（μm）")
    ax.legend(title="加密指数", ncol=3)
    export(fig, "process_q3_mesh_spacing")

    # process 2/3：空间和时间收敛性，两个面板避免双Y轴。
    with (RESULTS / "q3_convergence.csv").open(encoding="utf-8-sig") as stream:
        rows = list(csv.reader(stream))[1:]
    opt = [row for row in rows if row[0] == "优化加密"]
    fig, axes = plt.subplots(1, 2, layout="constrained")
    spatial = sorted([row for row in opt if float(row[3]) == 0.5], key=lambda x: int(float(x[1])))
    axes[0].plot([int(float(x[1])) for x in spatial], [float(x[4]) for x in spatial],
                 "o-", color=COLORS[0], lw=1.1)
    axes[0].set(xlabel="径向区间数 N", ylabel="临界烘干时间（h）")
    temporal = sorted([row for row in opt if int(float(row[1])) == 400],
                      key=lambda x: float(x[3]), reverse=True)
    axes[1].plot([float(x[3]) for x in temporal], [abs(float(x[5])) * 1000 for x in temporal],
                 "s--", color=COLORS[1], lw=1.1)
    axes[1].set(xlabel="时间步长（s）", ylabel="相对正式结果的时刻差（ms）")
    axes[0].ticklabel_format(axis="y", style="plain", useOffset=False)
    axes[1].set_ylim(bottom=-0.02)
    export(fig, "process_q3_convergence", (7.2, 3.2))

    # process 3/3：平均判据与全域判据对照。
    fig, ax = plt.subplots(layout="constrained")
    ax.plot(hours, cmax, color=COLORS[0], lw=1.2, label="全域最大含水率")
    ax.plot(hours, cavg, color=COLORS[1], ls="--", lw=1.2, label="截面平均含水率")
    ax.axhline(0.15, color="#333333", ls=":", lw=1.0, label="达标阈值")
    mean_h = summary["mean_threshold_time_s"] / 3600.0
    ax.axvspan(mean_h, event_h, color=COLORS[1], alpha=0.12,
               label=f"提前停机风险区间（{event_h-mean_h:.2f} h）")
    ax.set(xlabel="时间（h）", ylabel="水分浓度（kg/kg）", xlim=(0, 60), ylim=(0, 2.65))
    ax.legend(loc="upper right")
    export(fig, "process_q3_average_vs_max")

    # result 1/3：直接回答第三问的最大含水率时间曲线，使用局部纵轴突出交点。
    fig, ax = plt.subplots(layout="constrained")
    mask = hours >= 24.0
    ax.plot(hours[mask], cmax[mask], color=COLORS[0], lw=1.4, label="最大含水率（与中心重合）")
    ax.axhline(0.15, color=COLORS[1], ls="--", lw=1.1, label="阈值 0.15 kg/kg")
    ax.axvline(event_h, color=COLORS[2], ls=":", lw=1.1,
               label=f"临界时刻 {event_h:.4f} h")
    ax.scatter([event_h], [0.15], color=COLORS[2], s=24, zorder=4)
    ax.set(xlabel="时间（h）", ylabel="最大含水率（kg/kg）", xlim=(24, 59), ylim=(0.145, 0.245))
    ax.legend(loc="upper right")
    export(fig, "result_q3_Cmax_time")

    # result 2/3：多时刻径向剖面。
    profile_hours = [6, 18, 30, 42, 54]
    fig, axes = plt.subplots(2, 1, sharex=True, layout="constrained")
    for index, target_h in enumerate(profile_hours[:2]):
        row = int(np.argmin(np.abs(hours - target_h)))
        axes[0].plot(radius_cm, z["moisture_kg_kg"][row], color=COLORS[index],
                ls=("-", "--", "-.", ":", "-")[index], lw=1.15, label=f"{target_h:g} h")
    axes[0].set(ylabel="水分浓度（kg/kg）", xlim=(0, 2))
    axes[0].legend(ncol=2)
    for index, target_h in enumerate(profile_hours[2:], 2):
        row = int(np.argmin(np.abs(hours - target_h)))
        axes[1].plot(radius_cm, z["moisture_kg_kg"][row], color=COLORS[index],
                     ls=("-", "--", "-.", ":", "-")[index], lw=1.15,
                     label=f"{target_h:g} h")
    axes[1].plot(radius_cm, z["event_moisture_kg_kg"], color="#222222", lw=1.4,
                 label=f"结束 {event_h:.4f} h")
    axes[1].axhline(0.15, color="#777777", ls="--", lw=0.8, label="阈值 0.15")
    axes[1].set(xlabel="到药材中心的距离（cm）", ylabel="水分浓度（kg/kg）",
                xlim=(0, 2), ylim=(0.045, 0.215))
    axes[1].legend(ncol=3)
    export(fig, "result_q3_radial_profiles", (6.8, 5.6))

    # result 3/3：时空热力图，蓝—黄—红高对比色带，避免黑色起点。
    fig, ax = plt.subplots(layout="constrained")
    # 仅为矢量绘图降采样；数值表格和阈值事件仍使用全部节点与60 s数据。
    plot_hours = hours[::10]
    plot_radius = radius_cm[::2]
    field = z["moisture_kg_kg"][::10, ::2].T
    moisture_levels = np.linspace(0.05, 2.55, 65)
    mesh = ax.contourf(plot_hours, plot_radius, field, levels=moisture_levels,
                       cmap=MOISTURE_CMAP)
    cb = fig.colorbar(mesh, ax=ax)
    if cb.solids is not None:
        cb.solids.set_rasterized(False)
    cb.set_label("水分浓度（kg/kg）")
    ax.contour(plot_hours, plot_radius, field, levels=[0.15], colors="white",
               linewidths=0.9, linestyles="--")
    ax.set(xlabel="时间（h）", ylabel="到药材中心的距离（cm）",
           xlim=(0, event_h), ylim=(0, 2))
    export(fig, "result_q3_moisture_field")

    print(f"已生成 9 张中文论文候选图，目录：{FIGURES}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
