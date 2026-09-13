#!/usr/bin/env python3
"""从正式 CN CSV 生成问题一论文三联图"""

from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.hashsalt"] = "cumcm-2026-q1"
import matplotlib.pyplot as plt
import numpy as np

import q1_cn_comparison as solver
from utils.plot_style import COLOR_SEQUENCE, PALETTE, apply_publication_style


ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
SELECTED_TIMES = [100, 300, 600, 900, 1200, 1500, 1800]
SELECTED_RADII = [0, 0.5, 1.0, 1.5, 2.0]


# 读取数值结果矩阵并分离时间、半径和场变量
def load_matrix(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    raw = np.genfromtxt(path, delimiter=",", skip_header=1)
    header = path.read_text(encoding="utf-8-sig").splitlines()[0].split(",")[1:]
    radii = np.array([float(item.removeprefix("r_").removesuffix("_cm")) for item in header])
    return raw[:, 0], radii, raw[:, 1:]


# 按论文尺寸和分辨率导出当前图形
def export(fig: plt.Figure, name: str, size: tuple[float, float] = (7.2, 7.2)) -> None:
    fig.set_size_inches(*size)
    fig.savefig(FIGURES / f"{name}.svg", metadata={"Date": None})
    fig.savefig(FIGURES / f"{name}.png", dpi=300)
    plt.close(fig)


# 生成径向曲线、时空分布和三维曲面组成的三联图
def make_triptych(csv_name: str, variable: str, unit: str, cmap: str, output_name: str) -> None:
    time_s, radius_cm, field = load_matrix(RESULTS / csv_name)
    fig, (history, profiles, heatmap) = plt.subplots(
        1, 3,
        layout="constrained",
        gridspec_kw={"width_ratios": [1.0, 1.0, 1.15]},
    )

    for index, radius in enumerate(SELECTED_RADII):
        column = int(np.argmin(np.abs(radius_cm - radius)))
        history.plot(time_s, field[:, column], color=COLOR_SEQUENCE[index],
                     linestyle=("-", "--", "-.", ":", (0, (5, 2)))[index],
                     label=f"r={radius:g} cm")
    history.set(xlabel="时间（s）", ylabel=f"{variable}（{unit}）", title="（a）典型位置时间历程")
    history.legend(ncol=2, fontsize=6, frameon=False)

    for index, target_time in enumerate(SELECTED_TIMES):
        row = int(np.argmin(np.abs(time_s - target_time)))
        profiles.plot(radius_cm, field[row], color=COLOR_SEQUENCE[index], label=f"{target_time} s")
    profiles.set(xlabel="到药材中心的距离（cm）", ylabel=f"{variable}（{unit}）", title="（b）指定时刻径向分布")
    profiles.legend(ncol=2, fontsize=6, frameon=False)

    image = heatmap.imshow(field, origin="lower", aspect="auto", cmap=cmap,
                           extent=[radius_cm[0], radius_cm[-1], time_s[0], time_s[-1]])
    colorbar = fig.colorbar(image, ax=heatmap)
    colorbar.set_label(f"{variable}（{unit}）")
    heatmap.set(xlabel="到药材中心的距离（cm）", ylabel="时间（s）", title="（c）径向时空分布")
    export(fig, output_name, size=(12.0, 3.8))


# 生成环境数据插值效果及局部差异图
def make_interpolation_figures() -> None:
    raw = np.genfromtxt(RESULTS / "environment_q1.csv", delimiter=",", skip_header=1)
    time_s, temperature, moisture = raw[:, 0], raw[:, 1], raw[:, 2]
    styles = ("-", "--", "-.", ":")
    labels = {"linear": "线性", "pchip": "PCHIP", "cubic": "自然三次样条", "akima": "Akima"}
    for values, ylabel, difference_label, scale, output_name in (
        (temperature, "烘房温度（℃）", "相对线性插值差值（℃）", 1.0,
         "raw_q1_temperature_interpolation_detailed"),
        (moisture, "烘房水分浓度（kg/kg）", r"相对线性插值差值（$\times 10^{-5}$ kg/kg）", 1.0e5,
         "raw_q1_moisture_interpolation_detailed"),
    ):
        # 0.1 s 的致密采样只用于平滑显示，不改变求解器的边界输入或数值结果
        dense = np.arange(time_s[0], time_s[-1] + 0.05, 0.1)
        predictions = {
            method: np.asarray(solver.make_interpolator(time_s, values, method)(dense), dtype=float)
            for method in ("linear", "pchip", "cubic", "akima")
        }
        fig, (ax, difference_ax) = plt.subplots(
            2, 1, sharex=True, layout="constrained", gridspec_kw={"height_ratios": [2.15, 1.0]}
        )
        for index, method in enumerate(("linear", "pchip", "cubic", "akima")):
            ax.plot(dense, predictions[method], color=COLOR_SEQUENCE[index],
                    linestyle=styles[index], linewidth=1.35, label=labels[method])
        ax.scatter(time_s, values, s=12, color="black", zorder=5, label="附件实测点")
        ax.set(ylabel=ylabel, title="边界数据插值方法比较：整体曲线与局部差值")
        ax.legend(ncol=3, fontsize=7)

        baseline = predictions["linear"]
        for index, method in enumerate(("pchip", "cubic", "akima"), start=1):
            difference_ax.plot(
                dense, (predictions[method] - baseline) * scale,
                color=COLOR_SEQUENCE[index], linestyle=styles[index], linewidth=1.15,
                label=labels[method],
            )
        difference_ax.axhline(0.0, color="0.25", linewidth=0.7)
        difference_ax.set(xlabel="时间（s）", ylabel=difference_label,
                          title="局部放大：各平滑插值减去线性插值")
        difference_ax.legend(ncol=3, fontsize=7)
        difference_ax.margins(x=0)
        export_sized(fig, output_name, (7.2, 5.8))


# 绘制 Crank--Nicolson 与后向欧拉结果差异图
def make_method_difference_figure() -> None:
    be_t, _, be_temperature = load_matrix(RESULTS / "问题1_BE_指定时刻温度.csv")
    _, _, cn_temperature = load_matrix(RESULTS / "问题1_CN_指定时刻温度.csv")
    be_c, _, be_moisture = load_matrix(RESULTS / "问题1_BE_指定时刻水分浓度.csv")
    _, _, cn_moisture = load_matrix(RESULTS / "问题1_CN_指定时刻水分浓度.csv")
    fig, axes = plt.subplots(1, 2, layout="constrained")
    axes[0].plot(be_t, np.max(np.abs(be_temperature - cn_temperature), axis=1), marker="o", color=PALETTE["primary"])
    axes[0].set(xlabel="时间（s）", ylabel="最大绝对差（℃）", title="后向欧拉与 CN：温度")
    axes[1].plot(be_c, np.max(np.abs(be_moisture - cn_moisture), axis=1), marker="o", color=PALETTE["primary"])
    axes[1].set(xlabel="时间（s）", ylabel="最大绝对差（kg/kg）", title="后向欧拉与 CN：水分")
    export_sized(fig, "process_q1_be_cn_difference", (7.2, 3.6))


# 绘制网格与时间步长收敛性检验图
def make_convergence_figure() -> None:
    raw = np.genfromtxt(RESULTS / "问题1_CN收敛性.csv", delimiter=",", skip_header=1,
                        dtype=None, encoding="utf-8-sig")
    rows = [[row[0], int(row[1]), float(row[2]), float(row[3]), float(row[4])] for row in raw]
    grid_rows = [row for row in rows if row[0] == "grid"][:-1]
    time_rows = [row for row in rows if row[0] == "time"][:-1]
    fig, axes = plt.subplots(2, 2, layout="constrained")
    panels = (
        (axes[0, 0], grid_rows, 1, 3, "径向区间数 N", "温度最大差（℃）", "温度：空间网格"),
        (axes[1, 0], grid_rows, 1, 4, "径向区间数 N", "水分最大差（kg/kg）", "水分：空间网格"),
        (axes[0, 1], time_rows, 2, 3, "时间步长（s）", "温度最大差（℃）", "温度：时间步长"),
        (axes[1, 1], time_rows, 2, 4, "时间步长（s）", "水分最大差（kg/kg）", "水分：时间步长"),
    )
    for ax, values, xcol, ycol, xlabel, ylabel, title in panels:
        ax.loglog([row[xcol] for row in values], [row[ycol] for row in values], marker="o")
        ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
        if xcol == 2:
            ax.invert_xaxis()
    export_sized(fig, "process_q1_cn_convergence", (7.2, 5.4))


# 绘制表面与中心状态差异随时间的变化
def make_surface_center_gap_figure() -> None:
    """展示预热过程中表面与中心梯度的形成，作为求解过程证据"""
    time_s, _, temperature = load_matrix(RESULTS / "问题1_CN_完整温度.csv")
    _, _, moisture = load_matrix(RESULTS / "问题1_CN_完整水分浓度.csv")
    fig, axes = plt.subplots(1, 2, layout="constrained")
    axes[0].plot(time_s, temperature[:, -1] - temperature[:, 0], color=PALETTE["primary"])
    axes[0].set(xlabel="时间（s）", ylabel="表面温度−中心温度（℃）",
                title="径向温差的形成")
    axes[1].plot(time_s, moisture[:, 0] - moisture[:, -1], color=PALETTE["secondary"])
    axes[1].set(xlabel="时间（s）", ylabel="中心含水率−表面含水率（kg/kg）",
                title="径向含水率差的形成")
    export_sized(fig, "process_q1_surface_center_gap", (7.2, 3.6))


# 绘制温度场随时间和半径变化的三维曲面
def make_temperature_3d() -> None:
    time_s, radius_cm, field = load_matrix(RESULTS / "问题1_CN_完整温度.csv")
    stride = 10
    radius_grid, time_grid = np.meshgrid(radius_cm, time_s[::stride])
    fig = plt.figure(layout="constrained")
    ax = fig.add_subplot(111, projection="3d")
    surface = ax.plot_surface(radius_grid, time_grid, field[::stride], cmap="coolwarm",
                              linewidth=0, antialiased=True, rcount=120, ccount=21)
    ax.set(xlabel="到药材中心的距离（cm）", ylabel="时间（s）", zlabel="温度（℃）",
           title="药材温度场三维时空分布")
    ax.view_init(elev=28, azim=-128)
    colorbar = fig.colorbar(surface, ax=ax, shrink=0.68, pad=0.08)
    colorbar.set_label("温度（℃）")
    export_sized(fig, "result_q1_cn_temperature_3d", (7.2, 5.4))


# 按指定版面尺寸导出论文图形
def export_sized(fig: plt.Figure, name: str, size: tuple[float, float]) -> None:
    fig.set_size_inches(*size)
    fig.savefig(FIGURES / f"{name}.svg", metadata={"Date": None})
    fig.savefig(FIGURES / f"{name}.png", dpi=300)
    plt.close(fig)


# 组织当前脚本的完整执行流程并返回运行状态
def main() -> None:
    apply_publication_style()
    FIGURES.mkdir(exist_ok=True)
    make_triptych("问题1_CN_完整温度.csv", "温度", "℃", "coolwarm", "result_q1_cn_temperature")
    # plasma 在不扭曲数值映射的前提下增强低含水率紫色端与高含水率黄色端的视觉反差
    make_triptych("问题1_CN_完整水分浓度.csv", "干基含水率", "kg/kg", "plasma", "result_q1_cn_moisture")
    make_interpolation_figures()
    make_method_difference_figure()
    make_convergence_figure()
    make_surface_center_gap_figure()
    make_temperature_3d()


if __name__ == "__main__":
    main()
