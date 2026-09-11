#!/usr/bin/env python3
"""从正式 CN CSV 生成问题一论文三联图。"""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

import q1_cn_comparison as solver
from utils.plot_style import COLOR_SEQUENCE, PALETTE, apply_publication_style


ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
SELECTED_TIMES = [100, 300, 600, 900, 1200, 1500, 1800]
SELECTED_RADII = [0, 0.5, 1.0, 1.5, 2.0]


def load_matrix(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    raw = np.genfromtxt(path, delimiter=",", skip_header=1)
    header = path.read_text(encoding="utf-8-sig").splitlines()[0].split(",")[1:]
    radii = np.array([float(item.removeprefix("r_").removesuffix("_cm")) for item in header])
    return raw[:, 0], radii, raw[:, 1:]


def export(fig: plt.Figure, name: str) -> None:
    fig.set_size_inches(7.2, 7.2)
    fig.savefig(FIGURES / f"{name}.svg")
    fig.savefig(FIGURES / f"{name}.png", dpi=300)
    plt.close(fig)


def make_triptych(csv_name: str, variable: str, unit: str, cmap: str, output_name: str) -> None:
    time_s, radius_cm, field = load_matrix(RESULTS / csv_name)
    fig = plt.figure(layout="constrained")
    grid = fig.add_gridspec(2, 2)
    history = fig.add_subplot(grid[0, 0])
    profiles = fig.add_subplot(grid[0, 1])
    heatmap = fig.add_subplot(grid[1, :])

    for index, radius in enumerate(SELECTED_RADII):
        column = int(np.argmin(np.abs(radius_cm - radius)))
        history.plot(time_s, field[:, column], color=COLOR_SEQUENCE[index],
                     linestyle=("-", "--", "-.", ":", (0, (5, 2)))[index],
                     label=f"r={radius:g} cm")
    history.set(xlabel="时间（s）", ylabel=f"{variable}（{unit}）", title="（a）典型位置时间历程")
    history.legend(fontsize=6)

    for index, target_time in enumerate(SELECTED_TIMES):
        row = int(np.argmin(np.abs(time_s - target_time)))
        profiles.plot(radius_cm, field[row], color=COLOR_SEQUENCE[index], label=f"{target_time} s")
    profiles.set(xlabel="到药材中心的距离（cm）", ylabel=f"{variable}（{unit}）", title="（b）指定时刻径向分布")
    profiles.legend(ncol=2, fontsize=6)

    image = heatmap.imshow(field, origin="lower", aspect="auto", cmap=cmap,
                           extent=[radius_cm[0], radius_cm[-1], time_s[0], time_s[-1]])
    colorbar = fig.colorbar(image, ax=heatmap)
    colorbar.set_label(f"{variable}（{unit}）")
    heatmap.set(xlabel="到药材中心的距离（cm）", ylabel="时间（s）", title="（c）径向时空分布")
    export(fig, output_name)


def make_interpolation_figures() -> None:
    raw = np.genfromtxt(RESULTS / "environment_q1.csv", delimiter=",", skip_header=1)
    time_s, temperature, moisture = raw[:, 0], raw[:, 1], raw[:, 2]
    styles = ("-", "--", "-.", ":")
    labels = {"linear": "线性", "pchip": "PCHIP", "cubic": "自然三次样条", "akima": "Akima"}
    for values, ylabel, output_name in (
        (temperature, "烘房温度（℃）", "raw_q1_temperature_interpolation"),
        (moisture, "烘房水分浓度（kg/kg）", "raw_q1_moisture_interpolation"),
    ):
        fig, ax = plt.subplots(layout="constrained")
        dense = np.linspace(time_s[0], time_s[-1], 1201)
        for index, method in enumerate(("linear", "pchip", "cubic", "akima")):
            ax.plot(dense, solver.make_interpolator(time_s, values, method)(dense),
                    color=COLOR_SEQUENCE[index], linestyle=styles[index], label=labels[method])
        ax.scatter(time_s, values, s=12, color="black", zorder=5, label="附件实测点")
        ax.set(xlabel="时间（s）", ylabel=ylabel, title="边界数据插值方法比较")
        ax.legend(ncol=3)
        export_sized(fig, output_name, (6.3, 3.9))


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
    export_sized(fig, "paper_q1_cn_temperature_3d", (7.2, 5.4))


def export_sized(fig: plt.Figure, name: str, size: tuple[float, float]) -> None:
    fig.set_size_inches(*size)
    fig.savefig(FIGURES / f"{name}.svg")
    fig.savefig(FIGURES / f"{name}.png", dpi=300)
    plt.close(fig)


def main() -> None:
    apply_publication_style()
    FIGURES.mkdir(exist_ok=True)
    make_triptych("问题1_CN_完整温度.csv", "温度", "℃", "coolwarm", "paper_q1_cn_temperature")
    make_triptych("问题1_CN_完整水分浓度.csv", "干基含水率", "kg/kg", "viridis", "paper_q1_cn_moisture")
    make_interpolation_figures()
    make_method_difference_figure()
    make_convergence_figure()
    make_temperature_3d()


if __name__ == "__main__":
    main()
