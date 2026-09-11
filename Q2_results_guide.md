# 问题二结果文件说明

- `Q2_paper_section.md`：可直接交给论文手使用的第二问正文。
- `results/result2.xlsx`：题目要求的完整结果，两个工作表均为 10801 行×22 列。
- `results/q2_table3_temperature.csv`：论文表 3，时间单位为 h。
- `results/q2_table4_moisture.csv`：论文表 4，时间单位为 h。
- `results/q2_full_temperature.csv`、`q2_full_moisture.csv`：与 Excel 对应的纯文本备份。
- `results/q2_convergence.csv`：网格和时间步收敛性。
- `results/q2_nonlinear_solver_comparison.csv`：Picard 与 Newton 对照。
- `results/q2_validation_summary.json`：配置、关键结果和自动校验指标。
- `figures/*q2*.png`：可直接插入 Word 的中文图片。
- `figures/*q2*.svg`：可继续编辑的矢量版本；场图的颜色层已栅格化以控制文件大小。

复现命令：

```powershell
python code/q2_solver.py
```

仅重新生成图片：

```powershell
python code/q2_solver.py --figures-only
```

快速检查：

```powershell
python code/q2_solver.py --smoke
```
