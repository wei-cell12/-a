# 问题三：全域含水率达标时间

> **阶段性结果（2026-09-12）**：本目录保存第三问当前计算结果，已通过现阶段数值检验，但后续仍可能根据建模讨论、MATLAB 重绘和论文定稿继续修改；请勿视为最终提交版本。

正式结果：临界烘干时间为 **57.4740 h**；按整秒执行应取 **206907 s**。

## 目录

- `code/q3_solver.py`：有限体积 + Crank--Nicolson + Picard 主求解器。
- `code/export_q3_results.py`：生成 CSV、XLSX、MAT 和复现清单。
- `code/make_q3_figures.py`：生成 9 张中文论文候选图。
- `data/`：附件 1 与 result3 模板的只读副本。
- `results/`：正式结果、收敛计算与诊断。
- `figures/`：400 DPI PNG、SVG、PDF。
- `问题三_结果说明.md`：给论文手的完整结果解释。
- `问题三_论文正文.md`：按前两问风格组织、可直接并入总论文继续修改的第三问正文草稿。

## 复现

```powershell
python code/run_all_q3.py
```

该入口会依次执行三组原始 `p=2` 网格、三组优化 `p=1.75` 网格，随后生成表格和图件、刷新哈希清单并运行自动验证；因此可在没有历史结果文件的干净目录中完整复现。

本机 Python 3.9 若在中文路径下出现 `A题` 被显示为 `A��`，可先映射纯英文盘符再运行：

```powershell
subst Q: "当前仓库所在的 A题 目录"
Set-Location Q:\github_release_q1\question3
```

依赖见 `requirements.txt`。XLSX 使用题目模板的 OOXML 结构直接写入，不依赖 openpyxl。
