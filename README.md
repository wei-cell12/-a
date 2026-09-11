# 2026 数学建模 A 题计算仓库

问题一和问题二已拆分为两个可独立复现的目录：

| 目录 | 内容 | 主要入口 |
| --- | --- | --- |
| [`question1/`](question1/) | 预热平衡阶段：常物性温度场与水分场、CN/后向欧拉对照 | [`question1/README.md`](question1/README.md) |
| [`question2/`](question2/) | 整个烘干过程：变物性热湿耦合模型、收敛验证 | [`question2/README.md`](question2/README.md) |

## 快速入口

- 问题一论文正文：[question1/问题一_论文正文.md](question1/问题一_论文正文.md)
- 问题一提交结果：[question1/results/result1.xlsx](question1/results/result1.xlsx)
- 问题二论文正文：[question2/Q2_paper_section.md](question2/Q2_paper_section.md)
- 问题二提交结果：[question2/results/result2.xlsx](question2/results/result2.xlsx)

## 运行方式

安装公共依赖：

```powershell
python -m pip install -r requirements.txt
```

进入相应目录后独立运行：

```powershell
cd question1
python code/q1_solver.py

cd ../question2
python code/q2_solver.py
```

两个目录内均保留各自所需的输入模板、完整结果、论文图片和依赖清单。
