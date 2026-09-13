# 论文附录核心代码

本目录只保留第 1—4 问中适合放入论文附录的核心算法片段，不包含数据、图片、结果文件、测试脚本或完整工程代码。

## 文件内容

| 文件 | 附录重点 | 对应完整源文件 |
|---|---|---|
| 问题一_核心代码.py | 圆柱径向控制体、空间算子、Crank–Nicolson 三对角求解、扩散系数 | github_release_q1/question1/code/q1_cn_comparison.py |
| 问题二_核心代码.py | 变物性系数、调和平均、CN 离散、Picard 热湿耦合迭代 | github_release_q1/question2/code/q2_solver.py |
| 问题三_核心代码.py | 表层加密非均匀网格、长时程推进、全域达标时刻二分定位 | github_release_q1/question3/code/q3_solver.py |
| 问题四_核心代码.py | 收缩半径插值、参考坐标算子、移动边界 CN 推进、临界时刻定位 | github_release_q1/question4/code/q4_solver.py |

每个 Python 文件开头均明确标注“论文附录节选，不作为独立运行入口”。代码由实际求解器原样抽取，只删去了附件读取、批量导表、绘图和命令行部分。

## 放入 LaTeX

1. 把 附录_导言区配置.tex 中的内容复制到论文导言区。
2. 在论文正文末尾需要放附录的位置输入 附录_核心代码.tex。
3. 保持本目录与主 tex 文件的相对位置不变；若移动目录，同步修改四个 lstinputlisting 路径。

如果正文已经加载 listings 或 xcolor，只保留一次 usepackage，避免重复。

