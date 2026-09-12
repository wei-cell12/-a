# 第三问工作笔记

## 输入审计（2026-09-12）

- 题面要求：各处含水率低于 0.15 kg/kg；论文表 5 每 6 h、每 0.5 cm；`result3.xlsx` 每 60 s、每 0.1 cm。
- 物性：沿用附录 3，扩散系数中的温度换算为 K，指数为 `-0.45/C` 与 `-3850/(T+273.15)`。
- 边界：前 4 h 读取附件 1 并线性插值；4 h 后的 50 ℃、0.05 kg/kg 是稳定平台假设，需在结果说明中标明。
- 几何：第三问固定半径 2 cm，收缩属于第四问。
- 判据：程序检查 `max_i C_i < 0.15`，并独立记录径向单调性及最大值位置。
- 既有问题二代码只支持均匀网格和固定 3 h 终点；第三问单独实现非均匀几何、事件定位、60 s 导出及长期诊断。
- 发现外部临时代码曾用 30 s 长步长得到中间结果；本次正式收敛检验不直接采用该值。

## 环境

- Python 3.9.9；NumPy、SciPy、Pandas、Matplotlib 可用。
- 环境检查提示缺少 openpyxl；本项目沿用模板 ZIP/XML 写入方式，因此不把 openpyxl 设为运行依赖。
- Python 3.9.9 在原始中文绝对路径中把 `A题` 错解码为 `A��`；通过 `subst Q:` 将项目映射到纯 ASCII 路径后执行，计算与写出均正常。

## 最小真实算例

- 命令：`python code/q3_solver.py --smoke --quiet`（在 `Q:/github_release_q1/question3`）。
- 退出码：0。
- 配置：表面加密 N=20、dt=2 s、计算至 600 s。
- 最大代数残差：1.1703e-12；最大单步含水率积分闭合误差：2.5923e-16。
- Picard 平均/最大迭代次数：4/4；负含水率、NaN/Inf、径向单调性破坏均为 0。

## 问题二—问题三衔接检验

- 配置：表面加密 N=200、dt=0.5 s，连续从 0 推进到 3 h。
- 3 h 五点含水率：`(1.76618172, 1.71652724, 1.57019104, 1.33326347, 1.00813115)` kg/kg。
- 与交接值 `(1.7662, 1.7165, 1.5702, 1.3333, 1.0081)` 的最大绝对差为 3.65e-05 kg/kg，逐项四舍五入一致。
- 命令退出码：0；运行时间 27.32 s；全过程无负值、非有限值或径向单调性破坏。

## P1 与全量计算

- 独立 P1 两轮均 PASS，P0=0、P1=0；第二轮覆盖温度稳态锁定逻辑。
- 原交接网格 p=2：N=200、dt=0.5 s 为 206917.4673 s；N=400、dt=0.5 s 为 206904.5787 s；N=400、dt=0.25 s 为 206904.5785 s。
- p=2、N=400 的最外层节点间距仅 0.125 μm，CN 在最外层出现微小交替振荡，故不作为正式网格。
- 优化网格 p=1.75：N=200、dt=0.5 s 为 206924.8536 s；N=400、dt=0.5 s 为 206906.4326 s；N=400、dt=0.25 s 为 206906.4329 s。
- 正式采用 N=400、p=1.75、dt=0.5 s：空间加密差 18.4210 s（8.90e-5），时间步减半差 0.000305 s（1.47e-9），全过程单调性破坏为 0。
- 平均含水率在 36.0266 h 达标，比全域判据早 21.4474 h。

## 输出与可视化

- 题目结果 `result3.xlsx` 共 3450 行、22 列：表头 + 3448 个规则 60 s 时刻 + 事件末行；ZIP结构、首末行、列数均通过解析校验。
- 表5 `q3_result_table.xlsx` 共 11 行、6 列；收敛表共 7 行、8 列。
- 生成 raw/process/result 各 3 张中文候选图，均有 400 DPI PNG、SVG、PDF；关键三张与两张过程图已人工读图检查。
- 科研图机器检查：主 PNG/SVG 均 PASS；水分场 SVG 通过降采样保留矢量网格，无 base64 位图。

## 错误与处理

- `check_env.py` 初次遗漏 `--features`，退出码 2；补为 `--features data excel visualization` 后完成审计。
- `openpyxl` 未安装，沙箱内外两次 pip 安装均因网络/代理失败（退出码 1）；改用项目内 OOXML 模板写入，不影响交付。
- 单独工具调用中 `Q:` 映射不会跨进程保留，曾导致一次 `Set-Location` 失败（退出码 1）；后续每条命令内重新执行 `subst`。

## 输入资料、模型合同与验证状态

- 输入资料已审计：附件 1、第三问结果模板及两份交接文档均已核对；输入哈希见 `results/q3_repro_manifest.json`。
- 模型合同已冻结：固定半径圆柱、问题二变物性热湿模型、4 h 后稳定边界、全域最大含水率阈值事件。
- 运行与验证已完成：六组收敛计算、正式结果导出、图件审计及 `code/verify_q3.py` 自动检查均已通过。
# 最终交付记录

- 正式临界时刻：`206906.432617 s = 57.474009 h`；若按整数秒执行，取 `206907 s`。
- 正式网格：`N=400, p=1.75, dt=0.5 s`。空间细化差 `18.4210 s`（相对 `8.90e-5`），时间步细化差 `0.000305 s`。
- 全域最大含水率在中心取得；平均含水率判据会提前 `21.4474 h` 停机，因此不能代替全域判据。
- 自动验证 `code/verify_q3.py` 全部通过；9 张中文论文候选图通过 PNG/SVG 文件检查和严格图件审计，SVG 不含嵌入位图。
- `code/run_all_q3.py` 是干净目录下的唯一全流程复现入口，覆盖六组求解、导出、绘图、哈希刷新和验证。

## 论文手主张—证据映射（W1 输入）

| 论文主张 | 精确证据路径 | 章节与编号落点 | 正式图表及首次引用 |
|---|---|---|---|
| 第三问在固定半径 2 cm 下承接问题二变物性热湿模型 | 问题二权威正文 `../question2/Q2_paper_section.md` 第 3 节；问题二代码 `../question2/code/q2_solver.py`；第三问代码 `code/q3_solver.py::{volumetric_heat_capacity,conductivity,diffusivity}` | 第 2 节；式 (1)–(5) | 第 2 节不单独插图 |
| 4 h 后采用 50 ℃、0.05 kg/kg 稳定环境是延拓假设 | `code/q3_solver.py::environment` 与正式配置 `results/optimized400_dt05_summary.json::config::{post_temperature_c,post_moisture}` | 第 2.3 节；式 (8) | 第 2.3 节首次引用图 1 `raw_q3_boundary_temperature.png` 与图 2 `raw_q3_boundary_moisture.png` |
| 停止条件是全域最大含水率首次低于 0.15 kg/kg | `results/optimized400_dt05_summary.json::event`；`code/q3_solver.py::{simulate,locate_event}` | 第 2.2 节；式 (6)–(7) | 第 4.1 节首次引用图 4 `result_q3_Cmax_time.png` |
| 连续临界时刻为 206906.432617 s，即 57.474009 h；按整秒执行取 206907 s | `results/q3_diagnostics.json::{critical_time_s,critical_time_h,strict_integer_second_s,event_max_C}` | 第 4.1 节结论框 | 图 4；表 1 `q3_result_table.csv` |
| 全过程最湿点位于中心，中心湿芯控制总时长 | `results/optimized400_dt05_summary.json::event::argmax_node=0`、`diagnostics::{radial_monotonicity_violations,max_center_vs_global_gap}` | 第 4.2 节 | 图 5 `result_q3_radial_profiles.png`、图 6 `result_q3_moisture_field.png`，均在第 4.2 节首次引用 |
| 平均含水率不能代替全域判据 | `results/q3_diagnostics.json::{mean_threshold_time_h,premature_stop_gap_h}`；`code/q3_solver.py::area_average` | 第 4.3 节；式 (15) | 图 7 `process_q3_average_vs_max.png` 在第 4.3 节首次引用 |
| 正式离散对事件时刻已基本收敛 | `results/q3_convergence.csv`；`results/q3_final_summary.xlsx` | 第 5.1 节；表 2 | 图 8 `process_q3_convergence.png` 在第 5.1 节首次引用 |
| 在 N=400、dt=0.5 s 对照下，p=1.75 网格相较 p=2 消除可检出的表层交替振荡，事件时刻仅变 1.8539 s | `results/q3_convergence.csv` 对应两行；`results/optimized400_dt05_summary.json::diagnostics::minimum_dr_m` | 第 3.1 节；式 (9) | 图 3 `process_q3_mesh_spacing.png` 在第 3.1 节首次引用 |
| 非线性迭代全部收敛，未出现负值或非有限值，代数残差和累计质量闭合误差足够小 | `results/q3_diagnostics.json::formal_run`；`results/q3_verification_report.json::checks` | 第 5.2 节；表 3 | 不单独成图 |
| 正式解在 3 h 的五点与问题二四位小数表一致，最大绝对差 3.7116e-5 kg/kg | `results/optimized400_dt05_fields.npz` 中 `time_s=10800` 的五点线性提取值 `(1.76617802,1.71653107,1.57019723,1.33326288,1.00812900)`；参考 `../question2/Q2_paper_section.md` 第 5.2 节表；对四位小数参考值最大差 `3.7116e-5` | 第 5.3 节 | 表 4，首次且仅在第 5.3 节引用 |

拟用正文结构：问题分析→模型承接、环境与判据→非均匀有限体积/CN/Picard→结果→判据辨析→收敛与可靠性→结论与边界。正文只引用现有真实结果，不把数值检验写成实验验证。
