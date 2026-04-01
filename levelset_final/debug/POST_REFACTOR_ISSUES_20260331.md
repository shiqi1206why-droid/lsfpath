# 改完代码后出现的问题汇总（2026-03-31）

## 1. 结论摘要
- 问题 1（高优先级）：在严格守护阈值下（0 容差），`exact` 梯度链会导致迭代冻结（`accepted_steps=0`）。
- 问题 2（高优先级）：放宽守护阈值后，柔度会出现“先下降后上升”，优化后段依赖 `best-state rollback` 才返回最优结果。
- 问题 3（中优先级）：候选选择中 `theta_only` 在后段主导，HJ 候选长期未接管，几何更新贡献不足。
- 问题 4（中优先级）：历史序列语义存在“状态序列 vs 步序列”长度差，容易导致诊断误判。
- 问题 5（低优先级）：WSL/UNC 路径警告持续污染 MATLAB 批处理日志。

## 2. 运行证据

### 2.1 严格守护下迭代冻结（fast）
- 运行：`rerun_fast_20260331_161126`
- 证据文件：`/home/again/projects/recover/debug/refactor_artifacts/rerun_fast_20260331_161126/summary.txt`
- 关键结果：
  - `accepted_steps=0`
  - `rejected_steps=19`
  - `final_compliance=9.022712343230e-07`
  - `best_iter=1`
- 结论：优化没有实质推进，状态保持在初始附近。

### 2.2 放宽守护后可推进但出现后段回升（default）
- 运行：`rerun_default_20260331_170739`
- 证据文件：`/home/again/projects/recover/debug/refactor_artifacts/rerun_default_20260331_170739/summary.txt`
- 关键结果：
  - `accepted_steps=28`
  - `rejected_steps=0`
  - `best_iter=9`
  - `final_compliance=8.867999187144e-07`
- 进一步统计（来自 `rerun_default_results.mat`）：
  - `raw_final_compliance=9.194254437726e-07`
  - `final_to_best_gap_percent=3.679018`
  - `accepted_source`: `theta_only=17`, `reinit=11`, `hj=0`
- 结论：真实末态比历史最优差约 `3.68%`，最终输出依赖 rollback 回到最优。

### 2.3 放宽守护后 fast 同样恢复推进
- 运行：`rerun_fast_relaxed_20260331_165312`
- 证据文件：`/home/again/projects/recover/debug/refactor_artifacts/rerun_fast_relaxed_20260331_165312/summary.txt`
- 关键结果：
  - `accepted_steps=17`
  - `rejected_steps=2`
  - `final_compliance=8.921737585770e-07`
- 结论：冻结问题主要由“严格守护 + 新梯度分布”触发。

## 3. 根因定位（对应代码）

### 3.1 守护阈值与新梯度分布耦合
- 默认基础守护为 0 容差：
  - `/home/again/projects/recover/debug/config/get_fiber_optimization_params.m:42`
  - `/home/again/projects/recover/debug/config/get_fiber_optimization_params.m:44`
  - `/home/again/projects/recover/debug/config/get_fiber_optimization_params.m:47`
- 当前基础梯度链默认是 `exact`：
  - `/home/again/projects/recover/debug/config/get_fiber_optimization_params.m:95`
- 影响：候选步更容易触发 `next/current guard`，在 0 容差下表现为全拒绝冻结。

### 3.2 `theta_only` 接受规则允许“可接受恶化”
- 规则：`theta_only_compliance <= compliance * (1 + current_state_tol)`
  - `/home/again/projects/recover/debug/level_set_evolution/fiber_run_optimization_iterations.m:229`
- 现状：`default/fast` 被临时放宽到 `0.01`（1%）
  - `/home/again/projects/recover/debug/config/get_fiber_optimization_params.m:145`
  - `/home/again/projects/recover/debug/config/get_fiber_optimization_params.m:155`
- 影响：在局部最优之后，`theta_only` 仍可被接受，导致柔度后段上升。

### 3.3 候选选择对 HJ 的替换门槛较高
- 选择入口：`fiber_select_candidate_state`
  - `/home/again/projects/recover/debug/level_set_evolution/fiber_select_candidate_state.m:1`
- HJ 需显著优于已选候选才会替换：
  - `/home/again/projects/recover/debug/level_set_evolution/fiber_select_candidate_state.m:29`
- 本次 default 统计：`hj=0`，说明后段几何更新未主导。

### 3.4 末态质量由 rollback 兜底
- `raw_final_compliance` 与 `best_state` 比较后回滚：
  - `/home/again/projects/recover/debug/utilities/fiber_finalize_iteration_outputs.m:184`
  - `/home/again/projects/recover/debug/utilities/fiber_finalize_iteration_outputs.m:195`
  - `/home/again/projects/recover/debug/utilities/fiber_finalize_iteration_outputs.m:199`
  - 输出字段：`rollback_to_best/raw_final_compliance/final_to_best_gap_percent`
  - `/home/again/projects/recover/debug/utilities/fiber_finalize_iteration_outputs.m:348`

### 3.5 历史语义长度差（诊断风险）
- 观测：`history_count=29`，但 `accepted_source_len=28`。
- 解释：`compliance_history` 是状态序列（含初始状态），`accepted_source_history` 是步序列（每一步一次）。
- 风险：画图或审计时若不对齐索引，容易误判趋势。

### 3.6 环境日志噪声
- MATLAB `-batch` 下常见 WSL/UNC 路径清理警告（`run/rmpath` 相关）。
- 影响：不改变数值结果，但影响排障可读性。

## 4. 当前影响评估
- 功能可运行：是。
- 数值稳定性：中等风险（靠 rollback 保最终结果，过程末态可能退化）。
- 调参敏感性：高（守护阈值从 `0` 到 `0.01` 会显著改变轨迹）。
- 可观测性：中（已新增诊断字段，但历史口径需严格区分状态/步）。

## 5. 建议处理顺序
1. 先固定“默认配置策略”：明确 `default` 是否允许 1% 恶化（当前已放宽，属于行为改变）。
2. 给 `theta_only` 增加后段收敛约束（例如随迭代收紧 `current_state_tol`），避免最优后持续漂移。
3. 单独做 `HJ` 接管率诊断（持续记录 `theta_only/reinit/hj` 占比），确保几何步不是长期失效。
4. 在收敛图工具中强制区分“状态点”与“步事件”，避免索引错位。
5. 后续再处理 UNC 日志噪声（低优先级）。

## 6. 备注
- 本文档仅总结“改完后暴露的问题和现状”，不包含新的代码改动。
- 统计日期：2026-03-31。
