# MATLAB 工程重构分析与模块分类（主线版本）

## 1) 主函数职责识别

- 主入口：`/mnt/e/codex project/recover/debug/fiber_levelset.m`
- 主职责（流程总控）：
  - 加载参数与拓扑输入
  - 构建材料域掩膜与初始 level set
  - 迭代调度：灵敏度 -> 速度场 -> HJ 更新 -> 重初始化 -> 接受准则
  - 记录诊断、收敛判断、回滚与结果组装
  - 调用可视化与输出 `results`

## 2) 辅助函数职责识别

### 初始化模块（`initialization/`）
- `clean_material_mask`：清洗拓扑掩膜并输出掩膜诊断
- `construct_boundary_offset_levelset_with_parallel`：构造偏移初始化场
- `reconstruct_subcell_boundary_geometry` / `compute_signed_distance_to_boundary`：子像素边界重建与签名距离
- `compute_boundary_offset_stats` / `verify_path_spacing`：初始化质量与间距校验

### 演化模块（`level_set_evolution/`）
- `aggregate_node_sensitivity`：单元灵敏度到节点的链式聚合
- `build_velocity_field`：窄带速度构建（含掩膜传播约束）
- `compute_adaptive_timestep`：CFL 自适应时间步
- `compute_zero_mask_from_lsf` / `should_reinitialize`：重初始化触发依据

### 核心求解模块（`core_computation/` + 根目录重初始化）
- `update_levelset_HJ`：HJ 数值更新
- `fmm_reinitialize`：重初始化（几何优先 + FMM 回退）
- `FE_analysis_cantilever` / `compute_sensitivity_adjoint`：FE 与灵敏度主链路
- `evaluate_candidate_state`：候选状态一致性求值（compliance/FCS 等）

### 后处理与可视化模块
- `postprocess/export_printable_paths_from_lsf`：路径导出
- `postprocess/refine_lsf_for_printability`：打印友好补充优化
- `visualization/visualize_results_article`：论文主图输出

### 通用工具模块（`utilities/`）
- 等值线/几何提取、距离重建、质量指标计算、日志与检查点
- 本次新增：
  - `normalize_mask_to_lsf_grid`：统一核心/全域掩膜归一
  - `apply_neumann_boundary`：统一 Neumann ghost 边界扩展

## 3) 功能分类后的推荐分层

- `fiber_levelset` 仅承担编排，不内嵌通用数值细节
- `core_computation` 专注 FE、候选评估、HJ 更新
- `level_set_evolution` 专注速度场与时间推进决策
- `initialization` 专注初始几何/距离场构建
- `utilities` 专注纯工具能力（可复用、无业务语义）
- `postprocess` 仅做制造与导出补充，不反向影响主优化逻辑
- `visualization` 仅消费结果，不改状态

## 4) 推荐目录组织（兼容现有主线）

当前目录保持不变，建议后续迭代按以下方式细化（本轮不做激进迁移）：

- `core_computation/`
  - `fe/`（FE 装配、本构、应变应力）
  - `objective/`（候选状态求值、柔度/FCS）
  - `hj/`（HJ 更新）
- `level_set_evolution/`
  - `velocity/`（速度构建、灵敏度聚合）
  - `reinit/`（触发策略）
- `utilities/`
  - `mask/`（掩膜归一、边界扩展）
  - `geometry/`（轮廓提取、几何解析）
  - `metrics/`（路径质量指标）

## 5) 本轮重构落地点（不改变功能）

- 抽出并统一公共逻辑：
  - `utilities/normalize_mask_to_lsf_grid.m`
  - `utilities/apply_neumann_boundary.m`
- 调整调用方复用公共能力：
  - `core_computation/update_levelset_HJ.m`
  - `level_set_evolution/build_velocity_field.m`
  - `fmm_reinitialize.m`
- 在 `fiber_levelset.m` 内部将“历史收尾裁剪 / 诊断组装 / 结果结构体组装”改为独立局部函数，主流程可读性提升，算法逻辑不变。
