# Debug 笔记（2025-10-22 22:16）

- 目标：复现“水平集的悬臂梁纤维路径规划”论文方案，工程化实现（模块化、日志/诊断、可视化、基准）。
- 入口：`fiber_levelset.m`（支持 `default/fast/precise/debug`）。依赖 `topo_result.mat` 中的 `struc` 初始拓扑。

**要点映射（代码对应）**
- 初始化与锚定：`φ_target = φ_b + Δφ`，FMM 重初始化保持 |∇φ|≈1。
  - `fiber_levelset.m:145` 目标场；`initialization/construct_boundary_offset_levelset_with_parallel.m:1`、`fmm_reinitialize.m:1` 重初始化。
- 方向场与平滑：θ = π/2 + atan2(φ_y,φ_x)，范围 [0,π)；矢量化拉普拉斯平滑；每步限幅 Δθ_max。
  - `core_computation/compute_fiber_angles_from_lsf.m:1`，`level_set_evolution/angle_smooth_vectorized.m:1`，`fiber_levelset.m:224`。
- 力学与灵敏度：各向异性 Q→C(θ) 与 dC/dθ；S=C\I 稳定化；伴随灵敏度并映射到节点。
  - `core_computation/compute_element_strain_stress.m:1`，`core_computation/compute_sensitivity_adjoint.m:1`，`level_set_evolution/aggregate_node_sensitivity.m:1`。
- 速度场与演化：v_shape 去均值 + v_fid 锚定 + 曲率正则；HJ 上风更新；CFL 自适应步长；按策略重初始化。
  - `level_set_evolution/build_velocity_field.m:1`，`core_computation/update_levelset_HJ.m:1`，`level_set_evolution/compute_adaptive_timestep.m:1`，`level_set_evolution/should_reinitialize.m:1`。
- 诊断/指标/可视化：柔度 U'F、一致性校验、FCS 连续性、可视化与日志。
  - `fiber_levelset.m:266`，`core_computation/compute_fiber_continuity.m:1`，`visualization/visualize_results_article.m:1`，`utilities/diag_report.m:1`。

**发现的问题（精简）**
- 重初始化未传 `material_mask`：`fiber_levelset.m:470` 使用 `[]`，FMM 将在全域传播，可能削弱边界一致性。
- λ_fid 被固定值覆盖：`fiber_levelset.m:321` 直接设为 50，覆盖了前面的 dev95 自适应估计。
- CFL 因子保守：`level_set_evolution/compute_adaptive_timestep.m:1` 固定 `cfl_factor=0.08`，与要点建议（≈0.3）不一致。
- 曲率正则与诊断阈值硬编码：`fiber_levelset.m:323`（gamma_curv）、`build_velocity_field.m:1`（v_target）、`fiber_levelset.m:497`（验收阈值）。

**最小修复建议**
- FMM 传入结构掩码：`fiber_levelset.m:470` → `fmm_reinitialize(lsf, dx, dy, zero_mask_dynamic, material_mask);`
- 去除或加开关控制固定 λ_fid，保留 dev95 动态法并做上下限夹取。
- 将 `cfl_factor`、`v_target(early/late)`、`gamma_curv`、验收阈值参数化至 `config/get_fiber_optimization_params.m:1`。
- 在 `build_velocity_field.m:1` 为 v_fid 增加距离权重 `w_dist(|φ|, h)`（如线性衰减或 exp 衰减），仅窄带内生效。

**快速验证**
- 角度平滑与性能：运行 `tests/test_angle_smooth_vectorized.m:1`。
- 配置冒烟：运行 `tests/test_all_configurations.m:1`（建议 `debug`）。

