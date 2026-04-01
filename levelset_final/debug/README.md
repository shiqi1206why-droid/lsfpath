# Continuous Fiber Level-Set Project Handoff

## 1. What This Project Is

This is the maintained MATLAB project for optimizing continuous fiber paths inside a topology-optimized structure using a level-set field.

- Maintained working directory: `/home/again/projects/recover/debug`
- Main optimization entry: `fiber_levelset(config_name)`
- Topology source input: `topo_result.mat`
- Fiber path representation: zero level set of `lsf`
- Fiber direction representation: tangent direction of `lsf` contours
- Core objective: reduce compliance while keeping fiber paths continuous and restricted to the material domain

This README is meant as the primary handoff document for the next agent. It is more current than `PROJECT_MAP.md`.

## 2. Current Project State

### Mainline status

- The current maintained mainline is `fiber_levelset.m` plus the module directories:
  - `config/`
  - `initialization/`
  - `core_computation/`
  - `level_set_evolution/`
  - `utilities/`
  - `visualization/`
  - `postprocess/`
  - `tests/`
- Historical files exist but are not maintained:
  - `fiber_levelset_backup.m`
  - `fiber_levelset_v0.m`
  - `过往测试版本/`

### Current refactor status

The current mainline no longer treats the old local `dC/dphi` approximation as the only optimization path.

- `fiber_levelset(config_name)` is now only a wrapper entry.
- `run_fiber_levelset_pipeline` is now a 4-stage orchestrator:
  - `fiber_prepare_runtime_context`
  - `fiber_load_and_initialize_problem`
  - `fiber_run_optimization_iterations`
  - `fiber_finalize_pipeline_results`
- The heavy iteration logic was moved from `fiber_levelset.m` into `level_set_evolution/fiber_run_optimization_iterations.m`.
- Mid-level iteration submodules were added:
  - `fiber_run_hj_backtracking_step` (guard/backtracking logic)
  - `fiber_select_candidate_state` (theta/HJ/reinit candidate selection)
  - `fiber_record_history_entry` and `fiber_finalize_iteration_outputs` (history/diagnostics consolidation)
  - `fiber_log_iteration_snapshot` (iteration snapshot wrapper)
  - `fiber_log_velocity_path_diagnostics` / `fiber_log_sensitivity_velocity_summary` (iteration diagnostics wrappers)
- HJ RHS now supports explicit mode switch in params:
  - `params.levelset.hj_rhs_mode = legacy | indexed | vectorized_first_order`
  - default mode is `legacy` (strict-equivalence baseline); `indexed`/`vectorized_first_order` are opt-in
- `fiber_levelset.m` is now a compact flow controller.
- P1/P2治理（等价优先）已落地到主线，重点包括：
  - 候选语义拆分与可追踪：`hj_raw` / `hj_local_reinit` / `global_reinit` 分离，并保留 `accepted_source_history` 与 `accepted_source_detail_history`
  - 灵敏度聚合链式法则去冗余（在不改变数学表达和数值防护阈值前提下）
  - 迭代日志触发点包装化（仅搬迁输出，不改文案/时机）
  - 后处理拆分为职责函数：`extract_path_contours` / `smooth_path_candidates` / `select_best_smooth_path` / `write_path_to_file`
- 新的 `theta` 过渡链已统一为共享前向核心：
  - `phi -> theta_raw -> z_smooth -> theta_target -> theta_next`
  - `advance_theta_state` 与所有候选评估现在共享同一份前向离散实现
- 梯度链现支持三种模式：
  - `params.gradient.chain_mode = legacy | shadow | exact`
  - 当前配置默认：`default=legacy`、`fast=legacy`、`precise=exact`、`debug=shadow`
  - `shadow` 会并行计算 legacy / exact 并记录比较诊断，但仍用 legacy 驱动优化
- `differentiate_theta_transition_exact` 已加入主线：
  - 输入 `dC/dtheta_next` 和显式 `theta_transition_cache`
  - 沿真实离散链反传到 `dC/dphi`
- `aggregate_node_sensitivity` 已分为两条职责：
  - legacy 路径保留旧的局部近似装配
  - exact 路径只负责把 pullback 贡献装配到 `lsf` 节点
- exact 链离散模板已对齐为 Q4 四角节点（与 FE/HJ 一致）：
  - `compute_fiber_angles_from_lsf` 从中心差分切换为双线性 Q4 梯度
  - `differentiate_theta_transition_exact`/`assemble_exact_pullback` 改为四角节点贡献与装配
- exact 梯度到 HJ 速度新增了显式适配层：
  - `build_velocity_exact` 先做 `dC/dphi -> V_n` 适配（按 `|grad phi|` 归一）再进入统一速度整形
  - 2026-04-01 已修复 `build_velocity_exact` 的符号方向；`strict fast exact` 和 `test_exact_hj_descent_smoke` 现已通过
- `theta_only` 已新增后段防漂移机制：
  - 容差 schedule：前期放宽、中期线性收紧、后期严格不恶化
  - 连续接管熔断：`base.opt.theta_only_fuse_limit = 8`
- 候选质量诊断历史已接入结果：
  - `theta_only_vs_current_history`
  - `hj_raw_vs_theta_only_history`
  - `reinit_vs_theta_only_history`
- 优化态与评估态已经显式分离：
  - 评估态：`current_state` / `theta_only_state` / `hj_state` / `reinit_state`
  - 优化态：围绕 `theta_only_state.theta_transition_cache` 的真实离散导数链
- 制造约束附加项已接入优化态梯度入口：
  - `|grad phi|-1`
  - curvature proxy
  - gap/overlap proxy
  - 默认权重为 `0`，因此不会在默认配置下改变评估态行为

### Refactor equivalence check

Refactor equivalence was verified and passed at two levels:

- Historical refactor baseline:
  - `/home/again/projects/recover/debug/refactor_artifacts/refactor_comparison_report.txt`
- Latest P1/P2 strict equivalence (2026-03-30 14:12:44):
  - `/home/again/projects/recover/debug/refactor_artifacts/p1p2_equivalence_20260330_141244.txt`
  - `/home/again/projects/recover/debug/refactor_artifacts/p1p2_equivalence_20260330_141244.mat`

Latest strict report summary (`pass_all = 1`):

- `final_compliance = 1`
- `final_FCS = 1`
- `final_iter = 1`
- `compliance_history = 1`
- `FCS_history = 1`
- `lsf = 1`
- `theta_e = 1`
- `accepted_source_history = 1`
- key guard counters (`accepted/rejected/reinit-skip`) = 1
- `final_compliance_rel_diff = 0`
- `final_FCS_rel_diff = 0`
- `compliance_history_max_abs_diff = 0`
- `FCS_history_max_abs_diff = 0`
- `lsf_inf_diff = 0`

## 3. Latest Verified Runs

### Latest saved default rerun

- Artifact: `/home/again/projects/recover/debug/refactor_artifacts/rerun_default_20260329_004613.mat`
- Summary: `/home/again/projects/recover/debug/refactor_artifacts/rerun_default_20260329_004613_summary.txt`
- Image: `/home/again/projects/recover/debug/refactor_artifacts/rerun_default_20260329_004613_final_lsf.png`

Latest recorded metrics:

- `final_compliance = 8.5231331993422949e-07`
- `final_FCS = 8.7652681890600104e-01`
- `final_iter = 38`

### Latest HJ RHS benchmark (2026-03-30 14:13:33)

- Report: `/home/again/projects/recover/debug/refactor_artifacts/hj_rhs_benchmark_20260330_141333.txt`
- Data: `/home/again/projects/recover/debug/refactor_artifacts/hj_rhs_benchmark_20260330_141333.mat`
- `legacy_order2 mean_time = 0.058115 sec`
- `indexed_order2 mean_time = 0.049895 sec`
- `legacy_order1 mean_time = 0.017737 sec`
- `vectorized_first_order1 mean_time = 0.004615 sec`
- `eq_legacy_vs_indexed_order2 = 1`
- `eq_legacy_vs_vectorized_order1 = 1`
- `speedup_indexed_vs_legacy_order2 = 1.164736`
- `speedup_vectorized_vs_legacy_order1 = 3.843627`

### Latest paired default+precise rerun (2026-03-31 23:51:29)

- Artifact dir:
  - `/home/again/projects/recover/debug/refactor_artifacts/rerun_default_precise_compare_20260331_235129`
- Stitched comparison image:
  - `/home/again/projects/recover/debug/refactor_artifacts/rerun_default_precise_compare_20260331_235129/default_vs_precise_stitched.png`
- Summary files:
  - `/home/again/projects/recover/debug/refactor_artifacts/rerun_default_precise_compare_20260331_235129/summary.txt`
  - `/home/again/projects/recover/debug/refactor_artifacts/rerun_default_precise_compare_20260331_235129/summary.json`

Recorded metrics from this paired rerun:

- `default` (`chain_mode=legacy`)
  - `initial_compliance = 3.777978412500e-06`
  - `final_compliance = 3.771308144606e-06`
  - `best_compliance = 3.771306672634e-06`
  - `final_FCS = 0.769516728625`
  - `final_iter = 82` (`executed_iter = 82`)
  - `accepted_steps = 81`, `rejected_steps = 0`
- `precise` (`chain_mode=shadow`)
  - `initial_compliance = 3.875270585552e-06`
  - `final_compliance = 3.869509360502e-06`
  - `best_compliance = 3.869511289285e-06`
  - `final_FCS = 0.808019118428`
  - `final_iter = 99` (`executed_iter = 99`)
  - `accepted_steps = 0`, `rejected_steps = 98`

Interpretation for handoff:

- `default (legacy)` is still the stable optimization baseline with normal accepted-step progression.
- `precise (shadow)` currently behaves as a guard-dominated validation run (plateau + mostly rejected HJ updates), which matches the transition-stage rollout intent but is not yet a strong optimization driver.

### Latest exact-chain repair verification (2026-04-01)

- Stage 1 strict exact fast regression:
  - `/home/again/projects/recover/debug/refactor_artifacts/stage1_fix_exact_true_20260401_162840/results_fast_exact_stage1.mat`
  - `chain_mode = exact`
  - `accepted_steps = 38`, `rejected_steps = 6`
  - `final_iter = 45`
  - `final_compliance = 3.545347490606e-06`
- Stage 2 exact HJ smoke:
  - `/home/again/projects/recover/debug/refactor_artifacts/stage2_fix_20260401_163440/test_exact_hj_descent_smoke.log`
  - `C_current = 3.777978e-06`
  - `C_hj = 3.761182e-06`
  - `delta = -1.679682e-08`
  - `dt = 8.713035e-03`, `dt_angle = 8.713035e-03`
- Stage 2 strict exact fast regression:
  - `/home/again/projects/recover/debug/refactor_artifacts/stage2_fix_20260401_163440_strict/results_fast_exact_stage2_regression.mat`
  - `chain_mode = exact`
  - `accepted_steps = 38`, `rejected_steps = 6`
  - `final_iter = 45`
  - `final_compliance = 3.545347490606e-06`

Interpretation for handoff:

- Step1 / Step4 的直接失配点已经修复：`exact -> HJ` 的速度符号方向正确，测试口径也已对齐主循环。
- 这次修复证明 `exact` 链在严格 `fast` 条件下不再是“零接受步”。
- 但这不代表 exact 全流程已经完成标定；它只是从“明显错误”进入“可运行但仍需调参”的状态。

### Latest precise exact run (2026-04-01 15:27:43)

- Artifact dir:
  - `/home/again/projects/recover/debug/refactor_artifacts/precise_exact_rerun_20260401_152743`
- Files:
  - `/home/again/projects/recover/debug/refactor_artifacts/precise_exact_rerun_20260401_152743/results_precise_exact.mat`
  - `/home/again/projects/recover/debug/refactor_artifacts/precise_exact_rerun_20260401_152743/precise_exact_log.txt`
  - `/home/again/projects/recover/debug/refactor_artifacts/precise_exact_rerun_20260401_152743/precise_exact_like_example.png`

Recorded metrics:

- `chain_mode = exact`
- `final_compliance = 3.869509360502e-06`
- `final_FCS = 0.8080`
- `final_iter = 99`
- `accepted_steps = 0`, `rejected_steps = 98`
- `theta_only accepted = 86`, `theta_only rejected = 12`

Interpretation for handoff:

- `precise` 当前虽然已经配置为 `exact`，但整轮行为仍然是明显的 guard-dominated 模式。
- 也就是说，Step1 / Step4 修复解决了“方向错误”和“测试口径错误”，还没有解决 `exact` 在长程运行中的 acceptance 标定问题。

### Tests recently confirmed

The following were explicitly run and passed in this environment:

- `fiber_levelset('debug')` smoke run
- `tests/test_material_domain_restriction.m`
- `tests/test_second_order_hj_material_mask.m`
- `tests/test_postprocess_export_smoke.m`
- `tests/step_checks/test_stepA_acceptance_self_consistency.m`
- `tests/comparison/run_p1p2_strict_equivalence_check.m`
- `tests/comparison/run_hj_rhs_mode_benchmark.m`
- `tests/formula_audit/test_dE_dtheta_fd_audit.m`
- `tests/formula_audit/test_theta_transition_exact_fd_audit.m`
- `tests/formula_audit/test_theta_transition_directional_derivative_audit.m`
- `tests/formula_audit/test_exact_direct_phi_descent_smoke.m`
- `tests/formula_audit/test_exact_hj_descent_smoke.m`
- `tests/test_exact_chain_default_smoke.m`
- `tests/comparison/run_default_precise_rerun_and_stitch.m` (2026-03-31 23:51:29 run tag: `rerun_default_precise_compare_20260331_235129`)
- `strict fast exact` regression after Step1/Step4 fix (2026-04-01): `accepted_steps = 38`

Additional scripts in exact-chain rollout:

- `tests/formula_audit/test_theta_transition_forward_consistency.m`
- `tests/formula_audit/test_theta_transition_cache_selfcheck.m`
- `tests/formula_audit/test_theta_transition_exact_fd_audit.m`
- `tests/formula_audit/test_theta_transition_directional_derivative_audit.m`
- `tests/formula_audit/test_theta_transition_limiter_diagnostics.m`
- `tests/formula_audit/test_manufacturing_penalty_directional_derivatives.m`
- `tests/test_exact_chain_default_smoke.m`

## 4. Required Input / Output Contracts

### Required topology input

`fiber_levelset` expects `topo_result.mat` in the project root.

Minimum required field:

- `struc`

Optional but expected:

- `nelx`
- `nely`
- legacy `lsf` may exist but the current fiber optimizer rebuilds its own `lsf`

### Results struct

The returned `results` struct is the main programmatic output. Important fields include:

- `lsf`
- `theta_e`
- `compliance_history`
- `FCS_history`
- `final_compliance`
- `final_FCS`
- `material_mask_core`
- `material_mask_full`
- `path_quality_raw`
- `path_quality_history`
- `interface_diagnostics`
- `gradient_chain_history`
- `params`
- `init_info`

### Runtime path contract

The maintained mainline no longer assumes the current working directory is the project root.

- `params.runtime.project_root`: absolute project-root path
- `params.runtime.paths`: canonical absolute paths for topology, checkpoints, visualization, baseline, refactor, and tests

Maintained entrypoints should prefer `params.runtime.paths` or `build_project_paths(project_root)` and should not infer project locations from `pwd`.

### Material-mask semantics

- `material_mask_core`: `[nely, nelx]`
- `material_mask_full`: same field expanded to the `lsf` grid including ghost cells

All recent material-domain restrictions are supposed to use the cleaned mask as the single source of truth.

## 5. High-Level Workflow

The maintained optimizer does this:

1. Load parameters from `config/get_fiber_optimization_params.m`
2. Load `topo_result.mat`
3. Clean topology into a material mask
4. Build a boundary-offset initial level set
5. Evaluate current FE state
6. Per iteration:
   - evaluate `theta_only`
   - compute FE sensitivity
   - build optimization-state gradient chain (`legacy` / `shadow` / `exact`)
   - aggregate `dC/dtheta_next` or legacy `dC/dtheta` to node sensitivity
   - optionally add manufacturing penalty gradients
   - build a masked narrow-band velocity field
   - advance `lsf` using HJ update
   - optionally do local/global reinitialization
   - compare candidates under guard conditions
   - keep the best acceptable state
7. Apply best-state rollback if needed
8. Visualize and return `results`

This is **not** a plain gradient-descent loop. It is intentionally conservative and guard-heavy.

## 6. Directory Map

### Core project code

- `config/`
  - run modes and material database
- `initialization/`
  - material mask cleanup
  - boundary reconstruction
  - signed-distance initialization
- `core_computation/`
  - FE analysis
  - constitutive law
  - angle update
  - candidate-state evaluation
  - HJ update
- `level_set_evolution/`
  - sensitivity aggregation
  - velocity field
  - timestep control
  - reinit trigger logic
- `utilities/`
  - mask normalization
  - Neumann boundary handling
  - contour parsing
  - geometry extraction
  - path-quality assembly
  - logging/checkpoint helpers
- `visualization/`
  - main article plots
  - initialization diagnostics
- `postprocess/`
  - printable-path extraction and smoothing
- `tests/`
  - regression, formula audit, step checks, mechanism tests

### Important non-code directories

- `baseline_artifacts/`
  - many historical run outputs; large and noisy
- `refactor_artifacts/`
  - current handoff-relevant refactor/run evidence
- `checkpoints/`
  - optional checkpoint output
- `visualization_artifacts/`
  - plot outputs from earlier runs

## 7. Main Entry Points

### Main optimizer

- `fiber_levelset('default')`
- `fiber_levelset('fast')`
- `fiber_levelset('precise')`
- `fiber_levelset('debug')`

### Legacy topology generator

- `levelset_top.m`

This is legacy topology-generation code that writes `topo_result.mat`. It is not the modern fiber-path optimizer.

### Baseline scripts

- `tests/baseline/run_fast_baseline.m`
- `tmp_run_default_export.m`

### Strict equivalence check

- `tests/comparison/run_refactor_strict_equivalence_check.m`
  - Compares current run vs latest `refactor_artifacts/pre_default_*.mat`
  - Uses strict `isequaln` checks on key scalar fields, histories, state arrays, and guard counters
  - Writes report to `refactor_artifacts/strict_refactor_comparison_*.{mat,txt}`
- `tests/comparison/run_p1p2_strict_equivalence_check.m`
  - P1/P2治理后的严格等价入口（同样基于 `pre_default_*.mat`）
  - 报告写入 `refactor_artifacts/p1p2_equivalence_*.{mat,txt}`
- `tests/comparison/run_hj_rhs_mode_benchmark.m`
  - HJ RHS 模式性能对比（legacy/indexed/vectorized_first_order）
  - 报告写入 `refactor_artifacts/hj_rhs_benchmark_*.{mat,txt}`
- `tests/formula_audit/test_theta_transition_exact_fd_audit.m`
  - exact-chain 节点中心差分审计（当前已通过）
- `tests/formula_audit/test_theta_transition_directional_derivative_audit.m`
  - exact-chain 方向导数审计（当前已通过，方向采样支撑已收紧到有效梯度节点）
- `tests/formula_audit/test_exact_direct_phi_descent_smoke.m`
  - exact 梯度直接 `phi` 微步下降审计（当前已通过）
- `tests/formula_audit/test_exact_hj_descent_smoke.m`
  - exact 梯度驱动 HJ 微步下降审计（当前已通过）
- `tests/test_exact_chain_default_smoke.m`
  - 过渡期配置 smoke（校验默认链模式为 legacy 并跑通）

Important user preference:

- Do **not** create or update baseline artifacts unless the user explicitly asks.

## 8. Run Instructions

### From WSL shell in this environment

Use the Windows MATLAB binary against the WSL repo path:

```bash
PROJECT_ROOT="/home/again/projects/recover/debug"
MATLAB_WIN="/mnt/d/Program Files/matlab/R2024a/bin/matlab.exe"
MATLAB_CD_WIN=$(wslpath -w "$PROJECT_ROOT")
"$MATLAB_WIN" -batch "cd('$MATLAB_CD_WIN'); addpath(genpath(pwd), '-begin'); results = fiber_levelset('default');"
```

Or run the helper script:

```bash
/home/again/projects/recover/debug/run_matlab_default_from_wsl.sh
```

### Common run modes

- `default`: main acceptance run
- `fast`: quick iteration / lightweight checks
- `precise`: slower, stricter
- `debug`: high diagnostics, useful for smoke checks

## 9. Known Issues / Things The Next Agent Must Not Miss

### 1) Unit semantics are still deferred (explicitly out-of-scope in this round)

- `dx`, `dy`, `h` unit convention is not unified in code comments/log text.
- Any spacing or geometry interpretation still needs a dedicated unit-semantics cleanup task before metric-level optimization.

### 2) Spacing quality metric coherence still pending

- There is still a known mismatch risk between:
  - `path_quality_raw.parallel_spacing_measured`
  - spacing inferred from `path_quality_raw.parallel_spacing_error_percent`
- Do not claim spacing improvement until this metric pair is unified.

### 3) Complexity hotspot remains in iteration engine

- `fiber_levelset.m` is already orchestration-only.
- Main complexity concentration is still:
  - `level_set_evolution/fiber_run_optimization_iterations.m`

### 4) Manufacturing penalties are currently optimization-side only

- They are added to the optimization gradient path, not to candidate-performance evaluation.
- Current `curvature` / `gap_overlap` terms are rollout-stage proxies with default weight `0`.

### 5) Historical outputs remain noisy

- `baseline_artifacts/` and older refactor artifacts contain many legacy runs.
- Mainline maintenance target remains `/home/again/projects/recover/debug`.

### 6) Exact-chain rollout remains in transition (not default yet)

- `default` / `fast` remain the stable `legacy` baseline; `precise` is now configured as `exact`, but this should still be treated as a rollout-stage mode rather than a fully tuned production baseline.
- Step1 / Step4 failure has been fixed: `build_velocity_exact` sign and `test_exact_hj_descent_smoke` test path are now aligned with the main loop.
- Strict `fast exact` now passes with `accepted_steps = 38`, so exact is no longer failing at the “zero accepted steps under strict smoke” level.
- The remaining main risk is acceptance composition:
  - accepted steps are still dominated by `theta_only`
  - `HJ` and especially `reinit` are still frequently blocked by `next/current` guards
- Full-run acceptance behavior under `exact` still needs dedicated retuning of guard thresholds, velocity shaping, and possibly exact-only stabilization in non-smoke settings.
- Directional-derivative audits are now passing with differentiable-support sampling; this support definition must be preserved, otherwise audits can be falsely negative.
- Latest precise exact run on 2026-04-01 still confirms a guard-dominated long run (`accepted=0`, `rejected=98`) with most progress coming from `theta_only`, not from accepted HJ / reinit evolution.

## 10. User Preferences That Matter

These are important behavioral preferences from the current user and should be followed:

- Use absolute paths when referencing files.
- Keep logic correct end-to-end; do not hide incorrect process with cosmetic postprocessing.
- Prefer clear, modular MATLAB structure.
- Main function should remain a flow controller, not a utility dump.
- Do not run or generate baseline artifacts unless explicitly requested.
- When changing algorithms, preserve traceable logic and validate results rather than trusting plots.

## 11. Recommended Next Actions

If the next task is **continued refactor**:

1. Keep `fiber_levelset.m` as a pure orchestration file.
2. Continue slimming `fiber_run_optimization_iterations.m` (candidate/reinit/diagnostic branches remain the largest blocks).
3. Preserve the optimization-state vs evaluation-state separation.
4. Re-run `fast` smoke + exact-chain audits after each structural change.

If the next task is **spacing control / path-quality improvement**:

1. First unify unit semantics (`mm` vs `m`).
2. Then unify spacing measurement logic.
3. Only after that change target spacing or initialization offset.
4. Re-run `default` and compare `path_quality_raw`, not just plots.

If the next task is **numerical quality**:

1. Re-run the exact-chain audit suite and shadow diagnostics on real runs.
2. Then evaluate whether to switch default `hj_rhs_mode` from `legacy` to `indexed`.
3. Then focus on `update_levelset_HJ.m`, `fmm_reinitialize.m`, and `compute_lsf_path_quality.m` for deeper numerical tuning.

## 12. Useful Companion Files

- Current project map: `/home/again/projects/recover/debug/PROJECT_MAP.md`
- Refactor notes: `/home/again/projects/recover/debug/refactor_artifacts/refactor_change_notes.md`
- Refactor structure notes: `/home/again/projects/recover/debug/refactor_artifacts/refactor_structure_analysis.md`
- Refactor equivalence report: `/home/again/projects/recover/debug/refactor_artifacts/refactor_comparison_report.txt`
- Latest P1/P2 equivalence report: `/home/again/projects/recover/debug/refactor_artifacts/p1p2_equivalence_20260330_141244.txt`
- Latest HJ RHS benchmark report: `/home/again/projects/recover/debug/refactor_artifacts/hj_rhs_benchmark_20260330_141333.txt`

If this README and `PROJECT_MAP.md` disagree, trust this README first.
