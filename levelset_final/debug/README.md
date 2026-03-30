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

The recent refactor did **not** change numerical behavior.

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

### Tests recently confirmed

The following were explicitly run and passed in this environment:

- `fiber_levelset('debug')` smoke run
- `tests/test_material_domain_restriction.m`
- `tests/test_second_order_hj_material_mask.m`
- `tests/test_postprocess_export_smoke.m`
- `tests/step_checks/test_stepA_acceptance_self_consistency.m`
- `tests/comparison/run_p1p2_strict_equivalence_check.m`
- `tests/comparison/run_hj_rhs_mode_benchmark.m`

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
   - aggregate `dC/dtheta` to node sensitivity
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

### 4) Historical outputs remain noisy

- `baseline_artifacts/` and older refactor artifacts contain many legacy runs.
- Mainline maintenance target remains `/home/again/projects/recover/debug`.

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
3. Keep strict-equivalence gating as hard constraint (`run_p1p2_strict_equivalence_check`).
4. Re-run `debug` smoke + material-domain tests after each structural change.

If the next task is **spacing control / path-quality improvement**:

1. First unify unit semantics (`mm` vs `m`).
2. Then unify spacing measurement logic.
3. Only after that change target spacing or initialization offset.
4. Re-run `default` and compare `path_quality_raw`, not just plots.

If the next task is **numerical quality**:

1. Evaluate whether to switch default `hj_rhs_mode` from `legacy` to `indexed` (only after strict-equivalence + regression pass under target configs).
2. If pursuing more speed, keep `vectorized_first_order` opt-in and separately validate its use-case envelope.
3. Then focus on `update_levelset_HJ.m`, `fmm_reinitialize.m`, and `compute_lsf_path_quality.m` for deeper numerical tuning.

## 12. Useful Companion Files

- Current project map: `/home/again/projects/recover/debug/PROJECT_MAP.md`
- Refactor notes: `/home/again/projects/recover/debug/refactor_artifacts/refactor_change_notes.md`
- Refactor structure notes: `/home/again/projects/recover/debug/refactor_artifacts/refactor_structure_analysis.md`
- Refactor equivalence report: `/home/again/projects/recover/debug/refactor_artifacts/refactor_comparison_report.txt`
- Latest P1/P2 equivalence report: `/home/again/projects/recover/debug/refactor_artifacts/p1p2_equivalence_20260330_141244.txt`
- Latest HJ RHS benchmark report: `/home/again/projects/recover/debug/refactor_artifacts/hj_rhs_benchmark_20260330_141333.txt`

If this README and `PROJECT_MAP.md` disagree, trust this README first.
