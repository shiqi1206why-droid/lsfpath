# Project Map

## Purpose

This project optimizes continuous fiber paths inside a topology-optimized structure using a level set function `phi`.

- Objective: reduce compliance while keeping fiber paths smooth and continuous.
- Geometry representation: the zero level set of `lsf`.
- Fiber direction: tangent direction of `lsf` contours.
- Main entry: `fiber_levelset(config_name)`.

## Main Execution Flow

`fiber_levelset.m` is now a compact 4-stage orchestrator:

1. `fiber_prepare_runtime_context`
   - resolve project root/paths
   - load + validate configuration
   - build runtime context/options
2. `fiber_load_and_initialize_problem`
   - load `topo_result.mat`
   - clean material mask
   - build boundary-offset initial level set
   - preallocate histories and initialize first FE state
3. `fiber_run_optimization_iterations`
   - run the full iterative optimizer:
   - update fiber angles from `lsf`
   - run FE analysis
   - compute compliance and FCS
   - build optimization-state gradient chain (`legacy` / `shadow` / `exact`)
   - compute angle sensitivity `dC/dtheta_next`
   - aggregate to node sensitivity `dC/dphi`
   - optionally add optimization-side manufacturing penalties
   - build a narrow-band velocity field
   - advance `lsf` with HJ update
   - optionally reinitialize with FMM
   - compare `theta-only`, `HJ`, and `reinit` candidates
   - keep the best acceptable candidate
4. `fiber_finalize_pipeline_results`
   - visualize
   - emit summary logs
   - build and return the final `results` struct

## Module Layout

### Config

- `config/get_fiber_optimization_params.m`
  - Defines `default`, `fast`, `precise`, `debug`.
  - Central place for optimization controls, reinit logic, logging, smoothing, and velocity scaling.
- `config/get_material_params.m`
  - Material database.

### Initialization

- `initialization/clean_material_mask.m`
  - Removes small disconnected regions from topology output.
- `initialization/construct_boundary_offset_levelset_with_parallel.m`
  - Builds the initial signed distance field and parallel paths from the cleaned mask.
- `initialization/compute_boundary_offset_stats.m`
  - Checks offset quality against the target distance.

### Core Computation

- `core_computation/advance_theta_state.m`
  - Thin wrapper around the shared theta-transition forward core.
- `core_computation/compute_theta_transition_forward.m`
  - Shared forward discretization for `phi -> theta_raw -> z_smooth -> theta_target -> theta_next`.
- `core_computation/compute_fiber_angles_from_lsf.m`
  - Computes element fiber angle from `grad(lsf)` and returns explicit gradient caches.
- `core_computation/FE_analysis_cantilever.m`
  - Assembles and solves the FE system on a cantilever domain.
- `core_computation/compute_sensitivity_adjoint.m`
  - Computes `dC/dtheta` using `-Ue' * dKe/dtheta * Ue`.
- `core_computation/differentiate_theta_transition_exact.m`
  - Pulls `dC/dtheta_next` back through the exact shared theta-transition discretization.
- `core_computation/evaluate_state_with_theta.m`
  - Packs FE solve, compliance, strain energy, and FCS for a given `(lsf, theta)`.
- `core_computation/evaluate_candidate_state.m`
  - One-step candidate evaluation on a trial `lsf`, including `theta_transition_cache`.
- `core_computation/element_stiffness.m`
  - Builds `Ke` and optionally `dKe/dtheta` from the orthotropic constitutive law.
- `core_computation/orthotropic_constitutive_matrix.m`
  - Builds `C`, `dC/dtheta`, and optionally `S`.

### Level Set Evolution

- `level_set_evolution/aggregate_node_sensitivity.m`
  - Legacy local-chain aggregation plus exact pullback assembly.
- `level_set_evolution/build_velocity_field.m`
  - Builds the narrow-band normal velocity from node sensitivity.
- `level_set_evolution/compute_gradient_chain_sensitivity.m`
  - Central selector for legacy/shadow/exact gradient-chain evaluation.
- `level_set_evolution/compute_manufacturing_penalty_gradient.m`
  - Adds optimization-side manufacturing penalty gradients.
- `level_set_evolution/compute_adaptive_timestep.m`
  - CFL-based robust time-step selection.
- `level_set_evolution/should_reinitialize.m`
  - Decides when FMM reinitialization should run.
- `level_set_evolution/fiber_run_optimization_iterations.m`
  - Main iterative optimization loop (current largest module).
- `level_set_evolution/fiber_select_candidate_state.m`
  - Unified theta/HJ/reinit candidate selection.
- `core_computation/update_levelset_HJ.m`
  - Explicit upwind HJ update.
  - `rhs_mode` switch supports `legacy` / `indexed` / `vectorized_first_order`.
- `fmm_reinitialize.m`
  - Rebuilds the distance field while preserving the zero set.

### Refactor Support Helpers

- `core_computation/fiber_run_hj_backtracking_step.m`
  - Encapsulates HJ step acceptance, guard checks, and backtracking.
- `utilities/fiber_record_history_entry.m`
  - Centralized per-iteration history recording.
- `utilities/fiber_finalize_iteration_outputs.m`
  - Post-loop history trimming, rollback selection, and diagnostics assembly.
- `utilities/fiber_log_iteration_snapshot.m`
  - Wrapper for periodic iteration logging output.
- `utilities/fiber_log_velocity_path_diagnostics.m`
  - Wrapper for gradient/deviation/raw-path diagnostics.
- `utilities/fiber_log_sensitivity_velocity_summary.m`
  - Wrapper for sensitivity/velocity/timestep summary diagnostics.

### Visualization and Diagnostics

- `visualization/visualize_results_article.m`
- `visualization/enhanced_visualization_check.m`
- `utilities/log_message.m`
- `utilities/diag_report.m`
- `utilities/save_checkpoint.m`

### Postprocess

- `postprocess/export_printable_paths_from_lsf.m`
- `postprocess/refine_lsf_for_printability.m`
- `postprocess/run_printable_path_smoothing.m`
- `postprocess/extract_path_contours.m`
  - Contour extraction + material-domain clipping.
- `postprocess/smooth_path_candidates.m`
  - Candidate smoothing generation and metric evaluation.
- `postprocess/select_best_smooth_path.m`
  - Candidate ranking/selection.
- `postprocess/write_path_to_file.m`
  - Unified CSV export writer.

These focus on extracting printable paths and smoothing them after optimization.

## Important State Variables

Inside `fiber_levelset`, the most important evolving quantities are:

- `lsf`: signed distance-like level set field with ghost cells.
- `theta_e`: element fiber angle field.
- `theta_target`: smoothed target angle field derived from `lsf`.
- `theta_transition_cache`: explicit cache for the shared theta-transition forward chain.
- `U`, `K`, `F`: FE displacement, stiffness, load.
- `compliance`: objective value.
- `FCS`: fiber continuity score.
- `best_state`: historical best state used for rollback and early stop.

There are now two different roles for state:

- Evaluation state
  - explicit `(lsf, theta)` candidate with FE results used for acceptance and rollback.
- Optimization state
  - the shared theta-transition cache used to compute the real discrete `dC/dphi`.

## Candidate Selection Logic

Each iteration may produce up to three candidate states:

- `theta_only`
  - same `lsf`, only refresh `theta`.
- `hj`
  - advance `lsf` by one HJ step.
- `reinit`
  - rebuild the HJ candidate using FMM if reinit is triggered.

Acceptance is guarded by:

- `acceptance_tol`
  - candidate should not be worse than the reference next state.
- `current_state_tol`
  - candidate should not be worse than the current state.
- optional reinit-specific current-state guard.
- optional historical best-state patience for early stop.

This means the optimizer is deliberately conservative. It is not a pure gradient descent loop.

## Gradient Chain Modes

`params.gradient.chain_mode` supports:

- `legacy`
  - old local `dC/dtheta -> dC/dphi` approximation
- `shadow`
  - compute both legacy and exact chain, but still optimize with legacy
- `exact`
  - optimize with the exact shared-discretization pullback

`params.gradient.limiter_mode` supports:

- `hard`
  - piecewise-exact subgradient for the rate limiter
- `soft_experiment`
  - experimental smooth pullback for optimization only

## Numerical Intent

The current implementation is shaped by a few strong design choices:

- Keep `theta` evolution rate-limited with `delta_theta_max`.
- Restrict updates to a narrow band around the interface.
- Avoid every-step hard reinitialization.
- Use FMM as a recovery candidate, not as an unconditional overwrite.
- Prefer objective safety over aggressive movement.
- Track and potentially roll back to the best historical state.

## Existing Test Strategy

The project already contains useful verification layers:

- `tests/test_all_configurations.m`
  - config system, material database, logging, parameter validation.
- `tests/test_angle_smooth_vectorized.m`
  - speed and numerical consistency of vectorized angle smoothing.
- `tests/test_aggregate_sensitivity_narrowband.m`
  - verifies zero sensitivity leakage outside the narrow band.
- `tests/step_checks/*.m`
  - acceptance logic and step-guard behavior.
- `tests/formula_audit/test_dE_dtheta_fd_audit.m`
  - finite-difference audit of `dC/dtheta`.
- `tests/formula_audit/test_theta_transition_forward_consistency.m`
  - shared theta-transition forward consistency audit.
- `tests/formula_audit/test_theta_transition_cache_selfcheck.m`
  - cache layout / reconstruction self-check.
- `tests/formula_audit/test_theta_transition_exact_fd_audit.m`
  - nodewise finite-difference audit of exact `dC/dphi`.
- `tests/formula_audit/test_theta_transition_directional_derivative_audit.m`
  - random-direction audit of exact `dC/dphi`.
- `tests/formula_audit/test_theta_transition_limiter_diagnostics.m`
  - limiter saturation / zero-gradient diagnostics.
- `tests/formula_audit/test_manufacturing_penalty_directional_derivatives.m`
  - directional derivative audit for manufacturing penalty terms.
- `tests/baseline/run_fast_baseline.m`
  - captures a structured regression baseline with metrics, plots, log, and printability outputs.
- `tests/comparison/run_refactor_strict_equivalence_check.m`
  - strict `isequaln` pre/post comparison against `refactor_artifacts/pre_default_*.mat`.
- `tests/comparison/run_p1p2_strict_equivalence_check.m`
  - strict `isequaln` check entry for P1/P2治理.
- `tests/comparison/run_hj_rhs_mode_benchmark.m`
  - HJ RHS mode performance benchmark and report export.
- `tests/test_postprocess_export_smoke.m`
  - known-geometry postprocess smoke test (segment count/length/point validity).
- `tests/test_exact_chain_default_smoke.m`
  - smoke test for the default exact-chain optimization path.

## What The Project Currently Optimizes For

From code and tests, the project is balancing four things at once:

- lower compliance
- stable angle updates
- signed-distance quality of `lsf`
- fiber continuity / printability

This explains why there are multiple guards, reinit checks, and a best-state rollback.

## Likely Performance Hotspots

The main runtime cost is likely dominated by repeated FE work and per-element loops:

- `FE_analysis_cantilever`
  - full assembly every evaluation
- `compute_sensitivity_adjoint`
  - per-element derivative assembly
- `compute_strain_energy`
  - another element loop
- repeated candidate evaluation inside one iteration
  - current state
  - theta-only candidate
  - HJ candidate
  - possible reinit candidate

So one optimizer iteration can contain multiple full FE solves.

## Likely Optimization Entry Points

If the goal is runtime reduction without changing behavior first, the best entry points are:

- reduce redundant FE evaluations across candidates
- cache/reuse `Ke` and `dKe_dtheta` patterns where possible
- vectorize or batch parts of element-level loops
- profile how often `theta_only`, `HJ`, and `reinit` all trigger expensive reevaluation
- cut diagnostic overhead in non-debug modes

If the goal is optimization quality instead of runtime, the best entry points are:

- acceptance thresholds
- reinitialization criteria
- narrow-band definition
- velocity scaling / clipping
- angle update limit `delta_theta_max`

## Known Practical Constraints

- MATLAB MCP expects Windows absolute paths in this environment.
- WSL-style `/mnt/e/...` paths work in shell commands here, but not in MATLAB MCP.
- `test_all_configurations.m` has already been verified through MATLAB MCP in this environment.
