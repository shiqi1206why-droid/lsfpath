clear; clc; close all;
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(script_dir));
addpath(fullfile(project_root, 'utilities'), '-begin');
project_root = get_project_root(project_root);
cleanup_path = ensure_project_on_path(project_root); %#ok<NASGU>
paths = build_project_paths(project_root);

params = get_fiber_optimization_params('fast');
params.runtime = struct('project_root', project_root, 'paths', paths);
topo = load(paths.topology_file);
[material_mask, ~] = clean_material_mask(topo.struc, params.init.min_component_size, params.init.morph_radius);

delta_phi = params.levelset.delta_phi_factor * params.grid.h;
[lsf, ~, ~] = construct_boundary_offset_levelset_with_parallel( ...
    material_mask, params.grid.nelx, params.grid.nely, ...
    params.grid.dx, params.grid.dy, delta_phi, struct('morph_radius', params.init.morph_radius));

% Case A: theta已经处于目标状态时，固定lsf不应继续漂移。
num_steps = 30;
[theta_eq, theta_target] = advance_theta_state(lsf, [], params.opt.delta_theta_max, ...
    params.grid.dx, params.grid.dy, params.smooth.eta, params.smooth.iterations);
theta_prev = theta_eq;
C_eq = zeros(num_steps, 1);

for k = 1:num_steps
    [theta_k, theta_target_k] = advance_theta_state(lsf, theta_prev, params.opt.delta_theta_max, ...
        params.grid.dx, params.grid.dy, params.smooth.eta, params.smooth.iterations);
    state_k = evaluate_state_with_theta(lsf, theta_k, params.grid.nelx, params.grid.nely, material_mask, ...
        params.material.E_L, params.material.E_T, params.material.nu_LT, params.material.G_LT, ...
        params.material.thickness, params.load.F_mag, params.grid.dx, params.grid.dy);
    C_eq(k) = state_k.compliance;
    assert(max(abs(theta_target_k(:) - theta_target(:))) < 1e-12, '固定lsf时theta_target不应变化。');
    theta_prev = theta_k;
end

drift_eq_percent = (C_eq(end) - C_eq(1)) / max(C_eq(1), eps) * 100;
max_step_eq_percent = max(abs(diff(C_eq)) ./ max(abs(C_eq(1:end-1)), eps)) * 100;

fprintf('Case A (equilibrium): C1=%.10e Cend=%.10e drift=%.6f%% max_step=%.6f%%\n', ...
    C_eq(1), C_eq(end), drift_eq_percent, max_step_eq_percent);

assert(abs(drift_eq_percent) < 0.1, '固定lsf且theta已达目标时，多步柔度漂移过大。');
assert(max_step_eq_percent < 0.1, '固定lsf且theta已达目标时，单步柔度波动过大。');

% Case B: theta尚未追上目标时，theta-only候选会改变柔度；若状态被冻结，则柔度必须保持不变。
theta_hold = mod(theta_eq + deg2rad(10), pi);
state_hold = evaluate_state_with_theta(lsf, theta_hold, params.grid.nelx, params.grid.nely, material_mask, ...
    params.material.E_L, params.material.E_T, params.material.nu_LT, params.material.G_LT, ...
    params.material.thickness, params.load.F_mag, params.grid.dx, params.grid.dy);
state_hold_replay = evaluate_state_with_theta(lsf, theta_hold, params.grid.nelx, params.grid.nely, material_mask, ...
    params.material.E_L, params.material.E_T, params.material.nu_LT, params.material.G_LT, ...
    params.material.thickness, params.load.F_mag, params.grid.dx, params.grid.dy);
theta_only_state = evaluate_candidate_state(lsf, theta_hold, params.opt.delta_theta_max, ...
    params.grid.dx, params.grid.dy, params.grid.nelx, params.grid.nely, material_mask, ...
    params.material.E_L, params.material.E_T, params.material.nu_LT, params.material.G_LT, ...
    params.material.thickness, params.load.F_mag, params.smooth.eta, params.smooth.iterations);

theta_step_deg = max(abs(atan2(sin(theta_only_state.theta(:) - theta_hold(:)), ...
    cos(theta_only_state.theta(:) - theta_hold(:))))) * 180 / pi;
theta_only_delta_percent = (theta_only_state.compliance - state_hold.compliance) / max(state_hold.compliance, eps) * 100;
hold_replay_delta_percent = (state_hold_replay.compliance - state_hold.compliance) / max(state_hold.compliance, eps) * 100;

fprintf('Case B (offset hold): C_hold=%.10e C_theta_only=%.10e delta=%.6f%% theta_step=%.6f deg hold_replay=%.6f%%\n', ...
    state_hold.compliance, theta_only_state.compliance, theta_only_delta_percent, theta_step_deg, hold_replay_delta_percent);

assert(theta_step_deg > 1e-6, 'theta-only候选未发生角度推进，测试场景无效。');
assert(abs(hold_replay_delta_percent) < 1e-10, '冻结接受状态后，固定lsf的柔度仍发生漂移。');

fprintf('PASS: StepA theta idempotence / hold-state test.\n');
