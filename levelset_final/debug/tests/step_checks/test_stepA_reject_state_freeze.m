% StepA新增测试：当HJ/reinit候选被拒绝时，接受状态必须冻结，不能因theta继续推进而漂移。

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

[theta_eq, ~] = advance_theta_state(lsf, [], params.opt.delta_theta_max, ...
    params.grid.dx, params.grid.dy, params.smooth.eta, params.smooth.iterations);

theta_hold = mod(theta_eq + deg2rad(12), pi);
current_state = evaluate_state_with_theta(lsf, theta_hold, params.grid.nelx, params.grid.nely, material_mask, ...
    params.material.E_L, params.material.E_T, params.material.nu_LT, params.material.G_LT, ...
    params.material.thickness, params.load.F_mag, params.grid.dx, params.grid.dy);

theta_only_state = evaluate_candidate_state(lsf, theta_hold, params.opt.delta_theta_max, ...
    params.grid.dx, params.grid.dy, params.grid.nelx, params.grid.nely, material_mask, ...
    params.material.E_L, params.material.E_T, params.material.nu_LT, params.material.G_LT, ...
    params.material.thickness, params.load.F_mag, params.smooth.eta, params.smooth.iterations);

% 模拟 “HJ/reinit被拒绝，accepted_source='hold'” 的下一个接受状态。
accepted_state = current_state;
replayed_state = evaluate_state_with_theta(accepted_state.lsf, accepted_state.theta, ...
    params.grid.nelx, params.grid.nely, material_mask, params.material.E_L, params.material.E_T, ...
    params.material.nu_LT, params.material.G_LT, params.material.thickness, ...
    params.load.F_mag, params.grid.dx, params.grid.dy);

theta_only_delta_deg = max(abs(atan2(sin(theta_only_state.theta(:) - theta_hold(:)), ...
    cos(theta_only_state.theta(:) - theta_hold(:))))) * 180 / pi;
theta_only_delta_percent = (theta_only_state.compliance - current_state.compliance) / max(current_state.compliance, eps) * 100;
freeze_delta_percent = (replayed_state.compliance - current_state.compliance) / max(current_state.compliance, eps) * 100;

fprintf('C_current=%.10e\n', current_state.compliance);
fprintf('C_theta_only=%.10e delta=%.6f%% theta_step=%.6f deg\n', ...
    theta_only_state.compliance, theta_only_delta_percent, theta_only_delta_deg);
fprintf('C_replayed=%.10e freeze_delta=%.6f%%\n', ...
    replayed_state.compliance, freeze_delta_percent);

assert(theta_only_delta_deg > 1e-6, 'theta-only候选未推进，测试场景无效。');
assert(abs(freeze_delta_percent) < 1e-10, '拒绝步冻结后，柔度仍发生漂移。');

fprintf('PASS: StepA reject-state freeze test.\n');
