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

current_state = evaluate_candidate_state(lsf, [], params.opt.delta_theta_max, ...
    params.grid.dx, params.grid.dy, params.grid.nelx, params.grid.nely, material_mask, ...
    params.material.E_L, params.material.E_T, params.material.nu_LT, params.material.G_LT, ...
    params.material.thickness, params.load.F_mag, params.smooth.eta, params.smooth.iterations);
theta_e = current_state.theta;
U = current_state.U;
F = current_state.F;
C_current = current_state.compliance;
theta_only_state = evaluate_candidate_state(lsf, theta_e, params.opt.delta_theta_max, ...
    params.grid.dx, params.grid.dy, params.grid.nelx, params.grid.nely, material_mask, ...
    params.material.E_L, params.material.E_T, params.material.nu_LT, params.material.G_LT, ...
    params.material.thickness, params.load.F_mag, params.smooth.eta, params.smooth.iterations);

band_mask = abs(lsf) <= 1.5 * params.grid.h;
material_mask_full = expand_material_mask_to_full(material_mask);
primary_update_mask = band_mask & material_mask_full;
gradient_in = struct();
gradient_in.current_state = current_state;
gradient_in.theta_only_state = theta_only_state;
gradient_in.lsf = lsf;
gradient_in.nelx = params.grid.nelx;
gradient_in.nely = params.grid.nely;
gradient_in.dx = params.grid.dx;
gradient_in.dy = params.grid.dy;
gradient_in.material_mask_core = material_mask;
gradient_in.material_mask_full = material_mask_full;
gradient_in.primary_update_mask = primary_update_mask;
gradient_in.gradient_opts = params.gradient;
gradient_in.normalize_sensitivity = params.opt.normalize_sensitivity;
gradient_in.E_L = params.material.E_L;
gradient_in.E_T = params.material.E_T;
gradient_in.nu_LT = params.material.nu_LT;
gradient_in.G_LT = params.material.G_LT;
gradient_in.thickness = params.material.thickness;
gradient_out = compute_gradient_chain_sensitivity(gradient_in);
node_sensitivity = gradient_out.chosen_node_sensitivity;
sens_abs_band = abs(node_sensitivity(primary_update_mask));
if ~isempty(sens_abs_band)
    p95_norm = prctile(sens_abs_band, params.velocity.scale_quantile);
    if p95_norm > 0
        node_sensitivity = node_sensitivity / p95_norm;
    end
end
if isfinite(params.velocity.clip_abs)
    node_sensitivity = max(min(node_sensitivity, params.velocity.clip_abs), -params.velocity.clip_abs);
end

[velocity, velocity_stats] = build_velocity_field(node_sensitivity, lsf, ...
    params.grid.dx, params.grid.dy, 1.5 * params.grid.h, ...
    params.velocity.enable_bias_removal, 0, 1, primary_update_mask, params.velocity);
velocity(~primary_update_mask) = 0;
assert(nnz(abs(velocity(~primary_update_mask)) > 1e-14) == 0, ...
    '传播域外速度必须为零。');

dt_cfl = compute_adaptive_timestep(velocity, params.grid.dx, params.grid.dy);
if velocity_stats.max_band > 1e-12
    dt_angle = params.opt.delta_theta_max / velocity_stats.max_band;
else
    dt_angle = inf;
end
dt = min(dt_cfl, dt_angle);

C_ref_next = evaluate_candidate_compliance(lsf, theta_e, params.opt.delta_theta_max, ...
    params.grid.dx, params.grid.dy, params.grid.nelx, params.grid.nely, material_mask, ...
    params.material.E_L, params.material.E_T, params.material.nu_LT, params.material.G_LT, ...
    params.material.thickness, params.load.F_mag, params.smooth.eta, params.smooth.iterations);
assert(isfinite(C_ref_next) && C_ref_next > 0, 'C_ref_next无效。');

step_accepted = false;
accepted_dt = NaN;
accepted_compliance = NaN;
last_trial = NaN;

for bt = 0:params.opt.max_backtrack
    dt_trial = dt * (params.opt.backtrack_factor ^ bt);
    if dt_trial < params.opt.min_backtrack_dt
        break;
    end

    lsf_trial = update_levelset_HJ(lsf, velocity, dt_trial, params.grid.dx, params.grid.dy, primary_update_mask);
    assert(max(abs(lsf_trial(~material_mask_full) - lsf(~material_mask_full))) < 1e-14, ...
        'HJ更新不应改变材料域外的水平集。');

    C_trial = evaluate_candidate_compliance(lsf_trial, theta_e, params.opt.delta_theta_max, ...
        params.grid.dx, params.grid.dy, params.grid.nelx, params.grid.nely, material_mask, ...
        params.material.E_L, params.material.E_T, params.material.nu_LT, params.material.G_LT, ...
        params.material.thickness, params.load.F_mag, params.smooth.eta, params.smooth.iterations);
    last_trial = C_trial;

    trial_ok_next = isfinite(C_trial) && ...
        C_trial <= C_ref_next * (1 + params.opt.acceptance_tol);
    trial_ok_current = isfinite(C_trial) && ...
        C_trial <= C_current * (1 + params.opt.current_state_tol);

    if trial_ok_next && trial_ok_current
        step_accepted = true;
        accepted_dt = dt_trial;
        accepted_compliance = C_trial;
        break;
    end
end

fprintf('C_current=%.10e\nC_ref_next=%.10e\n', C_current, C_ref_next);
if step_accepted
    fprintf('accepted_dt=%.6e C_trial=%.10e\n', accepted_dt, accepted_compliance);
    fprintf('trial_vs_ref=%.6f%% trial_vs_current=%.6f%%\n', ...
        (accepted_compliance - C_ref_next) / max(C_ref_next, eps) * 100, ...
        (accepted_compliance - C_current) / max(C_current, eps) * 100);
    assert(accepted_compliance <= C_ref_next * (1 + params.opt.acceptance_tol + 1e-12), ...
        '接受的C_trial未满足相对C_ref_next的接受条件。');
    assert(accepted_compliance <= C_current * (1 + params.opt.current_state_tol + 1e-12), ...
        '接受的C_trial未满足相对C_current的接受条件。');
else
    fprintf('no acceptable trial found, last_C_trial=%.10e\n', last_trial);
    assert(~isfinite(last_trial) || ...
        last_trial > C_ref_next * (1 + params.opt.acceptance_tol) || ...
        last_trial > C_current * (1 + params.opt.current_state_tol), ...
        '未接受步却没有触发任何接受门禁。');
end

fprintf('PASS: StepA acceptance self-consistency test.\n');

function mask_full = expand_material_mask_to_full(mask_core)
mask_full = false(size(mask_core, 1) + 2, size(mask_core, 2) + 2);
mask_full(2:end-1, 2:end-1) = logical(mask_core);
mask_full(1, :) = mask_full(2, :);
mask_full(end, :) = mask_full(end-1, :);
mask_full(:, 1) = mask_full(:, 2);
mask_full(:, end) = mask_full(:, end-1);
end
