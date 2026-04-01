clear; clc; close all;
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(script_dir));
addpath(fullfile(project_root, 'utilities'), '-begin');
project_root = get_project_root(project_root);
cleanup_path = ensure_project_on_path(project_root); %#ok<NASGU>
paths = build_project_paths(project_root);

params = get_fiber_optimization_params('fast');
topo = load(paths.topology_file);
[material_mask, ~] = clean_material_mask(topo.struc, params.init.min_component_size, params.init.morph_radius);
material_mask_full = expand_material_mask_to_full(material_mask);

delta_phi = params.levelset.delta_phi_factor * params.grid.h;
[lsf, ~, init_info] = construct_boundary_offset_levelset_with_parallel( ...
    material_mask, params.grid.nelx, params.grid.nely, ...
    params.grid.dx, params.grid.dy, delta_phi, struct('morph_radius', params.init.morph_radius));

boundary_guard_band = abs(init_info.phi_boundary_full) <= params.levelset.boundary_guard_band_factor * params.grid.h;
boundary_guard_band = boundary_guard_band & material_mask_full;
primary_update_mask = (abs(lsf) <= 1.5 * params.grid.h) & material_mask_full & ~boundary_guard_band;
stencil_mask = dilate_binary_mask(primary_update_mask, params.levelset.stencil_buffer_cells) & material_mask_full;

current_state = evaluate_candidate_state(lsf, [], params.opt.delta_theta_max, ...
    params.grid.dx, params.grid.dy, params.grid.nelx, params.grid.nely, material_mask, ...
    params.material.E_L, params.material.E_T, params.material.nu_LT, params.material.G_LT, ...
    params.material.thickness, params.load.F_mag, params.smooth.eta, params.smooth.iterations);
theta_prev = current_state.theta;
theta_only_state = evaluate_candidate_state(lsf, theta_prev, params.opt.delta_theta_max, ...
    params.grid.dx, params.grid.dy, params.grid.nelx, params.grid.nely, material_mask, ...
    params.material.E_L, params.material.E_T, params.material.nu_LT, params.material.G_LT, ...
    params.material.thickness, params.load.F_mag, params.smooth.eta, params.smooth.iterations);

dC_dtheta_next = compute_sensitivity_adjoint(params.grid.nelx, params.grid.nely, ...
    theta_only_state.U, theta_only_state.theta, params.material.E_L, params.material.E_T, ...
    params.material.nu_LT, params.material.G_LT, params.material.thickness, ...
    params.grid.dx, params.grid.dy, material_mask, false);
[pullback, ~] = differentiate_theta_transition_exact(dC_dtheta_next, theta_only_state.theta_transition_cache, ...
    struct('limiter_mode', 'hard', ...
           'soft_limiter_beta', params.gradient.soft_limiter_beta, ...
           'grad_floor', params.gradient.theta_raw_grad_floor));
g_exact_opt = aggregate_node_sensitivity(pullback, [], lsf, params.grid.nelx, params.grid.nely, ...
    params.grid.dx, params.grid.dy, primary_update_mask);

node_sensitivity = g_exact_opt;
sens_scale_value = NaN;
sens_abs_band = abs(node_sensitivity(primary_update_mask));
if ~isempty(sens_abs_band) && params.velocity.scale_quantile > 0
    sens_scale_value = prctile(sens_abs_band, params.velocity.scale_quantile);
    if sens_scale_value > 0
        node_sensitivity = node_sensitivity / sens_scale_value;
    end
end
if isfinite(params.velocity.clip_abs)
    node_sensitivity = max(min(node_sensitivity, params.velocity.clip_abs), -params.velocity.clip_abs);
end

velocity_opts = params.velocity;
velocity_opts.exact_grad_floor = 1e-12;
[velocity_exact, velocity_stats] = build_velocity_exact( ...
    node_sensitivity, lsf, params.grid.dx, params.grid.dy, 1.5 * params.grid.h, ...
    velocity_opts.enable_bias_removal, 0, 1, primary_update_mask, velocity_opts);
band_velocity = velocity_exact;
band_velocity(~primary_update_mask) = 0;

dt_cfl = compute_adaptive_timestep(band_velocity, params.grid.dx, params.grid.dy);
if velocity_stats.max_band > 1e-12
    dt_angle = params.opt.delta_theta_max / velocity_stats.max_band;
else
    dt_angle = inf;
end
dt = min(dt_cfl, dt_angle);
hj_opts = struct( ...
    'advection_order', params.levelset.advection_order, ...
    'time_integrator', params.levelset.time_integrator, ...
    'fallback_first_order', logical(params.levelset.fallback_first_order), ...
    'freeze_on_incomplete_godunov', logical(params.levelset.freeze_on_incomplete_godunov), ...
    'stencil_mask', stencil_mask, ...
    'stencil_buffer_cells', params.levelset.stencil_buffer_cells, ...
    'eno_smoothness_factor', params.levelset.eno_smoothness_factor, ...
    'rhs_mode', char(params.levelset.hj_rhs_mode));
[phi_hj, ~] = update_levelset_HJ(lsf, velocity_exact, dt, params.grid.dx, params.grid.dy, primary_update_mask, hj_opts);

C_current = evaluate_candidate_compliance(lsf, theta_prev, params.opt.delta_theta_max, ...
    params.grid.dx, params.grid.dy, params.grid.nelx, params.grid.nely, material_mask, ...
    params.material.E_L, params.material.E_T, params.material.nu_LT, params.material.G_LT, ...
    params.material.thickness, params.load.F_mag, params.smooth.eta, params.smooth.iterations);
C_hj = evaluate_candidate_compliance(phi_hj, theta_prev, params.opt.delta_theta_max, ...
    params.grid.dx, params.grid.dy, params.grid.nelx, params.grid.nely, material_mask, ...
    params.material.E_L, params.material.E_T, params.material.nu_LT, params.material.G_LT, ...
    params.material.thickness, params.load.F_mag, params.smooth.eta, params.smooth.iterations);

fprintf('=== exact HJ descent smoke ===\n');
fprintf('C_current = %.6e\n', C_current);
fprintf('C_hj      = %.6e\n', C_hj);
fprintf('delta     = %.6e\n', C_hj - C_current);
fprintf('dt        = %.6e, dt_cfl = %.6e, dt_angle = %.6e, max_band = %.6e\n', ...
    dt, dt_cfl, dt_angle, velocity_stats.max_band);
fprintf('sens_scale = %.6e\n', sens_scale_value);

assert(C_hj < C_current, 'HJ exact 主循环口径下降验证失败。');
assert(C_hj > C_current * (1 - 0.10), 'HJ exact 单步下降幅度异常过大。');
fprintf('PASS: exact HJ descent smoke.\n');

function mask_full = expand_material_mask_to_full(mask_core)
mask_full = false(size(mask_core, 1) + 2, size(mask_core, 2) + 2);
mask_full(2:end-1, 2:end-1) = logical(mask_core);
mask_full(1, :) = mask_full(2, :);
mask_full(end, :) = mask_full(end-1, :);
mask_full(:, 1) = mask_full(:, 2);
mask_full(:, end) = mask_full(:, end-1);
end
