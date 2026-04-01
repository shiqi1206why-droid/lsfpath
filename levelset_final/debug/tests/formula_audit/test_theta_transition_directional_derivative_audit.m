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
material_mask_full = expand_material_mask_to_full(material_mask);

delta_phi = params.levelset.delta_phi_factor * params.grid.h;
[lsf, ~, ~] = construct_boundary_offset_levelset_with_parallel( ...
    material_mask, params.grid.nelx, params.grid.nely, ...
    params.grid.dx, params.grid.dy, delta_phi, struct('morph_radius', params.init.morph_radius));

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

g_exact = aggregate_node_sensitivity(pullback, [], lsf, params.grid.nelx, params.grid.nely, ...
    params.grid.dx, params.grid.dy, true(size(lsf)));

[audit_mask, audit_diag] = build_theta_transition_audit_support( ...
    theta_only_state.theta_transition_cache, ...
    struct('material_mask_full', material_mask_full, ...
           'support_mode', params.gradient.audit_support_mode, ...
           'grad_floor', params.gradient.theta_raw_grad_floor, ...
           'branch_cut_tol', params.gradient.theta_raw_branch_cut_tol));
assert(any(audit_mask(:)), 'exact directional derivative audit support为空。');

candidate_idx = find(audit_mask & abs(g_exact) > 1e-12);
if isempty(candidate_idx)
    candidate_idx = find(audit_mask);
end
[~, order] = sort(abs(g_exact(candidate_idx)), 'descend');
keep_count = min(max(32, min(64, numel(order))), numel(order));
direction_mask = false(size(audit_mask));
direction_mask(candidate_idx(order(1:keep_count))) = true;

eps_phi = params.gradient.shadow_fd_eps_factor * params.grid.h;
num_checks = params.gradient.shadow_num_directional_checks;
rel_err = zeros(num_checks, 1);

rng(23);
fprintf('=== exact directional derivative audit ===\n');
fprintf('eps_phi = %.6e, checks = %d\n', eps_phi, num_checks);
fprintf('audit_support_fraction = %.6f, branch_guard_ratio = %.6f\n', ...
    audit_diag.audit_support_fraction, audit_diag.branch_guard_ratio);
fprintf('direction_support_fraction = %.6f\n', nnz(direction_mask) / numel(direction_mask));

for k = 1:num_checks
    v = randn(size(lsf));
    v(~direction_mask) = 0;
    v = v / max(norm(v(:)), 1e-12);

    C_plus = evaluate_candidate_compliance(lsf + eps_phi * v, theta_prev, params.opt.delta_theta_max, ...
        params.grid.dx, params.grid.dy, params.grid.nelx, params.grid.nely, material_mask, ...
        params.material.E_L, params.material.E_T, params.material.nu_LT, params.material.G_LT, ...
        params.material.thickness, params.load.F_mag, params.smooth.eta, params.smooth.iterations);
    C_minus = evaluate_candidate_compliance(lsf - eps_phi * v, theta_prev, params.opt.delta_theta_max, ...
        params.grid.dx, params.grid.dy, params.grid.nelx, params.grid.nely, material_mask, ...
        params.material.E_L, params.material.E_T, params.material.nu_LT, params.material.G_LT, ...
        params.material.thickness, params.load.F_mag, params.smooth.eta, params.smooth.iterations);

    fd_dir = (C_plus - C_minus) / (2 * eps_phi);
    analytic_dir = sum(g_exact(:) .* v(:));
    rel_err(k) = abs(fd_dir - analytic_dir) / max([abs(fd_dir), abs(analytic_dir), 1e-12]);

    fprintf('check=%d: fd=%.6e analytic=%.6e rel_err=%.4f\n', ...
        k, fd_dir, analytic_dir, rel_err(k));
end

median_rel_err = median(rel_err, 'omitnan');
fprintf('median_rel_err=%.6f\n', median_rel_err);

assert(median_rel_err <= 5e-2, 'exact directional derivative 审计未通过。');
fprintf('PASS: exact directional derivative audit.\n');

function mask_full = expand_material_mask_to_full(mask_core)
mask_full = false(size(mask_core, 1) + 2, size(mask_core, 2) + 2);
mask_full(2:end-1, 2:end-1) = logical(mask_core);
mask_full(1, :) = mask_full(2, :);
mask_full(end, :) = mask_full(end-1, :);
mask_full(:, 1) = mask_full(:, 2);
mask_full(:, end) = mask_full(:, end-1);
end
