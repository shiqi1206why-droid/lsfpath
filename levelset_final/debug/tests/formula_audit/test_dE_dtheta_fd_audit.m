% 论文核心公式审查包（只审查，不改公式）
% 对 compute_sensitivity_adjoint 给出的 dE/dtheta 进行有限差分抽样比对。

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

state = evaluate_state_with_theta(lsf, theta_eq, params.grid.nelx, params.grid.nely, material_mask, ...
    params.material.E_L, params.material.E_T, params.material.nu_LT, params.material.G_LT, ...
    params.material.thickness, params.load.F_mag, params.grid.dx, params.grid.dy);

sensitivity = compute_sensitivity_adjoint(params.grid.nelx, params.grid.nely, state.U, theta_eq, ...
    params.material.E_L, params.material.E_T, params.material.nu_LT, params.material.G_LT, ...
    params.material.thickness, params.grid.dx, params.grid.dy, material_mask, false);

valid_idx = find(material_mask);
[~, order] = sort(abs(sensitivity(valid_idx)), 'descend');
sample_count = min(6, numel(order));
sample_linear = valid_idx(order(1:sample_count));
[sample_y, sample_x] = ind2sub(size(material_mask), sample_linear);

eps_theta = deg2rad(0.1);
rows = zeros(sample_count, 9);

fprintf('=== dE/dtheta finite-difference audit ===\n');
fprintf('Base compliance: %.10e\n', state.compliance);
fprintf('eps_theta = %.6e rad\n', eps_theta);

for k = 1:sample_count
    iy = sample_y(k);
    ix = sample_x(k);

    theta_plus = theta_eq;
    theta_minus = theta_eq;
    theta_plus(iy, ix) = mod(theta_plus(iy, ix) + eps_theta, pi);
    theta_minus(iy, ix) = mod(theta_minus(iy, ix) - eps_theta, pi);

    state_plus = evaluate_state_with_theta(lsf, theta_plus, params.grid.nelx, params.grid.nely, material_mask, ...
        params.material.E_L, params.material.E_T, params.material.nu_LT, params.material.G_LT, ...
        params.material.thickness, params.load.F_mag, params.grid.dx, params.grid.dy);
    state_minus = evaluate_state_with_theta(lsf, theta_minus, params.grid.nelx, params.grid.nely, material_mask, ...
        params.material.E_L, params.material.E_T, params.material.nu_LT, params.material.G_LT, ...
        params.material.thickness, params.load.F_mag, params.grid.dx, params.grid.dy);

    fd = (state_plus.compliance - state_minus.compliance) / (2 * eps_theta);
    adj = sensitivity(iy, ix);
    sign_agree = (abs(fd) < 1e-12 && abs(adj) < 1e-12) || sign(fd) == sign(adj);
    ratio = NaN;
    rel_err = NaN;
    if abs(adj) > 1e-12
        ratio = fd / adj;
    end
    rel_err = abs(fd - adj) / max([abs(fd), abs(adj), 1e-12]);

    rows(k, :) = [k, iy, ix, adj, fd, ratio, rel_err, sign_agree, theta_eq(iy, ix)];
    fprintf(['sample=%d (y=%d,x=%d): adj=%.6e fd=%.6e ratio=%.6f ' ...
        'rel_err=%.6f sign_agree=%d theta=%.3f deg\n'], ...
        k, iy, ix, adj, fd, ratio, rel_err, sign_agree, theta_eq(iy, ix) * 180 / pi);
end

sign_agreement_rate = mean(rows(:, 8));
median_rel_err = median(rows(:, 7), 'omitnan');
median_ratio = median(rows(:, 6), 'omitnan');
formula_review_needed = sign_agreement_rate < 0.8 || median_rel_err > 0.5;

fprintf('sign_agreement_rate=%.3f\n', sign_agreement_rate);
fprintf('median_rel_err=%.6f\n', median_rel_err);
fprintf('median_ratio=%.6f\n', median_ratio);
fprintf('FORMULA_REVIEW_NEEDED=%d\n', formula_review_needed);

assert(sign_agreement_rate >= 0.95, 'dE/dtheta 与有限差分符号一致率过低。');
assert(median_rel_err <= 1e-2, 'dE/dtheta 与有限差分量级不一致。');

audit = struct();
audit.base_compliance = state.compliance;
audit.eps_theta = eps_theta;
audit.sign_agreement_rate = sign_agreement_rate;
audit.median_rel_err = median_rel_err;
audit.median_ratio = median_ratio;
audit.formula_review_needed = formula_review_needed;
audit.samples = rows;

disp(audit);
