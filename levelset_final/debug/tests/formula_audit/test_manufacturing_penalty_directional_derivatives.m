clear; clc; close all;
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(script_dir));
addpath(fullfile(project_root, 'utilities'), '-begin');
project_root = get_project_root(project_root);
cleanup_path = ensure_project_on_path(project_root); %#ok<NASGU>

nelx = 20;
nely = 14;
dx = 0.04;
dy = 0.05;
[X, Y] = meshgrid(0:nelx+1, 0:nely+1);
lsf = 0.08 * sin(0.27 * X) + 0.06 * cos(0.19 * Y) + 0.02 * X;
lsf_target_global = lsf + 0.01 * sin(0.13 * X);
material_mask_full = true(size(lsf));
active_mask = abs(lsf) <= 1.5 * min(dx, dy);
eps_phi = 1e-4 * min(dx, dy);

rng(5);
v = randn(size(lsf));
v(~active_mask) = 0;
v = v / max(norm(v(:)), 1e-12);

cases = { ...
    struct('name', 'grad_norm', 'grad_norm_weight', 1e-2, 'curvature_weight', 0, 'gap_overlap_weight', 0), ...
    struct('name', 'curvature', 'grad_norm_weight', 0, 'curvature_weight', 1e-3, 'gap_overlap_weight', 0), ...
    struct('name', 'gap_overlap', 'grad_norm_weight', 0, 'curvature_weight', 0, 'gap_overlap_weight', 1e-2)};

fprintf('=== manufacturing penalty directional derivative audit ===\n');
for k = 1:numel(cases)
    opts = struct();
    opts.enable = true;
    opts.grad_norm_weight = cases{k}.grad_norm_weight;
    opts.curvature_weight = cases{k}.curvature_weight;
    opts.gap_overlap_weight = cases{k}.gap_overlap_weight;
    opts.curvature_radius_min = 2.0;
    opts.gap_overlap_target = 1.0;
    opts.penalty_band_factor = 1.5;

    [g0, diag0] = compute_manufacturing_penalty_gradient(lsf, dx, dy, material_mask_full, lsf_target_global, active_mask, opts);
    [~, diag_plus] = compute_manufacturing_penalty_gradient(lsf + eps_phi * v, dx, dy, material_mask_full, lsf_target_global, active_mask, opts);
    [~, diag_minus] = compute_manufacturing_penalty_gradient(lsf - eps_phi * v, dx, dy, material_mask_full, lsf_target_global, active_mask, opts);

    E0 = diag0.grad_norm_energy + diag0.curvature_energy + diag0.gap_overlap_energy;
    E_plus = diag_plus.grad_norm_energy + diag_plus.curvature_energy + diag_plus.gap_overlap_energy;
    E_minus = diag_minus.grad_norm_energy + diag_minus.curvature_energy + diag_minus.gap_overlap_energy;
    fd_dir = (E_plus - E_minus) / (2 * eps_phi);
    analytic_dir = sum(g0(:) .* v(:));
    rel_err = abs(fd_dir - analytic_dir) / max([abs(fd_dir), abs(analytic_dir), 1e-12]);

    fprintf('%s: E0=%.6e fd=%.6e analytic=%.6e rel_err=%.4f\n', ...
        cases{k}.name, E0, fd_dir, analytic_dir, rel_err);
    assert(rel_err <= 2e-1, '制造约束 %s 方向导数审计未通过。', cases{k}.name);
end

fprintf('PASS: manufacturing penalty directional derivative audit.\n');

