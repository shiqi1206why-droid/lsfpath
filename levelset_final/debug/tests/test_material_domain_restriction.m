clear; clc; close all;
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
addpath(fullfile(project_root, 'utilities'), '-begin');
project_root = get_project_root(project_root);
cleanup_path = ensure_project_on_path(project_root); %#ok<NASGU>

nelx = 20;
nely = 14;
dx = 0.05;
dy = 0.05;
h = min(dx, dy);

material_mask_core = false(nely, nelx);
material_mask_core(3:12, 4:17) = true;
material_mask_full = expand_material_mask_to_full(material_mask_core);

[lsf, ~, ~] = construct_boundary_offset_levelset_with_parallel( ...
    material_mask_core, nelx, nely, dx, dy, 0.8 * h, struct());

propagation_mask = (abs(lsf) <= 1.5 * h) & material_mask_full;
node_sensitivity = randn(size(lsf));
velocity_opts = struct('bias_beta', 0.10);

[velocity, velocity_stats] = build_velocity_field(node_sensitivity, lsf, dx, dy, ...
    1.5 * h, true, 0, 1, propagation_mask, velocity_opts);

outside_velocity_nonzero = nnz(abs(velocity(~propagation_mask)) > 1e-14);
fprintf('outside_velocity_nonzero=%d\n', outside_velocity_nonzero);
assert(outside_velocity_nonzero == 0, '传播域外速度必须为零。');

if velocity_stats.max_band > 1e-12
    dt = min(0.1, 0.25 * h / velocity_stats.max_band);
else
    dt = 0.05;
end

lsf_trial = update_levelset_HJ(lsf, velocity, dt, dx, dy, propagation_mask);
outside_change = max(abs(lsf_trial(~material_mask_full) - lsf(~material_mask_full)));
fprintf('outside_change=%e\n', outside_change);
assert(outside_change < 1e-14, 'HJ更新不应改变材料域外水平集。');

zero_mask_dynamic = compute_zero_mask_from_lsf(lsf_trial, h);
zero_mask_dynamic(1, 1) = true;
[lsf_reinit, reinit_diag] = fmm_reinitialize(lsf_trial, dx, dy, zero_mask_dynamic, ...
    material_mask_core, struct('method', 'subcell_signed_distance'));

zero_outside = nnz(abs(lsf_reinit) <= 1e-12 & ~material_mask_full);
nonpositive_outside = nnz(lsf_reinit(~material_mask_full) <= 0);
fprintf('zero_outside=%d\n', zero_outside);
fprintf('nonpositive_outside=%d\n', nonpositive_outside);
fprintf('used_method=%s\n', reinit_diag.used_method);

assert(zero_outside == 0, '重初始化后不应在void域产生零等值线。');
assert(nonpositive_outside == 0, '重初始化后void域必须保持严格正值。');
assert(~isempty(reinit_diag.used_method), '重初始化诊断应返回方法信息。');

fprintf('PASS: material-domain restriction test.\n');

function mask_full = expand_material_mask_to_full(mask_core)
mask_full = false(size(mask_core, 1) + 2, size(mask_core, 2) + 2);
mask_full(2:end-1, 2:end-1) = logical(mask_core);
mask_full(1, :) = mask_full(2, :);
mask_full(end, :) = mask_full(end-1, :);
mask_full(:, 1) = mask_full(:, 2);
mask_full(:, end) = mask_full(:, end-1);
end
