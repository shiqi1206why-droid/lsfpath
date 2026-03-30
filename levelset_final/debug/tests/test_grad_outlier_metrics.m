clear; clc; close all;
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))), '-begin');

nelx = 24;
nely = 18;
dx = 0.05;
dy = 0.05;
h = min(dx, dy);

material_mask_core = false(nely, nelx);
material_mask_core(3:16, 4:21) = true;
material_mask_full = expand_material_mask_to_full(material_mask_core);

[x_full, ~] = get_lsf_grid_coordinates([nely + 2, nelx + 2], dx, dy);
[X, ~] = meshgrid(x_full, 1:(nely + 2));
centered_x = X - mean(x_full(2:end-1));
lsf = 2.0 * centered_x;
lsf(~material_mask_full) = 0.8;
lsf(1, :) = lsf(2, :);
lsf(end, :) = lsf(end-1, :);
lsf(:, 1) = lsf(:, 2);
lsf(:, end) = lsf(:, end-1);

metrics = compute_raw_path_quality_metrics(lsf, dx, dy, material_mask_core, ...
    struct('parallel_spacing', h, 'resample_ds', h / 5, 'zero_bandwidth', h));

interior_material = false(size(lsf));
interior_material(2:end-1, 2:end-1) = true;
interior_material = interior_material & material_mask_full;
[grad_y, grad_x] = gradient(lsf, dy, dx);
grad_mag = hypot(grad_x, grad_y);
core_boundary_distance = bwdist(~material_mask_core) * h;
boundary_distance = inf(size(lsf));
boundary_distance(2:end-1, 2:end-1) = core_boundary_distance;
boundary_distance(1, :) = boundary_distance(2, :);
boundary_distance(end, :) = boundary_distance(end-1, :);
boundary_distance(:, 1) = boundary_distance(:, 2);
boundary_distance(:, end) = boundary_distance(:, end-1);
interior_high_mask = interior_material & grad_mag > 1.5;
expected_boundary_overlap = nnz(interior_high_mask & boundary_distance <= 2 * h) / nnz(interior_high_mask);

expected_outlier_ratio = metrics.gradient_scopes.interior_material_excluding_ghost.high_count / ...
    metrics.gradient_scopes.interior_material_excluding_ghost.valid_count;
assert(abs(metrics.grad_outlier_ratio_0p5_1p5 - expected_outlier_ratio) < 1e-12, ...
    '总 outlier 比例应与显式高梯度计数一致。');
assert(abs(metrics.near_zero_grad_outlier_ratio - 1) < 1e-12, '零线附近统一高梯度场应全部是 outlier。');
assert(abs(metrics.high_grad_outlier_ratio - expected_outlier_ratio) < 1e-12, ...
    'high ratio 应与显式高梯度计数一致。');
assert(abs(metrics.low_grad_outlier_ratio) < 1e-12, '统一高梯度场不应产生 low ratio。');
assert(abs(metrics.high_grad_boundary_overlap_ratio - expected_boundary_overlap) < 1e-12, ...
    '边界重叠比例应与显式计数一致。');
assert(abs(metrics.near_zero_grad_median - 2.0) < 1e-12, '零线附近梯度中位数应为 2。');

fprintf('near_zero_grad_median=%.6f\n', metrics.near_zero_grad_median);
fprintf('high_grad_boundary_overlap_ratio=%.6f\n', metrics.high_grad_boundary_overlap_ratio);
fprintf('PASS: gradient outlier metrics test.\n');
