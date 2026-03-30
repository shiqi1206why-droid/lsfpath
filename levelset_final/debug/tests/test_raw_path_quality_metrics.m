clear; clc; close all;
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))), '-begin');

nelx = 26;
nely = 18;
dx = 0.05;
dy = 0.05;

material_mask_core = false(nely, nelx);
material_mask_core(3:16, 4:23) = true;
material_mask_full = expand_material_mask_to_full(material_mask_core);

[x_full, y_full] = get_lsf_grid_coordinates([nely + 2, nelx + 2], dx, dy);
[X, Y] = meshgrid(x_full, y_full);
cx = mean(x_full(5:end-4));
cy = mean(y_full(5:end-4));
r = 0.20;
lsf = sqrt((X - cx).^2 + (Y - cy).^2) - r;
lsf(~material_mask_full) = 0.4;
lsf = impose_neumann(lsf);

metrics = compute_raw_path_quality_metrics(lsf, dx, dy, material_mask_core, ...
    struct('parallel_spacing', min(dx, dy), 'resample_ds', min(dx, dy) / 5));

assert(metrics.segment_count > 0, 'raw 路径质量评估未提取到零等值线。');
assert(isfinite(metrics.mean_abs_turn_deg), 'mean_abs_turn_deg 应为有限值。');
assert(isfinite(metrics.max_abs_kappa), 'max_abs_kappa 应为有限值。');
assert(isfinite(metrics.parallel_spacing_error_percent), 'parallel_spacing_error_percent 应为有限值。');
assert(isfinite(metrics.grad_dev_mean), 'grad_dev_mean 应为有限值。');
assert(isfinite(metrics.grad_outlier_ratio_0p5_1p5), 'grad_outlier_ratio_0p5_1p5 应为有限值。');
assert(isfinite(metrics.near_zero_grad_outlier_ratio), 'near_zero_grad_outlier_ratio 应为有限值。');
assert(isfinite(metrics.high_grad_boundary_overlap_ratio), 'high_grad_boundary_overlap_ratio 应为有限值。');
assert(metrics.outside_nonpositive_count == 0, '材料域外不应出现 phi<=0。');

fprintf('mean_abs_turn_deg=%.6f\n', metrics.mean_abs_turn_deg);
fprintf('max_abs_kappa=%.6e\n', metrics.max_abs_kappa);
fprintf('parallel_spacing_error_percent=%.6f\n', metrics.parallel_spacing_error_percent);
fprintf('grad_dev_mean=%.6e\n', metrics.grad_dev_mean);
fprintf('near_zero_grad_outlier_ratio=%.6f\n', metrics.near_zero_grad_outlier_ratio);
fprintf('PASS: raw path quality metrics test.\n');

function field = impose_neumann(field)
field(1, :) = field(2, :);
field(end, :) = field(end-1, :);
field(:, 1) = field(:, 2);
field(:, end) = field(:, end-1);
end
