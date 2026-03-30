clear; clc; close all;
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))), '-begin');

% 构造可控几何：圆形零等值线（含ghost cells）
nely = 40;
nelx = 60;
dx = 0.02;
dy = 0.02;
[x_full, y_full] = get_lsf_grid_coordinates([nely + 2, nelx + 2], dx, dy);
[X, Y] = meshgrid(x_full, y_full);
cx = mean(x_full);
cy = mean(y_full);
radius = 0.20;
lsf = hypot(X - cx, Y - cy) - radius;

material_mask_core = true(nely, nelx);
levels = 0;
out_dir = fullfile(tempdir, ['postprocess_export_smoke_' datestr(now, 'yyyymmdd_HHMMSSFFF')]);
mkdir(out_dir);

opts = struct();
opts.target_ds = min(dx, dy) / 3;
opts.smooth_window = 9;
opts.smooth_iters = 3;
opts.min_points = 20;
opts.auto_optimize = true;
opts.max_deviation_limit = min(dx, dy);
opts.candidate_windows = [7, 9];
opts.candidate_iters = [2, 3];
opts.candidate_methods = {'moving_average', 'chaikin'};
opts.chaikin_iters = [1, 2];

summary = export_printable_paths_from_lsf(lsf, dx, dy, levels, out_dir, material_mask_core, opts);

assert(summary.segment_count >= 1, '应至少导出一条路径段。');
assert(summary.sampled_point_count > 0, '应有采样点。');
assert(summary.sampled_point_violations == 0, '平滑路径应保持在材料域内。');
assert(isfile(summary.metrics_csv), 'metrics_csv应存在。');
assert(isfile(summary.overlay_png), 'overlay_png应存在。');
assert(isfile(summary.summary_txt), 'summary_txt应存在。');

T = readtable(summary.metrics_csv);
assert(height(T) >= 1, 'metrics表至少应包含一行。');
assert(all(T.smooth_n >= opts.min_points), '平滑路径点数应不小于最小点数阈值。');
assert(all(T.raw_len > 0) && all(T.smooth_len > 0), '路径长度应为正。');
assert(all(isfinite(T.smooth_max_abs_kappa)), '平滑曲率应为有限值。');

fprintf('segment_count=%d\n', summary.segment_count);
fprintf('sampled_point_count=%d\n', summary.sampled_point_count);
fprintf('PASS: postprocess export smoke test.\n');
