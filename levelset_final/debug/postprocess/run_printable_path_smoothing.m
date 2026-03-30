clear; clc; close all;

% 在debug副本中执行打印友好路径导出：
% 1) 跑一次fast得到最终lsf
% 2) 提取多条等值线路径
% 3) 平滑并导出CSV/图像/统计

script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
addpath(fullfile(project_root, 'utilities'), '-begin');
project_root = get_project_root(project_root);
cleanup_path = ensure_project_on_path(project_root); %#ok<NASGU>
paths = build_project_paths(project_root);

params = get_fiber_optimization_params('fast');
params.runtime = struct('project_root', project_root, 'paths', paths);
results = fiber_levelset('fast');

h = params.grid.h;
levels = (-2:2) * h;  % 5条主等值线: -2h,-h,0,h,2h

stamp = datestr(now, 'yyyymmdd_HHMMSS');
out_dir = fullfile(paths.baseline_dir, ['printable_paths_' stamp]);

opts = struct();
opts.target_ds = h / 4;
opts.smooth_window = 9;
opts.smooth_iters = 3;
opts.min_points = 20;
if isfield(results, 'path_quality_raw')
    opts.raw_path_quality = results.path_quality_raw;
end
if isfield(results, 'init_boundary_geometry')
    opts.init_boundary_geometry = results.init_boundary_geometry;
end

summary = export_printable_paths_from_lsf(results.lsf, params.grid.dx, params.grid.dy, ...
    levels, out_dir, results.material_mask_core, opts);

summary_json = fullfile(out_dir, 'printable_path_summary.json');
fid = fopen(summary_json, 'w');
if fid < 0
    error('无法写入: %s', summary_json);
end
cleanup_fid = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s', jsonencode(summary));

copyfile(summary.overlay_png, fullfile(paths.baseline_dir, 'printable_paths_latest_overlay.png'));
copyfile(summary.metrics_csv, fullfile(paths.baseline_dir, 'printable_paths_latest_metrics.csv'));
copyfile(summary.summary_txt, fullfile(paths.baseline_dir, 'printable_paths_latest_summary.txt'));
copyfile(summary_json, fullfile(paths.baseline_dir, 'printable_paths_latest_summary.json'));

fprintf('\n=== Printable Path Smoothing Complete ===\n');
fprintf('Output dir: %s\n', out_dir);
fprintf('Overlay:    %s\n', summary.overlay_png);
fprintf('Metrics:    %s\n', summary.metrics_csv);
fprintf('Summary:    %s\n', summary.summary_txt);
fprintf('Raw mean |dtheta| (deg): %.6f\n', summary.raw_mean_abs_turn_deg);
fprintf('Smooth mean |dtheta| (deg): %.6f\n', summary.smooth_mean_abs_turn_deg);
fprintf('Raw max kappa: %.6e\n', summary.raw_max_abs_kappa);
fprintf('Smooth max kappa: %.6e\n', summary.smooth_max_abs_kappa);
fprintf('Max deviation from raw (m): %.6e\n', summary.max_deviation_from_raw);
