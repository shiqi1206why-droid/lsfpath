% HJ RHS模式性能基准（legacy/indexed/vectorized_first_order）
% 输出报告到 refactor_artifacts

clc; clear; close all;
fprintf('=== HJ RHS Mode Benchmark ===\n');

script_path = mfilename('fullpath');
if isempty(script_path)
    script_path = which('run_hj_rhs_mode_benchmark.m');
end
if isempty(script_path)
    error('无法定位 run_hj_rhs_mode_benchmark.m');
end
script_dir = fileparts(script_path);
project_root = fileparts(fileparts(script_dir));

addpath(fullfile(project_root, 'utilities'), '-begin');
project_root = get_project_root(project_root);
cleanup_path = ensure_project_on_path(project_root); %#ok<NASGU>
paths = build_project_paths(project_root);

params = get_fiber_optimization_params('default');
ny = params.grid.nely + 2;
nx = params.grid.nelx + 2;
dx = params.grid.dx;
dy = params.grid.dy;

[x_full, y_full] = get_lsf_grid_coordinates([ny, nx], dx, dy);
[X, Y] = meshgrid(x_full, y_full);
lsf = 0.4 * sin(2*pi*X/max(x_full)) + 0.3 * cos(2*pi*Y/max(y_full)) + 0.1 * X;
velocity = 0.5 * sin(2*pi*Y/max(y_full)) + 0.2 * cos(2*pi*X/max(x_full));
dt = 0.02;

active_mask = abs(lsf) <= 2.0 * min(dx, dy);
active_mask(:, [1, end]) = false;
active_mask([1, end], :) = false;

base_opts = struct( ...
    'advection_order', 2, ...
    'time_integrator', 'ssprk2', ...
    'fallback_first_order', true, ...
    'freeze_on_incomplete_godunov', true, ...
    'stencil_buffer_cells', 2, ...
    'eno_smoothness_factor', 2.5, ...
    'stencil_mask', active_mask);

cases = {
    struct('name', 'legacy_order2', 'rhs_mode', 'legacy', 'advection_order', 2), ...
    struct('name', 'indexed_order2', 'rhs_mode', 'indexed', 'advection_order', 2), ...
    struct('name', 'legacy_order1', 'rhs_mode', 'legacy', 'advection_order', 1), ...
    struct('name', 'vectorized_first_order1', 'rhs_mode', 'vectorized_first_order', 'advection_order', 1)
};

n_repeat = 8;
results = struct('name', {}, 'rhs_mode', {}, 'advection_order', {}, 'elapsed_sec', {}, ...
    'elapsed_mean_sec', {}, 'diag', {});

for ci = 1:numel(cases)
    c = cases{ci};
    opts = base_opts;
    opts.rhs_mode = c.rhs_mode;
    opts.advection_order = c.advection_order;

    % warm-up
    update_levelset_HJ(lsf, velocity, dt, dx, dy, active_mask, opts);

    elapsed = nan(n_repeat, 1);
    diag_last = struct();
    for k = 1:n_repeat
        t0 = tic;
        [~, diag_last] = update_levelset_HJ(lsf, velocity, dt, dx, dy, active_mask, opts);
        elapsed(k) = toc(t0);
    end

    entry = struct();
    entry.name = c.name;
    entry.rhs_mode = c.rhs_mode;
    entry.advection_order = c.advection_order;
    entry.elapsed_sec = elapsed;
    entry.elapsed_mean_sec = mean(elapsed);
    entry.diag = diag_last;
    results(end+1) = entry; %#ok<AGROW>
end

opts_legacy2 = base_opts; opts_legacy2.rhs_mode = 'legacy'; opts_legacy2.advection_order = 2;
opts_index2 = base_opts; opts_index2.rhs_mode = 'indexed'; opts_index2.advection_order = 2;
[phi_legacy2, ~] = update_levelset_HJ(lsf, velocity, dt, dx, dy, active_mask, opts_legacy2);
[phi_index2, ~] = update_levelset_HJ(lsf, velocity, dt, dx, dy, active_mask, opts_index2);
eq_legacy_vs_indexed_order2 = isequaln(phi_legacy2, phi_index2);

opts_legacy1 = base_opts; opts_legacy1.rhs_mode = 'legacy'; opts_legacy1.advection_order = 1;
opts_vec1 = base_opts; opts_vec1.rhs_mode = 'vectorized_first_order'; opts_vec1.advection_order = 1;
[phi_legacy1, ~] = update_levelset_HJ(lsf, velocity, dt, dx, dy, active_mask, opts_legacy1);
[phi_vec1, ~] = update_levelset_HJ(lsf, velocity, dt, dx, dy, active_mask, opts_vec1);
eq_legacy_vs_vectorized_order1 = isequaln(phi_legacy1, phi_vec1);

legacy2_time = get_time(results, 'legacy_order2');
indexed2_time = get_time(results, 'indexed_order2');
legacy1_time = get_time(results, 'legacy_order1');
vec1_time = get_time(results, 'vectorized_first_order1');

report = struct();
report.timestamp = datestr(now, 'yyyymmdd_HHMMSS');
report.n_repeat = n_repeat;
report.results = results;
report.eq_legacy_vs_indexed_order2 = eq_legacy_vs_indexed_order2;
report.eq_legacy_vs_vectorized_order1 = eq_legacy_vs_vectorized_order1;
report.speedup_indexed_vs_legacy_order2 = legacy2_time / max(indexed2_time, eps);
report.speedup_vectorized_vs_legacy_order1 = legacy1_time / max(vec1_time, eps);

stamp = report.timestamp;
report_mat = fullfile(paths.refactor_dir, sprintf('hj_rhs_benchmark_%s.mat', stamp));
report_txt = fullfile(paths.refactor_dir, sprintf('hj_rhs_benchmark_%s.txt', stamp));
save(report_mat, 'report');

fid = fopen(report_txt, 'w');
if fid == -1
    warning('无法写入基准报告: %s', report_txt);
else
    fprintf(fid, 'HJ RHS benchmark report\n');
    fprintf(fid, 'timestamp: %s\n', stamp);
    fprintf(fid, 'n_repeat: %d\n\n', n_repeat);
    for i = 1:numel(results)
        r = results(i);
        fprintf(fid, '%s: rhs_mode=%s, advection_order=%d, mean_time=%.6f sec\n', ...
            r.name, r.rhs_mode, r.advection_order, r.elapsed_mean_sec);
    end
    fprintf(fid, '\neq_legacy_vs_indexed_order2 = %d\n', eq_legacy_vs_indexed_order2);
    fprintf(fid, 'eq_legacy_vs_vectorized_order1 = %d\n', eq_legacy_vs_vectorized_order1);
    fprintf(fid, 'speedup_indexed_vs_legacy_order2 = %.6f\n', report.speedup_indexed_vs_legacy_order2);
    fprintf(fid, 'speedup_vectorized_vs_legacy_order1 = %.6f\n', report.speedup_vectorized_vs_legacy_order1);
    fclose(fid);
end

fprintf('legacy(order2) mean: %.6f sec\n', legacy2_time);
fprintf('indexed(order2) mean: %.6f sec\n', indexed2_time);
fprintf('vectorized_first_order(order1) mean: %.6f sec\n', vec1_time);
fprintf('eq legacy/indexed (order2): %d\n', eq_legacy_vs_indexed_order2);
fprintf('eq legacy/vectorized_first_order (order1): %d\n', eq_legacy_vs_vectorized_order1);
fprintf('speedup indexed vs legacy (order2): %.3fx\n', report.speedup_indexed_vs_legacy_order2);
fprintf('speedup vectorized vs legacy (order1): %.3fx\n', report.speedup_vectorized_vs_legacy_order1);
fprintf('Saved MAT report: %s\n', report_mat);
fprintf('Saved TXT report: %s\n', report_txt);

function t = get_time(results, name)
idx = find(strcmp({results.name}, name), 1, 'first');
if isempty(idx)
    t = NaN;
else
    t = results(idx).elapsed_mean_sec;
end
end
