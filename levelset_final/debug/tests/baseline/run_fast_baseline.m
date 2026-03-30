% Baseline capture for debug copy only.
% This script runs fiber_levelset('fast') once, then stores:
% 1) machine-readable metrics JSON
% 2) markdown report
% 3) key plots and full command-window log

clc;
close all;

script_dir = fileparts(mfilename('fullpath'));
project_dir = fileparts(fileparts(script_dir));  % <project>/tests/baseline
addpath(fullfile(project_dir, 'utilities'), '-begin');
project_root = get_project_root(project_dir);
cleanup_path = ensure_project_on_path(project_root); %#ok<NASGU>
paths = build_project_paths(project_root);
set(0, 'DefaultFigureVisible', 'off');

output_root = paths.baseline_dir;
if ~exist(output_root, 'dir')
    mkdir(output_root);
end

run_stamp = datestr(now, 'yyyymmdd_HHMMSS');
run_dir = fullfile(output_root, ['baseline_' run_stamp]);
if ~exist(run_dir, 'dir')
    mkdir(run_dir);
end

if exist(paths.topology_file, 'file') ~= 2
    error('topo_result.mat not found in debug project root.');
end

function_paths = struct( ...
    'fiber_levelset', fullfile(project_root, 'fiber_levelset.m'), ...
    'verify_path_spacing', fullfile(project_root, 'initialization', 'verify_path_spacing.m'), ...
    'aggregate_node_sensitivity', fullfile(project_root, 'level_set_evolution', 'aggregate_node_sensitivity.m'));

previous_latest_metrics = load_previous_latest_metrics(output_root);

run_log = '';
runtime_seconds = NaN;
results = [];
max_attempts = 2;

for attempt = 1:max_attempts
    run_timer = tic;
    try
        run_log = evalc('results = fiber_levelset(''fast'');');
        runtime_seconds = toc(run_timer);
        break;
    catch ME
        runtime_seconds = toc(run_timer);
        err_text = compose_exception_text(ME);
        write_text_file(fullfile(output_root, 'baseline_probe_error.txt'), err_text);
        write_text_file(fullfile(run_dir, sprintf('run_error_attempt%d.txt', attempt)), err_text);
        if attempt < max_attempts && is_memory_error(ME)
            warning('run_fast_baseline:retryOnMemoryError', ...
                'Attempt %d/%d failed due to memory pressure. Retrying once after cleanup.', ...
                attempt, max_attempts);
            try
                close all force;
            catch
            end
            try
                clear mex;
            catch
            end
            try
                java.lang.System.gc;
            catch
            end
            pause(1.0);
            continue;
        end
        rethrow(ME);
    end
end

params = results.params;

lsf_change_series = extract_numeric_series(run_log, '水平集变化:\s*([-\d\.eE+]+)');
dev95_series = extract_numeric_series(run_log, 'dev95=([-\d\.eE+]+)');
mean_off_series = extract_numeric_series(run_log, 'mean_off=([-\d\.eE+]+)');
target_delta_phi = extract_first_numeric(run_log, '目标Δφ =\s*([-\d\.eE+]+)\s*m');
mean_offset = extract_first_numeric(run_log, '抽样均值 =\s*([-\d\.eE+]+)\s*m');
std_offset = extract_first_numeric(run_log, '标准差 =\s*([-\d\.eE+]+)\s*m');
num_samples = extract_first_numeric(run_log, '样本数 =\s*(\d+)');
max_offset_error = extract_first_numeric(run_log, '最大误差 =\s*([-\d\.eE+]+)\s*m');
thin_ratio_percent = extract_first_numeric(run_log, '薄壁警告：([-\d\.eE+]+)%');
spacing_error_percent = extract_first_numeric(run_log, '平行路径平均间距误差:\s*([-\d\.eE+]+)%');

if isnan(thin_ratio_percent)
    thin_ratio = 0;
else
    thin_ratio = thin_ratio_percent / 100;
end

target_spacing = params.grid.h;
if isnan(spacing_error_percent)
    measured_spacing = NaN;
else
    measured_spacing = target_spacing * (1 + spacing_error_percent / 100);
end

if numel(results.compliance_history) > 1
    compliance_diff = diff(results.compliance_history);
    compliance_up_steps = nnz(compliance_diff > 0);
    compliance_total_steps = numel(compliance_diff);
else
    compliance_up_steps = 0;
    compliance_total_steps = 0;
end

if isfield(results, 'best_compliance') && isfield(results, 'best_iter')
    best_compliance = results.best_compliance;
    best_iter = results.best_iter;
else
    [best_compliance, best_iter] = min(results.compliance_history);
end
rollback_to_best = isfield(results, 'rollback_to_best') && results.rollback_to_best;
early_stop_triggered = isfield(results, 'early_stop_triggered') && results.early_stop_triggered;
if isfield(results, 'reject_due_next_guard')
    reject_due_next_guard = results.reject_due_next_guard;
else
    reject_due_next_guard = NaN;
end
if isfield(results, 'reject_due_current_guard')
    reject_due_current_guard = results.reject_due_current_guard;
else
    reject_due_current_guard = NaN;
end
if isfield(results, 'reinit_skip_due_next_guard')
    reinit_skip_due_next_guard = results.reinit_skip_due_next_guard;
else
    reinit_skip_due_next_guard = NaN;
end
if isfield(results, 'reinit_skip_due_current_guard')
    reinit_skip_due_current_guard = results.reinit_skip_due_current_guard;
else
    reinit_skip_due_current_guard = NaN;
end
if isfield(results, 'executed_iter')
    executed_iter = results.executed_iter;
else
    executed_iter = results.final_iter;
end
if isfield(results, 'loop_iter_count')
    loop_iter_count = results.loop_iter_count;
else
    loop_iter_count = max(0, numel(results.compliance_history) - 1);
end
if isfield(results, 'raw_final_compliance')
    raw_final_compliance = results.raw_final_compliance;
else
    raw_final_compliance = results.final_compliance;
end
if isfield(results, 'raw_final_FCS')
    raw_final_FCS = results.raw_final_FCS;
else
    raw_final_FCS = results.final_FCS;
end
if isfield(results, 'raw_final_improvement_ratio')
    raw_final_improvement_ratio = results.raw_final_improvement_ratio;
else
    raw_final_improvement_ratio = results.improvement_ratio;
end
if isfield(results, 'final_to_best_gap_percent')
    final_to_best_gap_percent = results.final_to_best_gap_percent;
else
    final_to_best_gap_percent = (raw_final_compliance - best_compliance) / max(best_compliance, eps) * 100;
end
if isfield(results, 'theta_only_accept_count')
    theta_only_accept_count = results.theta_only_accept_count;
else
    theta_only_accept_count = NaN;
end
if isfield(results, 'theta_only_reject_count')
    theta_only_reject_count = results.theta_only_reject_count;
else
    theta_only_reject_count = NaN;
end
if isfield(results, 'best_FCS')
    best_FCS = results.best_FCS;
elseif best_iter >= 1 && best_iter <= numel(results.FCS_history)
    best_FCS = results.FCS_history(best_iter);
else
    best_FCS = NaN;
end

metrics = struct();
metrics.run_stamp = run_stamp;
metrics.config = 'fast';
metrics.runtime_seconds = runtime_seconds;
metrics.grid = struct('nelx', params.grid.nelx, 'nely', params.grid.nely, ...
    'dx', params.grid.dx, 'dy', params.grid.dy, 'h', params.grid.h);
metrics.function_paths = function_paths;
metrics.path_spacing = struct( ...
    'target', target_spacing, ...
    'measured_level0_to_level1', measured_spacing, ...
    'relative_error_percent', spacing_error_percent);
metrics.init_stats = struct( ...
    'target_delta_phi', target_delta_phi, ...
    'mean_offset', mean_offset, ...
    'std_offset', std_offset, ...
    'max_error', max_offset_error, ...
    'thin_ratio', thin_ratio, ...
    'num_samples', num_samples);
metrics.optimization = struct( ...
    'initial_compliance', results.compliance_history(1), ...
    'best_compliance', best_compliance, ...
    'best_iter', best_iter, ...
    'raw_final_compliance', raw_final_compliance, ...
    'final_compliance', results.final_compliance, ...
    'raw_final_improvement_ratio', raw_final_improvement_ratio, ...
    'improvement_ratio', results.improvement_ratio, ...
    'initial_FCS', results.FCS_history(1), ...
    'best_FCS', best_FCS, ...
    'raw_final_FCS', raw_final_FCS, ...
    'final_FCS', results.final_FCS, ...
    'final_iter', results.final_iter, ...
    'executed_iter', executed_iter, ...
    'loop_iter_count', loop_iter_count, ...
    'rollback_to_best', rollback_to_best, ...
    'final_to_best_gap_percent', final_to_best_gap_percent, ...
    'early_stop_triggered', early_stop_triggered, ...
    'theta_only_accept_count', theta_only_accept_count, ...
    'theta_only_reject_count', theta_only_reject_count, ...
    'reject_due_next_guard', reject_due_next_guard, ...
    'reject_due_current_guard', reject_due_current_guard, ...
    'reinit_skip_due_next_guard', reinit_skip_due_next_guard, ...
    'reinit_skip_due_current_guard', reinit_skip_due_current_guard, ...
    'compliance_up_steps', compliance_up_steps, ...
    'compliance_total_steps', compliance_total_steps);
metrics.raw_path_quality = get_raw_path_quality_metrics(results, params);
metrics.interface_diagnostics = get_interface_diagnostics_summary(results);
metrics.previous_latest = compare_with_previous(previous_latest_metrics, metrics.raw_path_quality);
metrics.series = struct( ...
    'lsf_change', lsf_change_series(:)', ...
    'dev95', dev95_series(:)', ...
    'mean_off', mean_off_series(:)');

log_path = fullfile(run_dir, 'run_log.txt');
write_text_file(log_path, run_log);

plot_paths = save_key_plots(run_dir, results, params);
if isfield(plot_paths, 'warnings')
    metrics.plot_warnings = plot_paths.warnings;
end
printable = run_printable_postprocess(run_dir, output_root, project_dir, results, params);
metrics.printable = printable;
printability_refinement = run_printability_refinement(run_dir, output_root, project_dir, results, params);
metrics.printability_refinement = printability_refinement;

json_path = fullfile(run_dir, 'baseline_metrics.json');
write_text_file(json_path, jsonencode(metrics));

report_path = fullfile(run_dir, 'baseline_report.md');
report_md = compose_report_md(metrics, log_path, json_path, plot_paths);
write_text_file(report_path, report_md);

copyfile(log_path, fullfile(output_root, 'baseline_latest_run_log.txt'));
copyfile(json_path, fullfile(output_root, 'baseline_latest_metrics.json'));
copyfile(report_path, fullfile(output_root, 'baseline_latest_report.md'));
copy_optional_artifact(plot_paths, 'convergence_png', output_root, 'baseline_latest_convergence.png');
copy_optional_artifact(plot_paths, 'final_lsf_png', output_root, 'baseline_latest_final_lsf.png');
copy_optional_artifact(plot_paths, 'full_domain_png', output_root, 'baseline_latest_full_domain.png');
copy_optional_artifact(plot_paths, 'boundary_compare_png', output_root, 'baseline_latest_boundary_compare.png');
copy_optional_artifact(plot_paths, 'raw_quality_png', output_root, 'baseline_latest_raw_quality.png');

if exist(paths.topology_file, 'file') == 2
    copyfile(paths.topology_file, fullfile(run_dir, 'topo_result_snapshot.mat'));
end
if exist('水平集路径规划.PDF', 'file') == 2
    copyfile('水平集路径规划.PDF', fullfile(run_dir, 'paper_snapshot.pdf'));
end

fprintf('\n=== Baseline Complete ===\n');
fprintf('Run dir: %s\n', run_dir);
fprintf('Report:  %s\n', report_path);
fprintf('JSON:    %s\n', json_path);
fprintf('Log:     %s\n', log_path);
fprintf('Initial compliance: %.10e\n', metrics.optimization.initial_compliance);
fprintf('Best compliance:    %.10e (iter=%d)\n', metrics.optimization.best_compliance, metrics.optimization.best_iter);
fprintf('Raw final compliance: %.10e\n', metrics.optimization.raw_final_compliance);
fprintf('Final compliance:   %.10e\n', metrics.optimization.final_compliance);
fprintf('Improvement ratio:  %.6f %%\n', metrics.optimization.improvement_ratio);
fprintf('Final FCS:          %.6f\n', metrics.optimization.final_FCS);
fprintf('Rollback to best:   %d\n', metrics.optimization.rollback_to_best);
fprintf('Final-to-best gap:  %.6f %%\n', metrics.optimization.final_to_best_gap_percent);
fprintf('Reject (next/current): %d / %d\n', metrics.optimization.reject_due_next_guard, metrics.optimization.reject_due_current_guard);
fprintf('Reinit rollback (next/current): %d / %d\n', ...
    metrics.optimization.reinit_skip_due_next_guard, metrics.optimization.reinit_skip_due_current_guard);
if isfield(metrics, 'printable') && isfield(metrics.printable, 'segment_count') && ~isnan(metrics.printable.segment_count)
    fprintf('Printable segments: %d\n', metrics.printable.segment_count);
    fprintf('Printable mean |dtheta| raw/smooth: %.6f / %.6f deg\n', ...
        metrics.printable.raw_mean_abs_turn_deg, metrics.printable.smooth_mean_abs_turn_deg);
end

function values = extract_numeric_series(raw_text, pattern)
tokens = regexp(raw_text, pattern, 'tokens');
values = nan(numel(tokens), 1);
for i = 1:numel(tokens)
    values(i) = str2double(tokens{i}{1});
end
values = values(isfinite(values));
end

function value = extract_first_numeric(raw_text, pattern)
all_values = extract_numeric_series(raw_text, pattern);
if isempty(all_values)
    value = NaN;
else
    value = all_values(1);
end
end

function write_text_file(path_str, text_content)
fid = fopen(path_str, 'w');
if fid == -1
    error('Failed to open file for writing: %s', path_str);
end
cleanup_fid = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s', text_content);
end

function plot_paths = save_key_plots(run_dir, results, params)
plot_paths = struct( ...
    'convergence_png', fullfile(run_dir, 'convergence.png'), ...
    'final_lsf_png', fullfile(run_dir, 'final_lsf.png'), ...
    'full_domain_png', fullfile(run_dir, 'full_domain_diagnostic.png'), ...
    'boundary_compare_png', fullfile(run_dir, 'boundary_vs_zero_contour.png'), ...
    'raw_quality_png', fullfile(run_dir, 'raw_path_quality_history.png'), ...
    'warnings', {{}});

try
    fig1 = figure('Visible', 'off', 'Position', [80, 80, 1200, 700]);
    yyaxis left;
    plot(results.compliance_history, 'b-', 'LineWidth', 1.6, 'DisplayName', 'Compliance');
    hold on;
    if isfield(results, 'best_iter') && isfield(results, 'best_compliance') && ...
            isfinite(results.best_iter) && results.best_iter >= 1 && results.best_iter <= numel(results.compliance_history)
        plot(results.best_iter, results.best_compliance, 'bo', 'MarkerSize', 7, ...
            'MarkerFaceColor', 'b', 'DisplayName', 'Best Compliance');
    end
    ylabel('Compliance');
    xlabel('Iteration');
    grid on;
    yyaxis right;
    plot(results.FCS_history * 100, 'r-', 'LineWidth', 1.2, 'DisplayName', 'FCS');
    hold on;
    if isfield(results, 'best_iter') && ...
            isfinite(results.best_iter) && results.best_iter >= 1 && results.best_iter <= numel(results.FCS_history)
        if isfield(results, 'best_FCS') && isfinite(results.best_FCS)
            best_fcs_value = results.best_FCS;
        else
            best_fcs_value = results.FCS_history(results.best_iter);
        end
        plot(results.best_iter, best_fcs_value * 100, 'ro', 'MarkerSize', 7, ...
            'MarkerFaceColor', 'r', 'DisplayName', 'Best FCS');
    end
    ylabel('FCS (%)');
    title('Fast Baseline: Compliance and FCS');
    legend('Location', 'best');
    save_png_figure(fig1, plot_paths.convergence_png, 120);
    close(fig1);
catch ME
    plot_paths.warnings{end + 1} = sprintf('convergence plot failed: %s', ME.message); %#ok<AGROW>
end

try
    [x_full, y_full] = build_lsf_grid_coordinates(size(results.lsf), params.grid.dx, params.grid.dy);
    fig2 = figure('Visible', 'off', 'Position', [100, 100, 900, 600]);
    contour(x_full, y_full, results.lsf, 25, 'LineWidth', 0.8);
    hold on;
    contour(x_full, y_full, results.lsf, [0 0], 'r', 'LineWidth', 2);
    axis equal tight;
    title('Fast Baseline: Final Level Set');
    xlabel('x (m)');
    ylabel('y (m)');
    save_png_figure(fig2, plot_paths.final_lsf_png, 120);
    close(fig2);
catch ME
    plot_paths.warnings{end + 1} = sprintf('final lsf plot failed: %s', ME.message); %#ok<AGROW>
end

try
    if ~exist('x_full', 'var') || ~exist('y_full', 'var')
        [x_full, y_full] = build_lsf_grid_coordinates(size(results.lsf), params.grid.dx, params.grid.dy);
    end
    fig3 = figure('Visible', 'off', 'Position', [120, 120, 900, 650]);
    contour(x_full, y_full, results.lsf, 30, 'LineWidth', 0.8);
    hold on;
    contour(x_full, y_full, results.lsf, [0 0], 'r', 'LineWidth', 2);
    draw_material_boundary(results.material_mask_core, params.grid.dx, params.grid.dy);
    axis equal tight;
    title('Full-domain diagnostic');
    xlabel('x (m)');
    ylabel('y (m)');
    save_png_figure(fig3, plot_paths.full_domain_png, 130);
    close(fig3);
catch ME
    plot_paths.warnings{end + 1} = sprintf('full domain plot failed: %s', ME.message); %#ok<AGROW>
end

try
    if isfield(results, 'init_boundary_geometry') && ~isempty(results.init_boundary_geometry)
        if ~exist('x_full', 'var') || ~exist('y_full', 'var')
            [x_full, y_full] = build_lsf_grid_coordinates(size(results.lsf), params.grid.dx, params.grid.dy);
        end
        fig4 = figure('Visible', 'off', 'Position', [140, 140, 900, 650]);
        imagesc([0, params.grid.Lx], [0, params.grid.Ly], flipud(results.material_mask_core));
        set(gca, 'YDir', 'normal');
        colormap(gca, gray);
        hold on;
        plot_geometry(results.init_boundary_geometry);
        zero_contour = extract_contour_segments_local(results.lsf, x_full, y_full, 0);
        for i = 1:numel(zero_contour.segments)
            seg = zero_contour.segments{i};
            plot(seg(:, 1), seg(:, 2), 'r-', 'LineWidth', 1.6);
        end
        title('Boundary reconstruction vs raw zero contour');
        xlabel('x (m)');
        ylabel('y (m)');
        save_png_figure(fig4, plot_paths.boundary_compare_png, 130);
        close(fig4);
    end
catch ME
    plot_paths.warnings{end + 1} = sprintf('boundary compare plot failed: %s', ME.message); %#ok<AGROW>
end

try
    if isfield(results, 'interface_diagnostics') && isfield(results.interface_diagnostics, 'path_quality_history')
        hst = results.interface_diagnostics.path_quality_history;
        fig5 = figure('Visible', 'off', 'Position', [160, 160, 1000, 650]);
        yyaxis left;
        plot(hst.mean_abs_turn_deg, 'LineWidth', 1.5);
        hold on;
        plot(hst.max_abs_kappa, 'LineWidth', 1.2);
        ylabel('turn / kappa');
        yyaxis right;
        plot(hst.parallel_spacing_error_percent, 'LineWidth', 1.2);
        hold on;
        plot(hst.grad_dev_mean, 'LineWidth', 1.2);
        if isfield(hst, 'near_zero_grad_outlier_ratio')
            plot(100 * hst.near_zero_grad_outlier_ratio, 'LineWidth', 1.2);
        end
        ylabel('spacing err / grad dev');
        xlabel('Iteration');
        grid on;
        legend('mean|turn| (deg)', 'max|kappa|', 'spacing err (%)', 'mean||grad|-1|', 'near-zero outlier (%)', ...
            'Location', 'best');
        title('Raw path quality history');
        save_png_figure(fig5, plot_paths.raw_quality_png, 130);
        close(fig5);
    end
catch ME
    plot_paths.warnings{end + 1} = sprintf('raw quality plot failed: %s', ME.message); %#ok<AGROW>
end
end

function save_png_figure(fig, fig_path, resolution)
if nargin < 3 || isempty(resolution)
    resolution = 120;
end

try
    exportgraphics(fig, fig_path, 'Resolution', resolution);
catch exportErr
    warning('run_fast_baseline:exportgraphicsFailed', ...
        'exportgraphics failed for %s (%s); falling back to print.', ...
        fig_path, exportErr.message);
    print(fig, fig_path, '-dpng', sprintf('-r%d', resolution));
end
end

function md = compose_report_md(metrics, log_path, json_path, plot_paths)
printable_block = '';
if isfield(metrics, 'printable') && isfield(metrics.printable, 'segment_count') && ~isnan(metrics.printable.segment_count)
    printable_block = sprintf([ ...
        '## Printable Postprocess\n\n', ...
        '- Segment count: `%d`\n', ...
        '- Raw mean abs turn (deg): `%.6f`\n', ...
        '- Smooth mean abs turn (deg): `%.6f`\n', ...
        '- Raw max kappa: `%.6e`\n', ...
        '- Smooth max kappa: `%.6e`\n', ...
        '- Max deviation from raw (m): `%.6e`\n', ...
        '- Overlay: `%s`\n', ...
        '- Summary: `%s`\n\n'], ...
        metrics.printable.segment_count, ...
        metrics.printable.raw_mean_abs_turn_deg, ...
        metrics.printable.smooth_mean_abs_turn_deg, ...
        metrics.printable.raw_max_abs_kappa, ...
        metrics.printable.smooth_max_abs_kappa, ...
        metrics.printable.max_deviation_from_raw, ...
        metrics.printable.overlay_png, ...
        metrics.printable.summary_txt);
end

raw_quality_block = sprintf([ ...
    '## Raw Path Quality\n\n', ...
    '- Primary gradient scope: `%s`\n', ...
    '- Mean abs turn (deg): `%.6f`\n', ...
    '- Max abs turn (deg): `%.6f`\n', ...
    '- Max abs kappa: `%.6e`\n', ...
    '- Parallel spacing measured: `%.6e`\n', ...
    '- Parallel spacing error (%%): `%.6f`\n', ...
    '- Mean ||grad|-1| near zero: `%.6e`\n', ...
    '- P95 ||grad|-1| near zero: `%.6e`\n', ...
    '- Segment count: `%d`\n\n'], ...
    metrics.raw_path_quality.gradient_primary_scope, ...
    metrics.raw_path_quality.mean_abs_turn_deg, ...
    metrics.raw_path_quality.max_abs_turn_deg, ...
    metrics.raw_path_quality.max_abs_kappa, ...
    metrics.raw_path_quality.parallel_spacing_measured, ...
    metrics.raw_path_quality.parallel_spacing_error_percent, ...
    metrics.raw_path_quality.grad_dev_mean, ...
    metrics.raw_path_quality.grad_dev_p95, ...
    metrics.raw_path_quality.segment_count);

gradient_outlier_block = sprintf([ ...
    '## Gradient Outlier Diagnostics\n\n', ...
    '- Near-zero grad median: `%.6f`\n', ...
    '- Near-zero grad p95: `%.6f`\n', ...
    '- Grad outlier ratio [0.5, 1.5]: `%.6f`\n', ...
    '- Near-zero grad outlier ratio: `%.6f`\n', ...
    '- High-grad outlier ratio: `%.6f`\n', ...
    '- Low-grad outlier ratio: `%.6f`\n', ...
    '- High-grad boundary overlap ratio: `%.6f`\n', ...
    '- Accepted-HJ local reinit count: `%.0f`\n', ...
    '- Refresh count: `%.0f`\n\n'], ...
    metrics.raw_path_quality.near_zero_grad_median, ...
    metrics.raw_path_quality.near_zero_grad_p95, ...
    metrics.raw_path_quality.grad_outlier_ratio_0p5_1p5, ...
    metrics.raw_path_quality.near_zero_grad_outlier_ratio, ...
    metrics.raw_path_quality.high_grad_outlier_ratio, ...
    metrics.raw_path_quality.low_grad_outlier_ratio, ...
    metrics.raw_path_quality.high_grad_boundary_overlap_ratio, ...
    metrics.interface_diagnostics.accepted_hj_reinit_count, ...
    metrics.interface_diagnostics.refresh_count);

interface_block = sprintf([ ...
    '## Material-domain checks\n\n', ...
    '- Outside `phi<=0` count: `%d`\n', ...
    '- Outside nonzero velocity count: `%d`\n', ...
    '- Boundary guard mean ratio: `%.6f`\n', ...
    '- Last reinit method: `%s`\n\n'], ...
    metrics.interface_diagnostics.outside_phi_nonpositive_count, ...
    metrics.interface_diagnostics.outside_velocity_nonzero_count, ...
    metrics.interface_diagnostics.boundary_guard_ratio_mean, ...
    metrics.interface_diagnostics.last_reinit_method);

comparison_block = '';
if isfield(metrics, 'previous_latest') && isstruct(metrics.previous_latest) && metrics.previous_latest.available
    comparison_block = sprintf([ ...
        '## Raw Quality Delta vs Previous Latest\n\n', ...
        '- Previous mean abs turn (deg): `%.6f`\n', ...
        '- Current mean abs turn (deg): `%.6f`\n', ...
        '- Improvement ratio (%%): `%.6f`\n\n'], ...
        metrics.previous_latest.previous_mean_abs_turn_deg, ...
        metrics.previous_latest.current_mean_abs_turn_deg, ...
        metrics.previous_latest.turn_improvement_percent);
end

plot_warning_block = '';
if isfield(metrics, 'plot_warnings') && ~isempty(metrics.plot_warnings)
    warning_lines = cell(1, numel(metrics.plot_warnings));
    for i = 1:numel(metrics.plot_warnings)
        warning_lines{i} = sprintf('- `%s`', metrics.plot_warnings{i});
    end
    plot_warning_block = sprintf('## Plot warnings\n\n%s\n\n', strjoin(warning_lines, newline));
end

md = sprintf([ ...
    '# Fast Baseline Report\n\n', ...
    '- Run stamp: `%s`\n', ...
    '- Config: `%s`\n', ...
    '- Runtime (s): `%.6f`\n\n', ...
    '## Core Metrics\n\n', ...
    '- Initial compliance: `%.10e`\n', ...
    '- Best compliance: `%.10e`\n', ...
    '- Raw final compliance: `%.10e`\n', ...
    '- Final compliance: `%.10e`\n', ...
    '- Raw final improvement ratio (%%): `%.6f`\n', ...
    '- Improvement ratio (%%): `%.6f`\n', ...
    '- Initial FCS: `%.6f`\n', ...
    '- Best FCS: `%.6f`\n', ...
    '- Raw final FCS: `%.6f`\n', ...
    '- Final FCS: `%.6f`\n', ...
    '- Final iter: `%d`\n', ...
    '- Loop iter count: `%d`\n', ...
    '- Rollback to best: `%d`\n', ...
    '- Final-to-best gap (%%): `%.6f`\n', ...
    '- Compliance up-steps: `%d / %d`\n\n', ...
    '## Initialization and Spacing\n\n', ...
    '- Target spacing: `%.10e`\n', ...
    '- Measured spacing (level0->level1): `%.10e`\n', ...
    '- Relative spacing error (%%): `%.6f`\n', ...
    '- Target delta phi: `%.10e`\n', ...
    '- Mean offset: `%.10e`\n', ...
    '- Std offset: `%.10e`\n', ...
    '- Max offset error: `%.10e`\n', ...
    '- Thin ratio: `%.10e`\n', ...
    '- Samples: `%d`\n\n', ...
    '## Diagnostics Series (from log parsing)\n\n', ...
    '- lsf_change samples: `%d`\n', ...
    '- dev95 samples: `%d`\n', ...
    '- mean_off samples: `%d`\n\n', ...
    '%s', ...
    '%s', ...
    '%s', ...
    '%s', ...
    '%s', ...
    '%s', ...
    '## Artifacts\n\n', ...
    '- Run log: `%s`\n', ...
    '- Metrics JSON: `%s`\n', ...
    '- Convergence plot: `%s`\n', ...
    '- Final LSF plot: `%s`\n', ...
    '- Full-domain plot: `%s`\n' ...
    ], ...
    metrics.run_stamp, metrics.config, metrics.runtime_seconds, ...
    metrics.optimization.initial_compliance, metrics.optimization.best_compliance, ...
    metrics.optimization.raw_final_compliance, metrics.optimization.final_compliance, ...
    metrics.optimization.raw_final_improvement_ratio, metrics.optimization.improvement_ratio, ...
    metrics.optimization.initial_FCS, metrics.optimization.best_FCS, ...
    metrics.optimization.raw_final_FCS, metrics.optimization.final_FCS, ...
    metrics.optimization.final_iter, metrics.optimization.loop_iter_count, ...
    metrics.optimization.rollback_to_best, metrics.optimization.final_to_best_gap_percent, ...
    metrics.optimization.compliance_up_steps, metrics.optimization.compliance_total_steps, ...
    metrics.path_spacing.target, metrics.path_spacing.measured_level0_to_level1, ...
    metrics.path_spacing.relative_error_percent, metrics.init_stats.target_delta_phi, ...
    metrics.init_stats.mean_offset, metrics.init_stats.std_offset, ...
    metrics.init_stats.max_error, metrics.init_stats.thin_ratio, ...
    metrics.init_stats.num_samples, ...
    numel(metrics.series.lsf_change), numel(metrics.series.dev95), ...
    numel(metrics.series.mean_off), printable_block, raw_quality_block, gradient_outlier_block, interface_block, ...
    comparison_block, plot_warning_block, log_path, json_path, ...
    get_field_or_default(plot_paths, 'convergence_png', ''), ...
    get_field_or_default(plot_paths, 'final_lsf_png', ''), ...
    get_field_or_default(plot_paths, 'full_domain_png', ''));
end

function printable = run_printable_postprocess(run_dir, output_root, project_dir, results, params)
printable = struct('segment_count', NaN);
try
    addpath(genpath(project_dir), '-begin');
    levels = (-2:2) * params.grid.h;
    printable_dir = fullfile(run_dir, 'printable_paths');
    opts = struct();
    opts.target_ds = params.grid.h / 4;
    opts.smooth_window = 9;
    opts.smooth_iters = 3;
    opts.min_points = 20;
    if isfield(results, 'path_quality_raw')
        opts.raw_path_quality = results.path_quality_raw;
    end
    if isfield(results, 'init_boundary_geometry')
        opts.init_boundary_geometry = results.init_boundary_geometry;
    end

    printable = export_printable_paths_from_lsf(results.lsf, params.grid.dx, params.grid.dy, ...
        levels, printable_dir, results.material_mask_core, opts);

    copyfile(printable.overlay_png, fullfile(output_root, 'printable_paths_latest_overlay.png'));
    copyfile(printable.metrics_csv, fullfile(output_root, 'printable_paths_latest_metrics.csv'));
    copyfile(printable.summary_txt, fullfile(output_root, 'printable_paths_latest_summary.txt'));
    write_text_file(fullfile(output_root, 'printable_paths_latest_summary.json'), jsonencode(printable));
catch ME
    printable.error = ME.message;
    printable.project_dir = project_dir; %#ok<STRNU>
end
end

function metrics = get_raw_path_quality_metrics(results, params)
if isfield(results, 'path_quality_raw') && ~isempty(results.path_quality_raw)
    metrics = results.path_quality_raw;
else
    metrics = compute_raw_path_quality_metrics(results.lsf, params.grid.dx, params.grid.dy, results.material_mask_core);
end
end

function summary = get_interface_diagnostics_summary(results)
summary = struct('outside_phi_nonpositive_count', NaN, ...
    'outside_velocity_nonzero_count', NaN, 'last_reinit_method', '', ...
    'boundary_guard_ratio_mean', NaN, 'accepted_hj_reinit_count', NaN, ...
    'refresh_count', NaN);
if isfield(results, 'interface_diagnostics') && ~isempty(results.interface_diagnostics)
    diag = results.interface_diagnostics;
    if isfield(diag, 'outside_phi_nonpositive_count')
        summary.outside_phi_nonpositive_count = diag.outside_phi_nonpositive_count;
    end
    if isfield(diag, 'outside_velocity_nonzero_count')
        summary.outside_velocity_nonzero_count = diag.outside_velocity_nonzero_count;
    end
    if isfield(diag, 'last_reinit_info') && isstruct(diag.last_reinit_info) && ...
            isfield(diag.last_reinit_info, 'used_method')
        summary.last_reinit_method = diag.last_reinit_info.used_method;
    end
    if isfield(diag, 'boundary_guard_ratio_history')
        summary.boundary_guard_ratio_mean = mean(diag.boundary_guard_ratio_history, 'omitnan');
    end
    if isfield(diag, 'accepted_hj_reinit_count')
        summary.accepted_hj_reinit_count = diag.accepted_hj_reinit_count;
    end
    if isfield(diag, 'refresh_count')
        summary.refresh_count = diag.refresh_count;
    end
end
end

function prev = load_previous_latest_metrics(output_root)
prev = struct('available', false);
json_path = fullfile(output_root, 'baseline_latest_metrics.json');
if exist(json_path, 'file') ~= 2
    return;
end
try
    prev = jsondecode(fileread(json_path));
    prev.available = true;
catch
    prev = struct('available', false);
end
end

function comparison = compare_with_previous(previous_latest_metrics, current_raw_metrics)
comparison = struct('available', false);
if ~isfield(previous_latest_metrics, 'available') || ~previous_latest_metrics.available
    return;
end
if ~isfield(previous_latest_metrics, 'raw_path_quality')
    return;
end
prev_turn = previous_latest_metrics.raw_path_quality.mean_abs_turn_deg;
curr_turn = current_raw_metrics.mean_abs_turn_deg;
if ~isfinite(prev_turn) || prev_turn <= 0 || ~isfinite(curr_turn)
    return;
end
comparison.available = true;
comparison.previous_mean_abs_turn_deg = prev_turn;
comparison.current_mean_abs_turn_deg = curr_turn;
comparison.turn_improvement_percent = (prev_turn - curr_turn) / prev_turn * 100;
end

function draw_material_boundary(material_mask_core, dx, dy)
boundary_mask = bwperim(material_mask_core);
[by, bx] = find(boundary_mask);
plot((bx - 0.5) * dx, (by - 0.5) * dy, 'k--', 'LineWidth', 1.2);
end

function plot_geometry(boundary_geometry)
if isfield(boundary_geometry, 'components')
    components = boundary_geometry.components;
elseif isfield(boundary_geometry, 'segments')
    components = boundary_geometry.segments;
else
    components = {};
end
for i = 1:numel(components)
    seg = components{i};
    plot(seg(:, 1), seg(:, 2), 'b-', 'LineWidth', 1.5);
end
end

function value = get_field_or_default(s, field_name, default_value)
if isfield(s, field_name)
    value = s.(field_name);
else
    value = default_value;
end
end

function copy_optional_artifact(plot_paths, field_name, output_root, dest_name)
if isfield(plot_paths, field_name) && exist(plot_paths.(field_name), 'file') == 2
    copyfile(plot_paths.(field_name), fullfile(output_root, dest_name));
end
end

function refinement = run_printability_refinement(run_dir, output_root, project_dir, results, params)
refinement = struct('refined_selected', false);
try
    addpath(genpath(project_dir), '-begin');
    refine_dir = fullfile(run_dir, 'printability_refinement');
    opts = struct();
    opts.initial_compliance = results.compliance_history(1);
    opts.max_rel_compliance_loss = 0.005;
    opts.theta_reference = results.theta_e;
    opts.baseline_theta = results.theta_e;
    opts.baseline_compliance = results.final_compliance;
    opts.baseline_FCS = results.final_FCS;
    opts.theta_adjust_limit_deg = 0.20;
    opts.levels = (-2:2) * params.grid.h;
    opts.target_ds = params.grid.h / 4;
    opts.path_smooth_window = 9;
    opts.path_smooth_iters = 3;
    opts.min_points = 20;
    opts.blend_values = [0.05, 0.10, 0.15];
    opts.lsf_smooth_iters = [1, 2];
    opts.bandwidth_factors = [2.0, 3.0];
    if isfield(results, 'init_boundary_geometry')
        opts.init_boundary_geometry = results.init_boundary_geometry;
    end

    refinement = refine_lsf_for_printability(results.lsf, params, refine_dir, ...
        results.material_mask_core, opts);

    copyfile(refinement.compare_png, fullfile(output_root, 'printability_refinement_latest_compare.png'));
    copyfile(refinement.best_lsf_png, fullfile(output_root, 'manufacturing_lsf_latest.png'));
    copyfile(refinement.best_overlay_png, fullfile(output_root, 'printability_refinement_latest_overlay.png'));
    copyfile(refinement.best_metrics_csv, fullfile(output_root, 'printability_refinement_latest_metrics.csv'));
    copyfile(refinement.summary_txt, fullfile(output_root, 'printability_refinement_latest_summary.txt'));
    write_text_file(fullfile(output_root, 'printability_refinement_latest_summary.json'), jsonencode(refinement));
catch ME
    refinement.error = ME.message;
    refinement.project_dir = project_dir; %#ok<STRNU>
end
end

function [x_full, y_full] = build_lsf_grid_coordinates(lsf_or_size, dx, dy)
if isnumeric(lsf_or_size) && isvector(lsf_or_size) && numel(lsf_or_size) == 2
    lsf_size = lsf_or_size;
else
    lsf_size = size(lsf_or_size);
end

ny = lsf_size(1);
nx = lsf_size(2);
x_full = ((1:nx) - 1.5) * dx;
y_full = ((1:ny) - 1.5) * dy;
end

function contour_data = extract_contour_segments_local(field, x_vec, y_vec, level)
segments = {};
C = contourc(x_vec, y_vec, field, [level level]);

idx = 1;
while idx < size(C, 2)
    contour_level = C(1, idx);
    npts = C(2, idx);
    next_idx = idx + npts;
    if next_idx <= size(C, 2) && abs(contour_level - level) <= 1e-12 && npts >= 2
        seg = C(:, idx + 1:next_idx)';
        if size(seg, 1) >= 2
            segments{end + 1} = seg; %#ok<AGROW>
        end
    end
    idx = next_idx + 1;
end

contour_data = struct('segments', {segments});
end

function text_content = compose_exception_text(ME)
lines = {ME.message};
for i = 1:numel(ME.stack)
    st = ME.stack(i);
    lines{end + 1} = sprintf('%s:%d %s', st.file, st.line, st.name); %#ok<AGROW>
end
text_content = strjoin(lines, newline);
end

function tf = is_memory_error(ME)
msg = lower(string(ME.message));
id = lower(string(ME.identifier));
tf = contains(msg, "内存不足") || ...
    contains(msg, "insufficient memory") || ...
    contains(msg, "out of memory") || ...
    contains(id, "nomem");
end
