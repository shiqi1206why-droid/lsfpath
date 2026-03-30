clc; close all;

script_dir = fileparts(mfilename('fullpath'));
project_dir = fileparts(fileparts(script_dir));
addpath(fullfile(project_dir, 'utilities'), '-begin');
project_root = get_project_root(project_dir);
cleanup_path = ensure_project_on_path(project_root); %#ok<NASGU>
set(0, 'DefaultFigureVisible', 'off');

paths = build_project_paths(project_root);
output_dir = fullfile(paths.tests_dir, 'comparison');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

configs = {'fast', 'default', 'precise'};
template = struct( ...
    'config', '', ...
    'initial_compliance', NaN, ...
    'improvement_ratio', NaN, ...
    'best_compliance', NaN, ...
    'final_compliance', NaN, ...
    'best_iter', NaN, ...
    'final_iter', NaN, ...
    'executed_iter', NaN, ...
    'rollback_to_best', false, ...
    'final_FCS', NaN, ...
    'raw_final_compliance', NaN, ...
    'final_to_best_gap_percent', NaN, ...
    'compliance_up_steps', NaN, ...
    'compliance_total_steps', NaN, ...
    'early_stop_triggered', false, ...
    'early_stop_reason', '', ...
    'compliance_history', [], ...
    'normalized_history', []);
summaries = repmat(template, numel(configs), 1);

for i = 1:numel(configs)
    cfg = configs{i};
    fprintf('=== Running %s ===\n', cfg);
    run_log = evalc('results = fiber_levelset(cfg);'); %#ok<NASGU>
    up_steps = nnz(diff(results.compliance_history) > 0);
    base_compliance = results.compliance_history(1);

    entry = template;
    entry.config = cfg;
    entry.initial_compliance = base_compliance;
    entry.improvement_ratio = results.improvement_ratio;
    entry.best_compliance = results.best_compliance;
    entry.final_compliance = results.final_compliance;
    entry.best_iter = results.best_iter;
    entry.final_iter = results.final_iter;
    entry.executed_iter = results.executed_iter;
    entry.rollback_to_best = results.rollback_to_best;
    entry.final_FCS = results.final_FCS;
    entry.raw_final_compliance = results.raw_final_compliance;
    entry.final_to_best_gap_percent = results.final_to_best_gap_percent;
    entry.compliance_up_steps = up_steps;
    entry.compliance_total_steps = numel(diff(results.compliance_history));
    entry.early_stop_triggered = results.early_stop_triggered;
    if isfield(results, 'early_stop_reason') && ~isempty(results.early_stop_reason)
        entry.early_stop_reason = results.early_stop_reason;
    end
    entry.compliance_history = results.compliance_history(:)';
    entry.normalized_history = results.compliance_history(:)' ./ base_compliance;
    summaries(i) = entry;

    log_path = fullfile(output_dir, sprintf('mode_compare_%s_run_log.txt', cfg));
    fid_log = fopen(log_path, 'w');
    if fid_log < 0
        error('无法写入日志文件: %s', log_path);
    end
    cleanup_log = onCleanup(@() fclose(fid_log)); %#ok<NASGU>
    fprintf(fid_log, '%s', run_log);
    clear cleanup_log
end

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
artifact_dir = fullfile(paths.baseline_dir, ['mode_compare_' timestamp]);
if ~exist(artifact_dir, 'dir')
    mkdir(artifact_dir);
end

fig = figure('Color', 'w', 'Position', [100, 100, 1200, 520]);
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
colors = lines(numel(configs));

nexttile;
hold on;
for i = 1:numel(summaries)
    s = summaries(i);
    plot(0:(numel(s.compliance_history) - 1), s.compliance_history, ...
        'LineWidth', 1.8, 'Color', colors(i, :), 'DisplayName', s.config);
    plot(s.best_iter, s.best_compliance, 'o', 'MarkerSize', 7, ...
        'LineWidth', 1.2, 'Color', colors(i, :), 'HandleVisibility', 'off');
end
grid on;
xlabel('Iteration');
ylabel('Compliance');
title('Raw Compliance History');
legend('Location', 'best');

nexttile;
hold on;
for i = 1:numel(summaries)
    s = summaries(i);
    plot(0:(numel(s.normalized_history) - 1), s.normalized_history, ...
        'LineWidth', 1.8, 'Color', colors(i, :), 'DisplayName', s.config);
    plot(s.best_iter, s.best_compliance / s.initial_compliance, 'o', ...
        'MarkerSize', 7, 'LineWidth', 1.2, 'Color', colors(i, :), ...
        'HandleVisibility', 'off');
end
grid on;
xlabel('Iteration');
ylabel('C / C_0');
title('Normalized Compliance History');
legend('Location', 'best');

sgtitle('Fiber Level-Set Mode Comparison');

fig_path = fullfile(artifact_dir, 'mode_comparison_convergence.png');
exportgraphics(fig, fig_path, 'Resolution', 200);
close(fig);

json_path = fullfile(artifact_dir, 'mode_comparison_summary.json');
fid_json = fopen(json_path, 'w');
if fid_json < 0
    error('无法写入JSON文件: %s', json_path);
end
cleanup_json = onCleanup(@() fclose(fid_json)); %#ok<NASGU>
fprintf(fid_json, '%s', jsonencode(summaries));
clear cleanup_json

mat_path = fullfile(artifact_dir, 'mode_comparison_results.mat');
save(mat_path, 'summaries');

table_path = fullfile(artifact_dir, 'mode_comparison_summary.txt');
fid_txt = fopen(table_path, 'w');
if fid_txt < 0
    error('无法写入文本文件: %s', table_path);
end
cleanup_txt = onCleanup(@() fclose(fid_txt)); %#ok<NASGU>
fprintf(fid_txt, 'config\timprovement_ratio\tbest_iter\tfinal_iter\tfinal_compliance\tfinal_FCS\trollback_to_best\tcompliance_up_steps\tcompliance_total_steps\n');
for i = 1:numel(summaries)
    s = summaries(i);
    fprintf(fid_txt, '%s\t%.6f\t%d\t%d\t%.10e\t%.6f\t%d\t%d\t%d\n', ...
        s.config, s.improvement_ratio, s.best_iter, s.final_iter, ...
        s.final_compliance, s.final_FCS, s.rollback_to_best, ...
        s.compliance_up_steps, s.compliance_total_steps);
end

fprintf('=== Mode Comparison Summary ===\n');
for i = 1:numel(summaries)
    s = summaries(i);
    fprintf('%s improvement=%.6f best=%.10e final=%.10e best_iter=%d final_iter=%d rollback=%d FCS=%.6f up=%d/%d\n', ...
        s.config, s.improvement_ratio, s.best_compliance, s.final_compliance, ...
        s.best_iter, s.final_iter, s.rollback_to_best, s.final_FCS, ...
        s.compliance_up_steps, s.compliance_total_steps);
end
fprintf('FIG=%s\n', fig_path);
fprintf('JSON=%s\n', json_path);
fprintf('MAT=%s\n', mat_path);
fprintf('TXT=%s\n', table_path);
