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

configs = {'default', 'precise'};
template = struct( ...
    'config', '', ...
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
    'early_stop_reason', '');
summaries = repmat(template, numel(configs), 1);

for i = 1:numel(configs)
    cfg = configs{i};
    run_log = evalc('results = fiber_levelset(cfg);'); %#ok<NASGU>
    up_steps = nnz(diff(results.compliance_history) > 0);

    entry = template;
    entry.config = cfg;
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
    summaries(i) = entry;

    log_path = fullfile(output_dir, sprintf('compare_%s_run_log.txt', cfg));
    fid_log = fopen(log_path, 'w');
    if fid_log < 0
        error('无法写入日志文件: %s', log_path);
    end
    cleanup_log = onCleanup(@() fclose(fid_log)); %#ok<NASGU>
    fprintf(fid_log, '%s', run_log);
    clear cleanup_log
end

json_path = fullfile(output_dir, 'default_precise_compare_latest.json');
fid = fopen(json_path, 'w');
if fid < 0
    error('无法写入JSON文件: %s', json_path);
end
cleanup_fid = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s', jsonencode(summaries));

fprintf('=== Comparison Summary ===\n');
for i = 1:numel(summaries)
    s = summaries(i);
    fprintf('%s improvement=%.6f best=%.10e final=%.10e best_iter=%d final_iter=%d executed_iter=%d rollback=%d FCS=%.6f up=%d/%d early_stop=%d\n', ...
        s.config, s.improvement_ratio, s.best_compliance, s.final_compliance, ...
        s.best_iter, s.final_iter, s.executed_iter, s.rollback_to_best, ...
        s.final_FCS, s.compliance_up_steps, s.compliance_total_steps, s.early_stop_triggered);
    if isfield(s, 'early_stop_reason') && ~isempty(s.early_stop_reason)
        fprintf('  reason=%s\n', s.early_stop_reason);
    end
end
fprintf('JSON=%s\n', json_path);
