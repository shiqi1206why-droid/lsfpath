% 严格等价对比：当前代码 vs 参考 pre_default 结果
% 用法：直接运行；可选修改 config_name / reference_file

clc; clear; close all;
fprintf('=== Strict Refactor Equivalence Check ===\n');

script_path = mfilename('fullpath');
if isempty(script_path)
    script_path = which('run_refactor_strict_equivalence_check.m');
end
if isempty(script_path)
    error('无法定位 run_refactor_strict_equivalence_check.m');
end
script_dir = fileparts(script_path);           % .../tests/comparison
project_root = fileparts(fileparts(script_dir));

addpath(fullfile(project_root, 'utilities'), '-begin');
project_root = get_project_root(project_root);
cleanup_path = ensure_project_on_path(project_root); %#ok<NASGU>
paths = build_project_paths(project_root);

config_name = 'default';
reference_file = find_latest_reference(paths.refactor_dir, 'pre_default_*.mat');
if isempty(reference_file)
    error('未找到参考文件 pre_default_*.mat，无法做严格等价对比。');
end
fprintf('Reference: %s\n', reference_file);

prev_plot_env = getenv('FIBER_FORCE_ENABLE_PLOTS');
prev_diag_env = getenv('FIBER_FORCE_ENABLE_DIAGNOSTICS');
setenv('FIBER_FORCE_ENABLE_PLOTS', '0');
setenv('FIBER_FORCE_ENABLE_DIAGNOSTICS', '0');

cleanup_env = onCleanup(@() restore_env(prev_plot_env, prev_diag_env));

ref_data = load(reference_file);
if ~isfield(ref_data, 'results')
    error('参考文件缺少 results 变量: %s', reference_file);
end
ref_results = ref_data.results;

fprintf('Running current pipeline: fiber_levelset(''%s'')\n', config_name);
cur_results = fiber_levelset(config_name);

checks = struct();
checks.final_compliance = isequaln(cur_results.final_compliance, ref_results.final_compliance);
checks.final_FCS = isequaln(cur_results.final_FCS, ref_results.final_FCS);
checks.final_iter = isequaln(cur_results.final_iter, ref_results.final_iter);
checks.compliance_history = isequaln(cur_results.compliance_history, ref_results.compliance_history);
checks.FCS_history = isequaln(cur_results.FCS_history, ref_results.FCS_history);
checks.lsf = isequaln(cur_results.lsf, ref_results.lsf);
checks.theta_e = isequaln(cur_results.theta_e, ref_results.theta_e);
checks.accepted_source_history = isequaln(cur_results.accepted_source_history, ref_results.accepted_source_history);
checks.accepted_steps = isequaln(cur_results.accepted_steps, ref_results.accepted_steps);
checks.rejected_steps = isequaln(cur_results.rejected_steps, ref_results.rejected_steps);
checks.reject_due_next_guard = isequaln(cur_results.reject_due_next_guard, ref_results.reject_due_next_guard);
checks.reject_due_current_guard = isequaln(cur_results.reject_due_current_guard, ref_results.reject_due_current_guard);
checks.reinit_trigger_count = isequaln(cur_results.reinit_trigger_count, ref_results.reinit_trigger_count);
checks.reinit_skip_due_reject = isequaln(cur_results.reinit_skip_due_reject, ref_results.reinit_skip_due_reject);
checks.reinit_skip_due_objective = isequaln(cur_results.reinit_skip_due_objective, ref_results.reinit_skip_due_objective);
checks.reinit_skip_due_next_guard = isequaln(cur_results.reinit_skip_due_next_guard, ref_results.reinit_skip_due_next_guard);
checks.reinit_skip_due_current_guard = isequaln(cur_results.reinit_skip_due_current_guard, ref_results.reinit_skip_due_current_guard);

check_names = fieldnames(checks);
check_values = struct2cell(checks);
pass_all = all(cellfun(@(x) logical(x), check_values));

fprintf('\n=== Strict Equivalence Summary ===\n');
for i = 1:numel(check_names)
    fprintf('  %-30s : %d\n', check_names{i}, checks.(check_names{i}));
end
fprintf('  %-30s : %d\n', 'pass_all', pass_all);

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
report = struct();
report.timestamp = timestamp;
report.config_name = config_name;
report.reference_file = reference_file;
report.pass_all = pass_all;
report.checks = checks;

report_mat = fullfile(paths.refactor_dir, sprintf('strict_refactor_comparison_%s.mat', timestamp));
report_txt = fullfile(paths.refactor_dir, sprintf('strict_refactor_comparison_%s.txt', timestamp));
save(report_mat, 'report');

fid = fopen(report_txt, 'w');
if fid == -1
    warning('无法写入文本报告: %s', report_txt);
else
    fprintf(fid, 'Strict refactor equivalence report\n');
    fprintf(fid, 'timestamp: %s\n', timestamp);
    fprintf(fid, 'config: %s\n', config_name);
    fprintf(fid, 'reference: %s\n\n', reference_file);
    for i = 1:numel(check_names)
        fprintf(fid, '%s = %d\n', check_names{i}, checks.(check_names{i}));
    end
    fprintf(fid, '\npass_all = %d\n', pass_all);
    fclose(fid);
end

fprintf('Saved MAT report: %s\n', report_mat);
fprintf('Saved TXT report: %s\n', report_txt);

if ~pass_all
    error('严格等价对比失败：至少一项检查不一致。');
end
fprintf('PASS: 严格等价对比通过。\n');

function latest_file = find_latest_reference(dir_path, pattern)
d = dir(fullfile(dir_path, pattern));
if isempty(d)
    latest_file = '';
    return;
end
[~, idx] = max([d.datenum]);
latest_file = fullfile(d(idx).folder, d(idx).name);
end

function restore_env(prev_plot_env, prev_diag_env)
setenv('FIBER_FORCE_ENABLE_PLOTS', prev_plot_env);
setenv('FIBER_FORCE_ENABLE_DIAGNOSTICS', prev_diag_env);
end
