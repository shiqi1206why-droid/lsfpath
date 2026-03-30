% P1/P2重构严格等价对比脚本
% 对比当前代码与参考pre_default结果，要求关键字段严格一致（isequaln）

clc; clear; close all;
fprintf('=== P1/P2 Strict Equivalence Check ===\n');

script_path = mfilename('fullpath');
if isempty(script_path)
    script_path = which('run_p1p2_strict_equivalence_check.m');
end
if isempty(script_path)
    error('无法定位 run_p1p2_strict_equivalence_check.m');
end
script_dir = fileparts(script_path);
project_root = fileparts(fileparts(script_dir));

addpath(fullfile(project_root, 'utilities'), '-begin');
project_root = get_project_root(project_root);
cleanup_path = ensure_project_on_path(project_root); %#ok<NASGU>
paths = build_project_paths(project_root);

config_name = 'default';
reference_file = find_latest_reference(paths.refactor_dir, 'pre_default_*.mat');
if isempty(reference_file)
    error('未找到参考文件 pre_default_*.mat，无法进行严格等价对比。');
end
fprintf('Reference: %s\n', reference_file);

prev_plot_env = getenv('FIBER_FORCE_ENABLE_PLOTS');
prev_diag_env = getenv('FIBER_FORCE_ENABLE_DIAGNOSTICS');
setenv('FIBER_FORCE_ENABLE_PLOTS', '0');
setenv('FIBER_FORCE_ENABLE_DIAGNOSTICS', '0');
cleanup_env = onCleanup(@() restore_env(prev_plot_env, prev_diag_env)); %#ok<NASGU>

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

if isfield(cur_results, 'hj_trial_compliance_history') && isfield(ref_results, 'hj_trial_compliance_history')
    checks.hj_trial_compliance_history = isequaln(cur_results.hj_trial_compliance_history, ref_results.hj_trial_compliance_history);
end
if isfield(cur_results, 'reinit_trial_compliance_history') && isfield(ref_results, 'reinit_trial_compliance_history')
    checks.reinit_trial_compliance_history = isequaln(cur_results.reinit_trial_compliance_history, ref_results.reinit_trial_compliance_history);
end

check_names = fieldnames(checks);
check_values = struct2cell(checks);
pass_all = all(cellfun(@(x) logical(x), check_values));

report = struct();
report.timestamp = datestr(now, 'yyyymmdd_HHMMSS');
report.config_name = config_name;
report.reference_file = reference_file;
report.pass_all = pass_all;
report.checks = checks;
report.final_compliance_rel_diff = rel_diff(cur_results.final_compliance, ref_results.final_compliance);
report.final_FCS_rel_diff = rel_diff(cur_results.final_FCS, ref_results.final_FCS);
report.final_iter_diff = cur_results.final_iter - ref_results.final_iter;
report.compliance_history_max_abs_diff = max_abs_diff(cur_results.compliance_history, ref_results.compliance_history);
report.FCS_history_max_abs_diff = max_abs_diff(cur_results.FCS_history, ref_results.FCS_history);
report.lsf_inf_diff = inf_norm_diff(cur_results.lsf, ref_results.lsf);
report.theta_inf_diff = inf_norm_diff(cur_results.theta_e, ref_results.theta_e);

fprintf('\n=== P1/P2 Strict Equivalence Summary ===\n');
for i = 1:numel(check_names)
    fprintf('  %-35s : %d\n', check_names{i}, checks.(check_names{i}));
end
fprintf('  %-35s : %d\n', 'pass_all', pass_all);
fprintf('  %-35s : %.16e\n', 'final_compliance_rel_diff', report.final_compliance_rel_diff);
fprintf('  %-35s : %.16e\n', 'final_FCS_rel_diff', report.final_FCS_rel_diff);
fprintf('  %-35s : %d\n', 'final_iter_diff', report.final_iter_diff);
fprintf('  %-35s : %.16e\n', 'compliance_history_max_abs_diff', report.compliance_history_max_abs_diff);
fprintf('  %-35s : %.16e\n', 'FCS_history_max_abs_diff', report.FCS_history_max_abs_diff);
fprintf('  %-35s : %.16e\n', 'lsf_inf_diff', report.lsf_inf_diff);
fprintf('  %-35s : %.16e\n', 'theta_inf_diff', report.theta_inf_diff);

stamp = report.timestamp;
report_mat = fullfile(paths.refactor_dir, sprintf('p1p2_equivalence_%s.mat', stamp));
report_txt = fullfile(paths.refactor_dir, sprintf('p1p2_equivalence_%s.txt', stamp));
save(report_mat, 'report');

fid = fopen(report_txt, 'w');
if fid == -1
    warning('无法写入文本报告: %s', report_txt);
else
    fprintf(fid, 'P1/P2 strict equivalence report\n');
    fprintf(fid, 'timestamp: %s\n', stamp);
    fprintf(fid, 'config: %s\n', config_name);
    fprintf(fid, 'reference: %s\n\n', reference_file);
    for i = 1:numel(check_names)
        fprintf(fid, '%s = %d\n', check_names{i}, checks.(check_names{i}));
    end
    fprintf(fid, '\nfinal_compliance_rel_diff = %.16e\n', report.final_compliance_rel_diff);
    fprintf(fid, 'final_FCS_rel_diff = %.16e\n', report.final_FCS_rel_diff);
    fprintf(fid, 'final_iter_diff = %d\n', report.final_iter_diff);
    fprintf(fid, 'compliance_history_max_abs_diff = %.16e\n', report.compliance_history_max_abs_diff);
    fprintf(fid, 'FCS_history_max_abs_diff = %.16e\n', report.FCS_history_max_abs_diff);
    fprintf(fid, 'lsf_inf_diff = %.16e\n', report.lsf_inf_diff);
    fprintf(fid, 'theta_inf_diff = %.16e\n', report.theta_inf_diff);
    fprintf(fid, '\npass_all = %d\n', pass_all);
    fclose(fid);
end

fprintf('Saved MAT report: %s\n', report_mat);
fprintf('Saved TXT report: %s\n', report_txt);

if ~pass_all
    error('P1/P2严格等价对比失败：至少一项检查不一致。');
end
fprintf('PASS: P1/P2严格等价对比通过。\n');

function latest_file = find_latest_reference(dir_path, pattern)
d = dir(fullfile(dir_path, pattern));
if isempty(d)
    latest_file = '';
    return;
end
[~, idx] = max([d.datenum]);
latest_file = fullfile(d(idx).folder, d(idx).name);
end

function d = rel_diff(a, b)
d = abs(a - b) / max(1e-30, abs(b));
end

function d = max_abs_diff(a, b)
d = max(abs(a(:) - b(:)));
end

function d = inf_norm_diff(a, b)
d = norm(a(:) - b(:), inf);
end

function restore_env(prev_plot_env, prev_diag_env)
setenv('FIBER_FORCE_ENABLE_PLOTS', prev_plot_env);
setenv('FIBER_FORCE_ENABLE_DIAGNOSTICS', prev_diag_env);
end
