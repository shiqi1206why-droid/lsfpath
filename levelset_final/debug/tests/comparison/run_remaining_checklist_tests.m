clc; close all;

script_path = mfilename('fullpath');
if isempty(script_path)
    script_path = which('run_remaining_checklist_tests.m');
end
script_dir = fileparts(script_path);
project_dir = fileparts(fileparts(script_dir));
addpath(fullfile(project_dir, 'utilities'), '-begin');
project_root = get_project_root(project_dir);
cleanup_path = ensure_project_on_path(project_root); %#ok<NASGU>
set(0, 'DefaultFigureVisible', 'off');

% Reduce UI overhead during batch checks
setenv('FIBER_FORCE_ENABLE_PLOTS', '0');
setenv('FIBER_FORCE_ENABLE_DIAGNOSTICS', '0');

paths = build_project_paths(project_root);
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
artifact_dir = fullfile(paths.refactor_dir, ['checklist_remaining_' timestamp]);
if ~exist(artifact_dir, 'dir')
    mkdir(artifact_dir);
end

log_path = fullfile(artifact_dir, 'run_log.txt');
diary(log_path);
diary on;

fprintf('=== CHECKLIST REMAINING TESTS START ===\n');

% Step 1/4: strict fast exact
setenv('FIBER_GRADIENT_CHAIN_MODE', 'exact');
setenv('FIBER_FAST_ACCEPTANCE_TOL', '0');
setenv('FIBER_FAST_CURRENT_TOL', '0');
setenv('FIBER_FAST_REINIT_CURRENT_TOL', '0');
setenv('FIBER_FAST_ENABLE_CURRENT_GUARD', '1');
setenv('FIBER_FAST_ENABLE_REINIT_CURRENT_GUARD', '1');

fprintf('RUN_FAST_STRICT_EXACT\n');
results_fast = fiber_levelset('fast');

% clear fast overrides
setenv('FIBER_GRADIENT_CHAIN_MODE', '');
setenv('FIBER_FAST_ACCEPTANCE_TOL', '');
setenv('FIBER_FAST_CURRENT_TOL', '');
setenv('FIBER_FAST_REINIT_CURRENT_TOL', '');
setenv('FIBER_FAST_ENABLE_CURRENT_GUARD', '');
setenv('FIBER_FAST_ENABLE_REINIT_CURRENT_GUARD', '');

% Step 2: shadow cosine check
fprintf('RUN_DEBUG_SHADOW\n');
results_debug = fiber_levelset('debug');

% Step 5/6: default run
fprintf('RUN_DEFAULT\n');
results_default = fiber_levelset('default');

% Step 7 (user-modified): run precise once
fprintf('RUN_PRECISE\n');
results_precise = fiber_levelset('precise');

summary = struct();
summary.artifact_dir = artifact_dir;
summary.log_path = log_path;

% Step 1/4
summary.step1_pass = results_fast.accepted_steps > 0;
summary.step4_pass = summary.step1_pass;
summary.step1_accepted_steps = results_fast.accepted_steps;
summary.step1_rejected_steps = results_fast.rejected_steps;

% Step 2
cos_vals = [];
if isfield(results_debug, 'gradient_chain_history') && isstruct(results_debug.gradient_chain_history) && ...
        isfield(results_debug.gradient_chain_history, 'cosine_similarity_history')
    cos_vals = results_debug.gradient_chain_history.cosine_similarity_history(:);
elseif isfield(results_debug, 'interface_diagnostics') && isstruct(results_debug.interface_diagnostics) && ...
        isfield(results_debug.interface_diagnostics, 'gradient_chain') && ...
        isfield(results_debug.interface_diagnostics.gradient_chain, 'cosine_similarity_history')
    cos_vals = results_debug.interface_diagnostics.gradient_chain.cosine_similarity_history(:);
end
cos_vals = cos_vals(isfinite(cos_vals));
if isempty(cos_vals)
    summary.step2_pass = false;
    summary.step2_cos_median = NaN;
    summary.step2_cos_min = NaN;
    summary.step2_cos_mean = NaN;
else
    summary.step2_cos_median = median(cos_vals);
    summary.step2_cos_min = min(cos_vals);
    summary.step2_cos_mean = mean(cos_vals);
    summary.step2_pass = summary.step2_cos_median > 0.7;
end

% Step 5
hist_default = results_default.compliance_history(:);
N = min(50, numel(hist_default));
dh = diff(hist_default(1:N));
summary.step5_firstN = N;
summary.step5_up_steps = sum(dh > 0);
summary.step5_pass = (summary.step5_up_steps == 0);

% Step 6
summary.step6_rollback_count = double(results_default.rollback_to_best);
summary.step6_pass_rollback_lt3 = summary.step6_rollback_count < 3;

fuse_limit = results_default.params.opt.theta_only_fuse_limit;
src = results_default.accepted_source_history;
if isstring(src)
    src = cellstr(src);
end
if ~iscell(src)
    src = {};
end
max_streak = 0;
cur_streak = 0;
for k = 1:numel(src)
    s = src{k};
    if isstring(s)
        s = char(s);
    end
    if ischar(s) && strcmpi(strtrim(s), 'theta_only')
        cur_streak = cur_streak + 1;
        max_streak = max(max_streak, cur_streak);
    else
        cur_streak = 0;
    end
end
summary.step6_fuse_limit = fuse_limit;
summary.step6_theta_only_max_streak = max_streak;
summary.step6_fuse_condition_seen = max_streak >= fuse_limit;

% Step 7
summary.step7_precise_executed = true;
summary.step7_precise_final_iter = results_precise.final_iter;
summary.step7_precise_executed_iter = results_precise.executed_iter;
summary.step7_precise_final_compliance = results_precise.final_compliance;
summary.step7_precise_final_FCS = results_precise.final_FCS;
summary.step7_precise_chain_mode = char(results_precise.params.gradient.chain_mode);
summary.step7_precise_accepted_steps = results_precise.accepted_steps;
summary.step7_precise_rejected_steps = results_precise.rejected_steps;

summary.fast = struct( ...
    'final_compliance', results_fast.final_compliance, ...
    'final_iter', results_fast.final_iter, ...
    'executed_iter', results_fast.executed_iter, ...
    'accepted_steps', results_fast.accepted_steps, ...
    'rejected_steps', results_fast.rejected_steps, ...
    'chain_mode', char(results_fast.params.gradient.chain_mode));
summary.debug = struct( ...
    'final_compliance', results_debug.final_compliance, ...
    'final_iter', results_debug.final_iter, ...
    'executed_iter', results_debug.executed_iter, ...
    'accepted_steps', results_debug.accepted_steps, ...
    'rejected_steps', results_debug.rejected_steps, ...
    'chain_mode', char(results_debug.params.gradient.chain_mode));
summary.default = struct( ...
    'final_compliance', results_default.final_compliance, ...
    'final_iter', results_default.final_iter, ...
    'executed_iter', results_default.executed_iter, ...
    'accepted_steps', results_default.accepted_steps, ...
    'rejected_steps', results_default.rejected_steps, ...
    'rollback_to_best', logical(results_default.rollback_to_best), ...
    'chain_mode', char(results_default.params.gradient.chain_mode));
summary.precise = struct( ...
    'final_compliance', results_precise.final_compliance, ...
    'final_iter', results_precise.final_iter, ...
    'executed_iter', results_precise.executed_iter, ...
    'accepted_steps', results_precise.accepted_steps, ...
    'rejected_steps', results_precise.rejected_steps, ...
    'chain_mode', char(results_precise.params.gradient.chain_mode));

save(fullfile(artifact_dir, 'summary.mat'), 'summary', 'results_fast', 'results_debug', 'results_default', 'results_precise');

summary_txt = fullfile(artifact_dir, 'summary.txt');
fid = fopen(summary_txt, 'w');
if fid < 0
    error('无法写入摘要文件: %s', summary_txt);
end
cleanup_summary = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'artifact_dir=%s\n', artifact_dir);
fprintf(fid, 'log_path=%s\n', log_path);

fprintf(fid, '\n[step1]\npass=%d\naccepted_steps=%d\nrejected_steps=%d\n', ...
    summary.step1_pass, summary.step1_accepted_steps, summary.step1_rejected_steps);

fprintf(fid, '\n[step2]\npass=%d\ncos_median=%.6f\ncos_min=%.6f\ncos_mean=%.6f\n', ...
    summary.step2_pass, summary.step2_cos_median, summary.step2_cos_min, summary.step2_cos_mean);

fprintf(fid, '\n[step4]\npass=%d\n', summary.step4_pass);

fprintf(fid, '\n[step5]\npass=%d\nfirstN=%d\nup_steps=%d\n', ...
    summary.step5_pass, summary.step5_firstN, summary.step5_up_steps);

fprintf(fid, '\n[step6]\npass_rollback_lt3=%d\nrollback_count=%d\nfuse_limit=%d\nmax_theta_only_streak=%d\nfuse_condition_seen=%d\n', ...
    summary.step6_pass_rollback_lt3, summary.step6_rollback_count, ...
    summary.step6_fuse_limit, summary.step6_theta_only_max_streak, summary.step6_fuse_condition_seen);

fprintf(fid, '\n[step7]\nprecise_executed=%d\nchain_mode=%s\nfinal_iter=%d\nexecuted_iter=%d\nfinal_compliance=%.12e\nfinal_FCS=%.12f\naccepted_steps=%d\nrejected_steps=%d\n', ...
    summary.step7_precise_executed, summary.step7_precise_chain_mode, ...
    summary.step7_precise_final_iter, summary.step7_precise_executed_iter, ...
    summary.step7_precise_final_compliance, summary.step7_precise_final_FCS, ...
    summary.step7_precise_accepted_steps, summary.step7_precise_rejected_steps);

clear cleanup_summary

json_path = fullfile(artifact_dir, 'summary.json');
fid_json = fopen(json_path, 'w');
if fid_json < 0
    error('无法写入JSON文件: %s', json_path);
end
cleanup_json = onCleanup(@() fclose(fid_json)); %#ok<NASGU>
fprintf(fid_json, '%s', jsonencode(summary));
clear cleanup_json

diary off;

fprintf('CHECKLIST_ARTIFACT=%s\n', artifact_dir);
fprintf('CHECKLIST_SUMMARY=%s\n', summary_txt);
