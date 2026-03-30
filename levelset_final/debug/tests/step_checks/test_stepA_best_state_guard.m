% StepA补充测试：历史最优状态守护应保证输出不劣于历史最优柔度

clc; clear; close all;
fprintf('=== StepA best-state guard test ===\n');

script_path = mfilename('fullpath');
if isempty(script_path)
    script_path = which('test_stepA_best_state_guard.m');
end
if isempty(script_path)
    error('无法定位 test_stepA_best_state_guard.m');
end

script_dir = fileparts(script_path);             % .../tests/step_checks
project_dir = fileparts(fileparts(script_dir));  % .../debug

addpath(fullfile(project_dir, 'utilities'), '-begin');
project_root = get_project_root(project_dir);
cleanup_path = ensure_project_on_path(project_root); %#ok<NASGU>

results = fiber_levelset('fast');
hist = results.compliance_history(:);
hist = hist(isfinite(hist) & hist > 0);

assert(~isempty(hist), 'compliance_history为空。');
[min_hist, min_idx] = min(hist);

tol = max(1e-12, 1e-8 * abs(min_hist));
fprintf('executed_iter=%d final_iter=%d best_iter=%d\n', ...
    results.executed_iter, results.final_iter, results.best_iter);
fprintf('min(history)=%.10e (iter=%d) best=%.10e final=%.10e rollback=%d\n', ...
    min_hist, min_idx, results.best_compliance, results.final_compliance, results.rollback_to_best);

assert(results.final_compliance <= min_hist + tol, ...
    '输出最终柔度高于历史最优柔度，守护失效。');
assert(results.best_iter >= 1 && results.best_iter <= results.final_iter, ...
    'best_iter越界。');
assert(abs(results.best_compliance - min_hist) <= max(1e-12, 1e-6 * abs(min_hist)), ...
    'best_compliance与历史最小柔度不一致。');

fprintf('PASS: StepA best-state guard test.\n');
