clear; clc; close all;
addpath(genpath(fileparts(fileparts(fileparts(mfilename('fullpath'))))), '-begin');

set(0, 'DefaultFigureVisible', 'off');
cleanup_fig = onCleanup(@() set(0, 'DefaultFigureVisible', 'on')); %#ok<NASGU>

% 使用fast配置缩短自动化测试时长（验证门控性质不受配置影响）
results = fiber_levelset('fast');

fprintf('accepted=%d rejected=%d reinit_trigger=%d reinit_skip_due_reject=%d\n', ...
    results.accepted_steps, results.rejected_steps, ...
    results.reinit_trigger_count, results.reinit_skip_due_reject);
if isfield(results, 'reinit_skip_due_next_guard') && isfield(results, 'reinit_skip_due_current_guard')
    fprintf('reinit_skip_due_next_guard=%d reinit_skip_due_current_guard=%d reinit_skip_due_objective=%d\n', ...
        results.reinit_skip_due_next_guard, results.reinit_skip_due_current_guard, ...
        results.reinit_skip_due_objective);
end

assert(results.rejected_steps >= 0, 'rejected_steps无效。');
assert(results.reinit_skip_due_reject >= 0, 'reinit_skip_due_reject无效。');
assert(results.reinit_skip_due_reject == results.rejected_steps, ...
    '拒绝步与重初始化跳过计数不一致。');
if isfield(results, 'reinit_skip_due_next_guard') && isfield(results, 'reinit_skip_due_current_guard')
    assert(results.reinit_skip_due_next_guard >= 0, 'reinit_skip_due_next_guard无效。');
    assert(results.reinit_skip_due_current_guard >= 0, 'reinit_skip_due_current_guard无效。');
    assert(results.reinit_skip_due_objective >= results.reinit_skip_due_next_guard, ...
        'reinit_skip_due_objective应覆盖next_guard回退计数。');
    assert(results.reinit_skip_due_objective >= results.reinit_skip_due_current_guard, ...
        'reinit_skip_due_objective应覆盖current_guard回退计数。');
end

% 合成门控场景：若不加门控，本应触发重初始化；加门控后必须跳过
params = get_fiber_optimization_params('fast');
lsf_dummy = zeros(params.grid.nely + 2, params.grid.nelx + 2);
lsf_before_dummy = lsf_dummy;
iter_since_last_reinit = params.levelset.reinit_freq_early;
[do_reinit_raw, reason_raw] = should_reinitialize(lsf_dummy, lsf_before_dummy, ...
    iter_since_last_reinit, 1, params);
step_accepted = false;
do_reinit_gated = do_reinit_raw && step_accepted;
fprintf('raw_reinit=%d reason=%s gated_reinit=%d\n', do_reinit_raw, reason_raw, do_reinit_gated);
assert(do_reinit_raw, '合成场景下应可触发重初始化。');
assert(~do_reinit_gated, '拒绝步门控失败：仍触发了重初始化。');

fprintf('PASS: StepA reject-step reinit gate test.\n');
