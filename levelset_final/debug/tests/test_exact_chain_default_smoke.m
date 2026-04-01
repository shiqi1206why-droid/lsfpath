clear; clc; close all;
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
addpath(fullfile(project_root, 'utilities'), '-begin');
project_root = get_project_root(project_root);
cleanup_path = ensure_project_on_path(project_root); %#ok<NASGU>

params = get_fiber_optimization_params('fast');
params_default = get_fiber_optimization_params('default');
assert(strcmpi(params_default.gradient.chain_mode, 'legacy'), ...
    'default 配置在 exact 审计通过前应保持 legacy。');
assert(strcmpi(params.gradient.chain_mode, 'legacy'), ...
    'fast 配置在 exact 审计通过前应保持 legacy。');

results = fiber_levelset('fast');
assert(isfield(results, 'final_compliance') && isfinite(results.final_compliance), ...
    'baseline smoke 运行未返回有效 final_compliance。');
assert(isfield(results, 'interface_diagnostics') && ...
    isfield(results.interface_diagnostics, 'gradient_chain'), ...
    'baseline smoke 运行缺少 gradient_chain 诊断。');

fprintf('PASS: baseline chain config smoke test.\n');
