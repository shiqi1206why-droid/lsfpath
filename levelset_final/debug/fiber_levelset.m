function results = fiber_levelset(config_name)
    % 流程入口：主函数仅保留流程调用
    results = run_fiber_levelset_pipeline(config_name);
end

function results = run_fiber_levelset_pipeline(config_name)
    % 纤维路径优化主流程编排（等价重构版）

    if nargin < 1
        config_name = 'fast';
    end

    clc; close all;
    clearvars -except config_name;

    bootstrap_root = fileparts(mfilename('fullpath'));
    addpath(fullfile(bootstrap_root, 'utilities'), '-begin');

    runtime_ctx = fiber_prepare_runtime_context(config_name, bootstrap_root);
    problem_ctx = fiber_load_and_initialize_problem(runtime_ctx);
    iter_out = fiber_run_optimization_iterations(runtime_ctx, problem_ctx);
    results = fiber_finalize_pipeline_results(runtime_ctx, problem_ctx, iter_out);
end
