% 测试所有配置模式
% 验证配置系统、参数验证、日志系统的正确性

clc; clear; close all;

fprintf('=== 测试所有配置模式 ===\n\n');

% 添加路径
cd ..
script_dir = pwd;
addpath(genpath(script_dir));

%% 测试1: 默认配置
fprintf('【测试1】默认配置 (default)\n');
try
    params_default = get_fiber_optimization_params('default');
    validate_params(params_default);
    fprintf('  ✅ 默认配置验证通过\n');
    fprintf('     - max_iter: %d\n', params_default.opt.max_iter);
    fprintf('     - log_level: %s\n', params_default.debug.log_level);
    fprintf('     - enable_plots: %d\n', params_default.debug.enable_plots);
catch ME
    fprintf('  ❌ 失败: %s\n', ME.message);
end

%% 测试2: 快速配置
fprintf('\n【测试2】快速配置 (fast)\n');
try
    params_fast = get_fiber_optimization_params('fast');
    validate_params(params_fast);
    fprintf('  ✅ 快速配置验证通过\n');
    fprintf('     - max_iter: %d\n', params_fast.opt.max_iter);
    fprintf('     - transition_iter: %d\n', params_fast.levelset.transition_iter);
    fprintf('     - enable_plots: %d\n', params_fast.debug.enable_plots);
catch ME
    fprintf('  ❌ 失败: %s\n', ME.message);
end

%% 测试3: 精确配置
fprintf('\n【测试3】精确配置 (precise)\n');
try
    params_precise = get_fiber_optimization_params('precise');
    validate_params(params_precise);
    fprintf('  ✅ 精确配置验证通过\n');
    fprintf('     - max_iter: %d\n', params_precise.opt.max_iter);
    fprintf('     - save_checkpoints: %d\n', params_precise.debug.save_checkpoints);
    fprintf('     - checkpoint_interval: %d\n', params_precise.debug.checkpoint_interval);
catch ME
    fprintf('  ❌ 失败: %s\n', ME.message);
end

%% 测试4: 调试配置
fprintf('\n【测试4】调试配置 (debug)\n');
try
    params_debug = get_fiber_optimization_params('debug');
    validate_params(params_debug);
    fprintf('  ✅ 调试配置验证通过\n');
    fprintf('     - max_iter: %d\n', params_debug.opt.max_iter);
    fprintf('     - log_level: %s\n', params_debug.debug.log_level);
    fprintf('     - log_interval: %d\n', params_debug.debug.log_interval);
catch ME
    fprintf('  ❌ 失败: %s\n', ME.message);
end

%% 测试5: 材料参数库
fprintf('\n【测试5】材料参数库\n');
try
    mat_carbon = get_material_params('carbon_fiber');
    fprintf('  ✅ 碳纤维材料加载成功\n');
    fprintf('     - 名称: %s\n', mat_carbon.name);
    fprintf('     - E_L/E_T: %.1f\n', mat_carbon.E_L/mat_carbon.E_T);
    
    mat_glass = get_material_params('glass_fiber');
    fprintf('  ✅ 玻璃纤维材料加载成功\n');
    fprintf('     - 名称: %s\n', mat_glass.name);
    fprintf('     - E_L/E_T: %.1f\n', mat_glass.E_L/mat_glass.E_T);
catch ME
    fprintf('  ❌ 失败: %s\n', ME.message);
end

%% 测试6: 日志系统
fprintf('\n【测试6】日志系统\n');
try
    % 测试不同日志级别
    params_info = get_fiber_optimization_params('default');
    params_info.debug.log_level = 'INFO';
    
    fprintf('  测试INFO级别（应显示3条）:\n');
    log_message('DEBUG', params_info, '    这条DEBUG不应显示');
    log_message('INFO', params_info, '    这条INFO应该显示');
    log_message('WARN', params_info, '    这条WARN应该显示');
    log_message('ERROR', params_info, '    这条ERROR应该显示');
    
    fprintf('  ✅ 日志系统工作正常\n');
catch ME
    fprintf('  ❌ 失败: %s\n', ME.message);
end

%% 测试7: 参数验证（错误情况）
fprintf('\n【测试7】参数验证（错误检测）\n');
try
    params_bad = get_fiber_optimization_params('default');
    
    % 故意设置错误参数
    params_bad.opt.max_iter = -10;  % 错误：负迭代次数
    validate_params(params_bad);
    fprintf('  ❌ 应该捕获错误但未捕获\n');
catch ME
    fprintf('  ✅ 成功捕获错误参数: %s\n', ME.message);
end

%% 测试总结
fprintf('\n=== 测试完成 ===\n');
fprintf('所有配置系统测试已完成！\n');
fprintf('如果所有测试都显示✅，说明配置系统工作正常。\n');

cd tests

