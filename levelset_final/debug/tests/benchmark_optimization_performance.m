% 性能基准测试 - 对比优化前后的性能提升
% 注意：完整测试需要较长时间（约20-30分钟）

clc; close all;
clearvars -except config enable_plots;

fprintf('=== 纤维路径优化 - 性能基准测试 ===\n\n');
fprintf('⚠️  本测试将运行多个完整优化，预计需要20-30分钟\n');
fprintf('    建议使用fast模式进行快速验证\n\n');

% 稳定路径初始化（不依赖当前工作目录）
script_path = mfilename('fullpath');
if isempty(script_path)
    script_path = which('benchmark_optimization_performance.m');
end
if isempty(script_path)
    error('无法定位 benchmark_optimization_performance.m');
end
script_dir = fileparts(script_path);
project_dir = fileparts(script_dir);
addpath(fullfile(project_dir, 'utilities'), '-begin');
project_root = get_project_root(project_dir);
cleanup_path = ensure_project_on_path(project_root); %#ok<NASGU>

% 非交互模式选择（优先变量，其次环境变量，最后默认fast）
if exist('config', 'var') && ~isempty(config)
    cfg = lower(string(config));
else
    env_mode = string(getenv('BENCHMARK_MODE'));
    if strlength(env_mode) == 0
        cfg = "fast";
    else
        cfg = lower(strtrim(env_mode));
    end
end

switch char(cfg)
    case {'1', 'default'}
        config = 'default';
    case {'2', 'fast'}
        config = 'fast';
    case {'3', 'debug'}
        config = 'debug';
    otherwise
        config = 'fast';
end
fprintf('运行模式: %s\n', config);

% 图像开关（默认关闭，适合批量/CI）
if ~exist('enable_plots', 'var') || isempty(enable_plots)
    enable_plots = false;
end
env_plots = string(getenv('BENCHMARK_PLOTS'));
if strlength(env_plots) > 0
    enable_plots = any(strcmpi(strtrim(env_plots), ["1", "true", "on", "yes"]));
end

%% 运行优化版本
fprintf('\n=== 运行优化版本 (%s配置) ===\n', config);
tic;
results_optimized = fiber_levelset(config);
time_optimized = toc;

fprintf('\n优化版本完成！\n');
fprintf('  耗时: %.2f 秒 (%.2f 分钟)\n', time_optimized, time_optimized/60);
fprintf('  实际迭代: %d 次\n', results_optimized.final_iter);
fprintf('  最终柔度: %.6e\n', results_optimized.final_compliance);
fprintf('  最终FCS: %.2f%%\n', results_optimized.final_FCS * 100);
fprintf('  柔度改善: %.2f%%\n', results_optimized.improvement_ratio);

%% 性能分析
fprintf('\n=== 性能统计 ===\n');
fprintf('总耗时: %.2f 秒\n', time_optimized);
fprintf('平均每次迭代: %.3f 秒\n', time_optimized / results_optimized.final_iter);
fprintf('配置模式: %s\n', config);

%% 收敛性分析
initial = results_optimized.compliance_history(1);
improvement = (initial - results_optimized.compliance_history) / initial * 100;

if enable_plots
    fig = figure('Name', '优化收敛分析', 'Position', [100, 100, 1200, 800], 'Visible', 'on');

    % 柔度历史
    subplot(2,2,1);
    plot(results_optimized.compliance_history, 'b-', 'LineWidth', 1.5);
    grid on;
    title('柔度收敛历史');
    xlabel('迭代次数');
    ylabel('柔度');

    % FCS历史
    subplot(2,2,2);
    plot(results_optimized.FCS_history * 100, 'r-', 'LineWidth', 1.5);
    grid on;
    title('纤维连续性评分');
    xlabel('迭代次数');
    ylabel('FCS (%)');

    % 柔度改善率
    subplot(2,2,3);
    plot(improvement, 'g-', 'LineWidth', 1.5);
    grid on;
    title('累积柔度改善率');
    xlabel('迭代次数');
    ylabel('改善率 (%)');
    yline(0, 'k--', 'LineWidth', 1);

    % 统计摘要
    subplot(2,2,4);
    axis off;
    text(0.1, 0.9, '优化统计摘要', 'FontSize', 14, 'FontWeight', 'bold');
    text(0.1, 0.75, sprintf('配置模式: %s', config), 'FontSize', 11);
    text(0.1, 0.65, sprintf('总耗时: %.2f 秒', time_optimized), 'FontSize', 11);
    text(0.1, 0.55, sprintf('迭代次数: %d', results_optimized.final_iter), 'FontSize', 11);
    text(0.1, 0.45, sprintf('初始柔度: %.4e', initial), 'FontSize', 11);
    text(0.1, 0.35, sprintf('最终柔度: %.4e', results_optimized.final_compliance), 'FontSize', 11);
    text(0.1, 0.25, sprintf('柔度改善: %.2f%%', results_optimized.improvement_ratio), 'FontSize', 11);
    text(0.1, 0.15, sprintf('最终FCS: %.2f%%', results_optimized.final_FCS * 100), 'FontSize', 11);

    drawnow;
    if exist('fig', 'var') && isvalid(fig)
        close(fig);
    end
else
    fprintf('已跳过收敛图绘制（BENCHMARK_PLOTS 未开启）。\n');
end

%% 保存结果
save_file = sprintf('benchmark_results_%s_%s.mat', config, datestr(now, 'yyyymmdd_HHMMSS'));
save(save_file, 'results_optimized', 'time_optimized', 'config');
fprintf('\n结果已保存至: %s\n', save_file);

%% 性能评估
fprintf('\n=== 性能评估 ===\n');

% 根据配置评估
switch config
    case 'default'
        expected_time_per_iter = 10;  % 秒/次（预估）
        fprintf('预期每次迭代耗时: ~%.1f 秒\n', expected_time_per_iter);
    case 'fast'
        expected_time_per_iter = 8;
        fprintf('预期每次迭代耗时: ~%.1f 秒\n', expected_time_per_iter);
    case 'debug'
        expected_time_per_iter = 12;
        fprintf('预期每次迭代耗时: ~%.1f 秒\n', expected_time_per_iter);
end

actual_time_per_iter = time_optimized / results_optimized.final_iter;
fprintf('实际每次迭代耗时: %.2f 秒\n', actual_time_per_iter);

if actual_time_per_iter < expected_time_per_iter * 0.7
    fprintf('✅ 性能优秀！（比预期快%.0f%%）\n', ...
        (1 - actual_time_per_iter/expected_time_per_iter) * 100);
elseif actual_time_per_iter < expected_time_per_iter * 1.2
    fprintf('✅ 性能正常\n');
else
    fprintf('⚠️  性能偏慢（比预期慢%.0f%%）\n', ...
        (actual_time_per_iter/expected_time_per_iter - 1) * 100);
end

%% 优化质量评估
fprintf('\n=== 优化质量评估 ===\n');

if results_optimized.improvement_ratio > 10
    fprintf('✅ 优秀：柔度降低 > 10%%\n');
elseif results_optimized.improvement_ratio > 5
    fprintf('✅ 良好：柔度降低 > 5%%\n');
elseif results_optimized.improvement_ratio > 0
    fprintf('⚠️  一般：柔度降低 > 0%%\n');
else
    fprintf('❌ 警告：柔度增加了！可能需要更多迭代或调整参数\n');
end

if results_optimized.final_FCS > 0.85
    fprintf('✅ FCS优秀 (%.1f%% > 85%%)\n', results_optimized.final_FCS * 100);
elseif results_optimized.final_FCS > 0.75
    fprintf('✅ FCS良好 (%.1f%% > 75%%)\n', results_optimized.final_FCS * 100);
else
    fprintf('⚠️  FCS偏低 (%.1f%% < 75%%)\n', results_optimized.final_FCS * 100);
end

fprintf('\n=== 测试完成 ===\n');

