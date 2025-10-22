% 性能基准测试 - 对比优化前后的性能提升
% 注意：完整测试需要较长时间（约20-30分钟）

clc; clear; close all;

fprintf('=== 纤维路径优化 - 性能基准测试 ===\n\n');
fprintf('⚠️  本测试将运行多个完整优化，预计需要20-30分钟\n');
fprintf('    建议使用fast模式进行快速验证\n\n');

% 用户选择
choice = input('选择测试模式:\n  1) 完整测试 (default配置, ~15分钟)\n  2) 快速测试 (fast配置, ~3分钟)\n  3) 调试测试 (debug配置, <1分钟)\n请输入 (1/2/3): ', 's');

switch choice
    case '1'
        config = 'default';
        fprintf('\n选择: 完整测试模式\n');
    case '2'
        config = 'fast';
        fprintf('\n选择: 快速测试模式\n');
    case '3'
        config = 'debug';
        fprintf('\n选择: 调试测试模式\n');
    otherwise
        config = 'fast';
        fprintf('\n默认: 快速测试模式\n');
end

% 添加路径
cd ..
script_dir = pwd;
addpath(genpath(script_dir));

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
figure('Name', '优化收敛分析', 'Position', [100, 100, 1200, 800]);

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
initial = results_optimized.compliance_history(1);
improvement = (initial - results_optimized.compliance_history) / initial * 100;
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

cd tests

