% 测试矢量化角度平滑函数
% 验证数值精度和性能提升

clc; clear; close all;

fprintf('=== 矢量化角度平滑函数测试 ===\n\n');

% 测试参数
nely = 50;
nelx = 80;
eta = 0.10;
num_iters = 2;

% 创建测试角度场（随机角度）
rng(42);  % 固定随机种子确保可重复
theta_e = rand(nely, nelx) * pi;

%% 1. 原始实现（三重嵌套循环）
fprintf('1. 运行原始实现（三重嵌套循环）...\n');
tic;
z_old = cos(2*theta_e) + 1i*sin(2*theta_e);
for k = 1:num_iters
    z_smooth = zeros(size(z_old));
    for i = 2:nely-1
        for j = 2:nelx-1
            z_smooth(i,j) = z_old(i,j) + eta * (z_old(i-1,j) + z_old(i+1,j) + z_old(i,j-1) + z_old(i,j+1) - 4*z_old(i,j));
        end
    end
    % 边界外推
    z_smooth(1,:) = z_smooth(2,:);
    z_smooth(end,:) = z_smooth(end-1,:);
    z_smooth(:,1) = z_smooth(:,2);
    z_smooth(:,end) = z_smooth(:,end-1);
    % 归一化
    z_old = z_smooth ./ max(abs(z_smooth), 1e-12);
end
time_old = toc;
result_old = 0.5 * angle(z_old);
result_old = mod(result_old, pi);

fprintf('  耗时: %.6f 秒\n', time_old);

%% 2. 矢量化实现
fprintf('\n2. 运行矢量化实现（conv2）...\n');

% 添加路径
script_dir = fileparts(mfilename('fullpath'));
project_dir = fileparts(script_dir);
addpath(genpath(project_dir));

tic;
z_new = angle_smooth_vectorized(theta_e, eta, num_iters);
time_new = toc;
result_new = 0.5 * angle(z_new);
result_new = mod(result_new, pi);

fprintf('  耗时: %.6f 秒\n', time_new);

%% 3. 性能对比
speedup = time_old / time_new;
fprintf('\n=== 性能对比 ===\n');
fprintf('原始耗时: %.6f 秒\n', time_old);
fprintf('优化耗时: %.6f 秒\n', time_new);
fprintf('性能提升: %.1f%%\n', (1 - time_new/time_old) * 100);
fprintf('加速比: %.2fx\n', speedup);

%% 4. 数值误差验证
max_error = max(abs(result_old(:) - result_new(:)));
mean_error = mean(abs(result_old(:) - result_new(:)));
relative_error = max_error / (max(abs(result_old(:))) + 1e-12);

fprintf('\n=== 数值一致性验证 ===\n');
fprintf('最大绝对误差: %.2e (目标 < 1e-6)\n', max_error);
fprintf('平均绝对误差: %.2e\n', mean_error);
fprintf('最大相对误差: %.2e\n', relative_error);

%% 5. 验收判定
fprintf('\n=== 验收判定 ===\n');

% 性能要求：加速比 >= 3倍（预期80%加速 = 5倍）
if speedup >= 3.0
    fprintf('✅ 性能测试通过（加速比 %.2fx >= 3.0x）\n', speedup);
else
    fprintf('❌ 性能测试未通过（加速比 %.2fx < 3.0x）\n', speedup);
end

% 精度要求：最大误差 < 1e-6
if max_error < 1e-6
    fprintf('✅ 数值精度测试通过（误差 %.2e < 1e-6）\n', max_error);
else
    fprintf('❌ 数值精度测试未通过（误差 %.2e >= 1e-6）\n', max_error);
end

%% 6. 可视化对比
figure('Name', '角度平滑对比', 'Position', [100, 100, 1200, 400]);

subplot(1,3,1);
imagesc(theta_e * 180/pi);
colorbar;
title('原始角度场（度）');
axis equal tight;

subplot(1,3,2);
imagesc(result_old * 180/pi);
colorbar;
title('原始实现平滑结果（度）');
axis equal tight;

subplot(1,3,3);
imagesc(result_new * 180/pi);
colorbar;
title('矢量化实现平滑结果（度）');
axis equal tight;

% 差异图
figure('Name', '数值误差分析', 'Position', [150, 150, 800, 600]);
imagesc(abs(result_old - result_new) * 180/pi);
colorbar;
title(sprintf('绝对误差分布（度）\\newline最大误差=%.2e rad', max_error));
axis equal tight;

fprintf('\n测试完成！\n');

