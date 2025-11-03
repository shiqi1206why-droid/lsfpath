% 验证脚本：测试节点灵敏度聚合的窄带过滤功能
% 修改时间：2025-10-31
% 目的：验证窄带外的 node_sensitivity 是否为零

clear; clc;

%% 1. 构造测试案例
nelx = 20;
nely = 20;
dx = 0.1;
dy = 0.1;

% 构造简单的水平集场（圆形）
[X, Y] = meshgrid(0:dx:nelx*dx, 0:dy:nely*dy);
center_x = nelx * dx / 2;
center_y = nely * dy / 2;
radius = 0.3;
lsf = sqrt((X - center_x).^2 + (Y - center_y).^2) - radius;

% 随机单元灵敏度
element_sensitivity = randn(nely, nelx);

% 虚拟角度场（不会被使用）
theta_field = zeros(nely, nelx);

% 定义窄带
h = min(dx, dy);
band_mask = abs(lsf) <= 1.0 * h;

%% 2. 调用聚合函数
fprintf('=== 测试节点灵敏度聚合（窄带过滤）===\n');
node_sensitivity = aggregate_node_sensitivity(element_sensitivity, theta_field, lsf, nelx, nely, dx, dy, band_mask);

%% 3. 验证窄带外灵敏度为零
% 扩展窄带掩码到节点（ghost cells）
band_mask_extended = false(size(lsf));
band_mask_extended(2:end-1, 2:end-1) = band_mask(2:end-1, 2:end-1);

% 检查窄带外的灵敏度
outside_band = ~band_mask_extended;
sensitivity_outside = node_sensitivity(outside_band);

max_outside = max(abs(sensitivity_outside));
nonzero_outside = nnz(sensitivity_outside);

fprintf('\n=== 验证结果 ===\n');
fprintf('窄带内节点数: %d\n', nnz(band_mask));
fprintf('窄带外节点数: %d\n', nnz(outside_band));
fprintf('窄带外最大灵敏度: %.6e\n', max_outside);
fprintf('窄带外非零节点数: %d\n', nonzero_outside);

% 可视化窄带外的灵敏度分布
tolerance = 1e-14;  % 数值精度容忍度
if nonzero_outside > 0 && max_outside > tolerance
    fprintf('⚠️  警告：窄带外存在非零灵敏度！\n');
    
    % 显示非零分布
    figure('Name', '窄带外灵敏度分布');
    subplot(1,2,1);
    imagesc(abs(lsf));
    colorbar;
    title('水平集场 |φ|');
    hold on;
    contour(lsf, [0 0], 'r', 'LineWidth', 2);
    
    subplot(1,2,2);
    imagesc(abs(node_sensitivity));
    colorbar;
    title('节点灵敏度 |∂E/∂φ_i|');
    hold on;
    contour(lsf, [0 0], 'r', 'LineWidth', 2);
    contour(lsf, [-h, h], 'g--', 'LineWidth', 1);
else
    fprintf('✓  通过：窄带外灵敏度全部为零！\n');
end

%% 4. 验证窄带内灵敏度非零
inside_band = band_mask_extended;
sensitivity_inside = node_sensitivity(inside_band);
nonzero_inside = nnz(sensitivity_inside);
ratio_nonzero = nonzero_inside / nnz(inside_band) * 100;

fprintf('\n=== 窄带内统计 ===\n');
fprintf('窄带内非零节点数: %d / %d (%.1f%%)\n', nonzero_inside, nnz(inside_band), ratio_nonzero);
fprintf('窄带内最大灵敏度: %.6e\n', max(abs(sensitivity_inside)));
fprintf('窄带内平均灵敏度: %.6e\n', mean(abs(sensitivity_inside)));

if nonzero_inside > 0
    fprintf('✓  通过：窄带内存在非零灵敏度！\n');
else
    fprintf('⚠️  警告：窄带内全部为零，可能存在问题！\n');
end

fprintf('\n=== 测试完成 ===\n');


