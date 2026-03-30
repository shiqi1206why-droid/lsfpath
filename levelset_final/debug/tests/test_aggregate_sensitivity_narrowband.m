% 验证脚本：测试节点灵敏度聚合的窄带过滤功能
% 修改时间：2025-10-31
% 目的：验证窄带外的 node_sensitivity 是否为零

clear; clc;

% 稳定路径初始化，确保调用 debug 副本函数
script_path = mfilename('fullpath');
if isempty(script_path)
    script_path = which('test_aggregate_sensitivity_narrowband.m');
end
if isempty(script_path)
    error('无法定位 test_aggregate_sensitivity_narrowband.m');
end
script_dir = fileparts(script_path);
project_dir = fileparts(script_dir);
addpath(fullfile(project_dir, 'utilities'), '-begin');
project_root = get_project_root(project_dir);
cleanup_path = ensure_project_on_path(project_root); %#ok<NASGU>

which_agg = which('aggregate_node_sensitivity');
fprintf('aggregate_node_sensitivity path: %s\n', which_agg);

% 1. 构造测试案例
nelx = 20;
nely = 20;
dx = 0.1;
dy = 0.1;

% 构造简单的水平集场（圆形，含 ghost cells，尺寸为 (nely+2) x (nelx+2)）
[X, Y] = meshgrid(0:dx:(nelx+1)*dx, 0:dy:(nely+1)*dy);
center_x = (nelx + 1) * dx / 2;
center_y = (nely + 1) * dy / 2;
radius = 0.3;
lsf = sqrt((X - center_x).^2 + (Y - center_y).^2) - radius;

% 随机单元灵敏度
element_sensitivity = randn(nely, nelx);

% 虚拟角度场（不会被使用）
theta_field = zeros(nely, nelx);

% 定义窄带
h = min(dx, dy);
band_mask = abs(lsf) <= 1.0 * h;
material_mask_core = false(nely, nelx);
material_mask_core(4:end-3, 5:end-4) = true;
material_mask_full = false(nely + 2, nelx + 2);
material_mask_full(2:end-1, 2:end-1) = material_mask_core;
material_mask_full(1, :) = material_mask_full(2, :);
material_mask_full(end, :) = material_mask_full(end-1, :);
material_mask_full(:, 1) = material_mask_full(:, 2);
material_mask_full(:, end) = material_mask_full(:, end-1);
active_mask = band_mask & material_mask_full;

% 2. 调用聚合函数
fprintf('=== 测试节点灵敏度聚合（窄带过滤）===\n');
node_sensitivity = aggregate_node_sensitivity(element_sensitivity, theta_field, lsf, nelx, nely, dx, dy, active_mask);

% 3. 验证窄带外灵敏度为零
outside_band = ~active_mask;
sensitivity_outside = node_sensitivity(outside_band);

max_outside = max(abs(sensitivity_outside));
nonzero_outside = nnz(sensitivity_outside);

fprintf('\n=== 验证结果 ===\n');
fprintf('窄带内节点数: %d\n', nnz(active_mask));
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

% 4. 验证窄带内灵敏度非零
inside_band = active_mask;
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

% 5. 验证材料域外窄带无灵敏度泄漏
void_band = band_mask & ~material_mask_full;
void_leak = nnz(abs(node_sensitivity(void_band)) > tolerance);

fprintf('\n=== 材料域外泄漏检查 ===\n');
fprintf('材料域外窄带节点数: %d\n', nnz(void_band));
fprintf('材料域外非零灵敏度节点数: %d\n', void_leak);
assert(void_leak == 0, '材料域外窄带存在灵敏度泄漏。');
fprintf('✓  通过：材料域外窄带无灵敏度泄漏！\n');

fprintf('\n=== 测试完成 ===\n');

