function enhanced_visualization_check(lsf, material_mask, parallel_paths, struc, nelx, nely, Lx, Ly, dx, dy, delta_phi, init_info)
    % 对初始化的水平集配置进行扩展诊断

    fprintf('正在执行强化初始化诊断...\n');

    figure('Name', '初始化诊断', 'Position', [50, 50, 1600, 1000]);

    % 子图1：拓扑与掩膜
    subplot(2,3,1);
    imagesc(struc);
    colormap(gray);
    axis equal; axis tight;
    title('拓扑与掩膜');
    xlabel('x方向单元索引');
    ylabel('y方向单元索引');
    hold on;
    boundary_mask = bwperim(material_mask);
    [by, bx] = find(boundary_mask);
    plot(bx, by, 'r.', 'MarkerSize', 2);
    legend({'拓扑', '边界'}, 'Location', 'best');

    % 子图2：水平集等值线
    subplot(2,3,2);
    % 使用索引坐标系统，与子图1保持一致（Codex方案）
    [X_idx, Y_idx] = meshgrid(0:nelx+1, 0:nely+1);
    contour(X_idx, Y_idx, lsf, 20, 'LineWidth', 0.5, 'DisplayName', '等值线');
    hold on;
    contour(X_idx, Y_idx, lsf, [0 0], 'r', 'LineWidth', 2, 'DisplayName', '主路径 \phi=0');
    % Parallel offset paths and material boundary removed - not needed
    axis equal; axis tight;
    set(gca, 'YDir', 'reverse');  % Y轴向下，与imagesc一致
    title('水平集等值线');
    xlabel('x方向单元索引');
    ylabel('y方向单元索引');
    legend('Location', 'best');

    % 子图3：梯度模值
    subplot(2,3,3);
    [grad_y, grad_x] = gradient(lsf, dy, dx);
    grad_magnitude = hypot(grad_x, grad_y);
    imagesc(0:nelx+1, 0:nely+1, grad_magnitude);
    set(gca, 'YDir', 'reverse');
    axis equal; axis tight;
    colorbar;
    title('|∇φ|');
    xlabel('x方向单元索引');
    ylabel('y方向单元索引');

    % 子图4：符号分布
    subplot(2,3,4);
    imagesc(0:nelx+1, 0:nely+1, sign(lsf));
    set(gca, 'YDir', 'reverse');
    axis equal; axis tight;
    colormap(gca, jet);
    colorbar;
    title('φ 符号分布');
    xlabel('x方向单元索引');
    ylabel('y方向单元索引');

    % 子图5：窄带梯度直方图
    subplot(2,3,5);
    band_limit = 2 * min(dx, dy);
    band_mask = abs(lsf) <= band_limit;
    band_values = grad_magnitude(band_mask);
    if isempty(band_values)
        band_values = grad_magnitude(:);
    end
    histogram(band_values, 50);
    grid on;
    title('梯度模直方图');
    xlabel('|∇φ|');
    ylabel('频次');

    % 子图6：偏移距离样本
    subplot(2,3,6);
    if isfield(init_info, 'distance_samples') && ~isempty(init_info.distance_samples)
        plot(init_info.distance_samples, 'bo', 'MarkerSize', 4, 'DisplayName', '测得值');
        hold on;
        yline(delta_phi, 'r--', 'LineWidth', 1.5, 'DisplayName', '目标值');
        grid on;
        xlabel('样本索引');
        ylabel('到边界距离 (m)');
        title('偏移距离样本');
        legend('Location', 'best');
    else
        text(0.5, 0.5, '暂无样本', 'HorizontalAlignment', 'center');
        title('偏移距离样本');
        axis off;
    end

    sgtitle('初始化诊断');

    fprintf('\n=== 初始化质量报告 ===\n');
    fprintf('目标Δφ = %.4f m，抽样均值 = %.4f m，最大误差 = %.4f m\n', ...
        delta_phi, init_info.mean_offset, init_info.max_error);
    sign_consistency = check_sign_consistency(lsf);
    fprintf('符号一致性: %.2f%%\n', sign_consistency * 100);
    fprintf('薄壁单元比例: %.2f%%\n', init_info.thin_ratio * 100);

    if isfield(parallel_paths, 'positive_levels') && ~isempty(parallel_paths.positive_levels)
        spacing_error = 0;
        for level = parallel_paths.positive_levels
            actual = verify_path_spacing(lsf, 0, level, dx, dy);
            spacing_error = spacing_error + abs(actual - level) / max(level, eps);
        end
        spacing_error = spacing_error / numel(parallel_paths.positive_levels);
        fprintf('平行路径平均间距误差: %.2f%%\n', spacing_error * 100);
    end
    fprintf('====================================\n\n');
end

