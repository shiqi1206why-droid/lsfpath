function visualize_results_article(lsf, theta_e, strain_energy, compliance_history, FCS_history, nelx, nely, Lx, Ly, dx, dy)
    % 按参考文章风格展示优化结果

    figure('Position', [50, 50, 1600, 900]);

    subplot(2, 3, 1);
    % 使用索引坐标系统，与其他imagesc子图保持一致（Codex方案）
    [X_idx, Y_idx] = meshgrid(0:nelx+1, 0:nely+1);
    contour(X_idx, Y_idx, lsf, 30, 'LineWidth', 1);
    hold on;
    contour(X_idx, Y_idx, lsf, [0 0], 'r', 'LineWidth', 2);
    title('纤维路径（水平集）');
    xlabel('x方向单元索引'); ylabel('y方向单元索引'); 
    axis equal; axis tight;
    set(gca, 'YDir', 'reverse');  % Y轴向下，与imagesc一致
    colorbar;

    subplot(2, 3, 2);
    [X_e, Y_e] = meshgrid(linspace(dx/2, Lx-dx/2, nelx), linspace(dy/2, Ly-dy/2, nely));
    U_dir = cos(theta_e);
    V_dir = sin(theta_e);
    skip = 2;
    quiver(X_e(1:skip:end, 1:skip:end), Y_e(1:skip:end, 1:skip:end), ...
        U_dir(1:skip:end, 1:skip:end), V_dir(1:skip:end, 1:skip:end), 0.5, 'b');
    axis equal; axis([0 Lx 0 Ly]);
    title('纤维方向'); xlabel('x 坐标'); ylabel('y 坐标');

    subplot(2, 3, 3);
    imagesc(theta_e * 180/pi);
    colormap(hsv);
    colorbar;
    title('纤维角度（度）');
    xlabel('单元 x'); ylabel('单元 y'); axis equal tight;

    subplot(2, 3, 4);
    imagesc(strain_energy);
    colormap(jet);
    colorbar;
    title('应变能密度');
    xlabel('单元 x'); ylabel('单元 y'); axis equal tight;

    subplot(2, 3, 5);
    yyaxis left;
    semilogy(compliance_history, 'b-', 'LineWidth', 2);
    ylabel('柔度');
    xlabel('迭代步');
    grid on;
    yyaxis right;
    plot(FCS_history * 100, 'r-', 'LineWidth', 1);
    ylabel('FCS (%)');
    title('收敛历史');
    legend('柔度', 'FCS', 'Location', 'best');

    subplot(2, 3, 6);
    hold on;
    for i = 1:nely
        for j = 1:nelx
            x_center = (j-0.5) * dx;
            y_center = (i-0.5) * dy;
            angle = theta_e(i, j);
            fibre_len = 0.4 * min(dx, dy);
            x_fibre = x_center + fibre_len * cos(angle) * [-1, 1];
            y_fibre = y_center + fibre_len * sin(angle) * [-1, 1];
            plot(x_fibre, y_fibre, 'b-', 'LineWidth', 0.5);
        end
    end
    axis equal; axis([0 Lx 0 Ly]);
    title('单元纤维走向');
    xlabel('x 坐标'); ylabel('y 坐标');

    sgtitle('边界偏移纤维路径优化结果');
end

