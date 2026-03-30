function visualize_results_article(lsf, theta_e, strain_energy, compliance_history, FCS_history, ...
    nelx, nely, ~, ~, dx, dy, material_mask_core, material_mask_full, diag_info)
    % 绘制论文主结果图：raw路径质量为主，printable smoothing仅作补充分析

    if nargin < 12 || isempty(material_mask_core)
        material_mask_core = true(nely, nelx);
    end
    if nargin < 13 || isempty(material_mask_full)
        material_mask_full = expand_material_mask_to_full_local(material_mask_core);
    end
    if nargin < 14 || isempty(diag_info)
        diag_info = struct();
    end

    [raw_history, raw_metrics] = resolve_raw_quality_info(diag_info);
    plot_paths = resolve_plot_paths(diag_info);

    [x_full, y_full, x_core, y_core] = get_lsf_grid_coordinates(size(lsf), dx, dy);
    [X_core, Y_core] = meshgrid(x_core, y_core);

    lsf_masked = lsf;
    lsf_masked(~material_mask_full) = NaN;

    fig = figure('Position', [40, 40, 1680, 980]);
    tiledlayout(2, 4, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    hold on;
    contour(x_full, y_full, lsf_masked, 24, 'LineWidth', 0.6, 'Color', [0.60, 0.60, 0.60]);
    contour(x_full, y_full, lsf_masked, [0, 0], 'r', 'LineWidth', 2);
    plot_material_boundary(material_mask_core, dx, dy, '--', [0.10, 0.10, 0.10], 1.1);
    title('材料域内水平集');
    xlabel('x (m)');
    ylabel('y (m)');
    axis equal tight;

    nexttile;
    hold on;
    U_dir = cos(theta_e);
    V_dir = sin(theta_e);
    U_dir(~material_mask_core) = NaN;
    V_dir(~material_mask_core) = NaN;
    skip = 2;
    quiver(X_core(1:skip:end, 1:skip:end), Y_core(1:skip:end, 1:skip:end), ...
        U_dir(1:skip:end, 1:skip:end), V_dir(1:skip:end, 1:skip:end), 0.55, ...
        'Color', [0.10, 0.30, 0.75]);
    plot_material_boundary(material_mask_core, dx, dy, '-', [0.15, 0.15, 0.15], 0.8);
    title('纤维方向场');
    xlabel('x (m)');
    ylabel('y (m)');
    axis equal tight;

    nexttile;
    theta_plot = theta_e * 180 / pi;
    theta_plot(~material_mask_core) = NaN;
    imagesc(x_core, y_core, theta_plot);
    axis equal tight;
    set(gca, 'YDir', 'normal');
    colormap(gca, hsv);
    colorbar;
    title('纤维角度 (deg)');
    xlabel('x (m)');
    ylabel('y (m)');

    nexttile;
    strain_plot = strain_energy;
    strain_plot(~material_mask_core) = NaN;
    imagesc(x_core, y_core, strain_plot);
    axis equal tight;
    set(gca, 'YDir', 'normal');
    colormap(gca, turbo);
    colorbar;
    title('应变能密度');
    xlabel('x (m)');
    ylabel('y (m)');

    nexttile;
    yyaxis left;
    semilogy(compliance_history, 'LineWidth', 1.8, 'Color', [0.12, 0.34, 0.72]);
    ylabel('Compliance');
    yyaxis right;
    plot(FCS_history * 100, 'LineWidth', 1.3, 'Color', [0.78, 0.18, 0.15]);
    ylabel('FCS (%)');
    grid on;
    xlabel('Iteration');
    title('收敛历史');

    nexttile;
    hold on;
    if ~isempty(raw_history)
        plot(raw_history.mean_abs_turn_deg, 'LineWidth', 1.5, 'DisplayName', 'mean |turn| (deg)');
        plot(raw_history.parallel_spacing_error_percent, 'LineWidth', 1.3, 'DisplayName', 'spacing err (%)');
        plot(raw_history.grad_dev_mean, 'LineWidth', 1.3, 'DisplayName', 'mean ||grad|-1|');
        if isfield(raw_history, 'max_abs_kappa')
            plot(raw_history.max_abs_kappa, 'LineWidth', 1.2, 'DisplayName', 'max curvature');
        end
        legend('Location', 'best');
    else
        text(0.5, 0.5, 'No path-quality history', 'HorizontalAlignment', 'center');
        axis off;
    end
    grid on;
    xlabel('Iteration');
    title('Raw 路径质量历史');

    nexttile;
    hold on;
    contour(x_full, y_full, lsf, 24, 'LineWidth', 0.6, 'Color', [0.60, 0.60, 0.75]);
    contour(x_full, y_full, lsf, [0, 0], 'r', 'LineWidth', 2);
    plot_material_boundary(material_mask_core, dx, dy, '--', [0.10, 0.10, 0.10], 1.1);
    title('全域未裁剪 LSF 诊断');
    xlabel('x (m)');
    ylabel('y (m)');
    axis equal tight;

    nexttile;
    hold on;
    if isfield(diag_info, 'init_boundary_geometry') && ~isempty(diag_info.init_boundary_geometry)
        geom = diag_info.init_boundary_geometry;
        if isfield(geom, 'pixel_boundary') && ~isempty(geom.pixel_boundary)
            plot(geom.pixel_boundary(:, 1), geom.pixel_boundary(:, 2), '.', ...
                'Color', [0.35, 0.35, 0.35], 'MarkerSize', 7, 'DisplayName', 'pixel boundary');
        end
        boundary_segments = collect_boundary_segments(geom);
        for k = 1:numel(boundary_segments)
            seg = boundary_segments{k};
            plot(seg(:, 1), seg(:, 2), '-', 'Color', [0.10, 0.45, 0.20], ...
                'LineWidth', 1.1, 'DisplayName', ternary_label(k == 1, 'subcell boundary', ''));
        end
    else
        plot_material_boundary(material_mask_core, dx, dy, '-', [0.35, 0.35, 0.35], 1.0);
    end
    zero_geom = extract_levelset_geometry(lsf, dx, dy, 0, material_mask_full);
    for k = 1:numel(zero_geom.segments)
        seg = zero_geom.segments{k};
        plot(seg(:, 1), seg(:, 2), 'r-', 'LineWidth', 1.6, ...
            'DisplayName', ternary_label(k == 1, 'zero contour', ''));
    end
    title('像素边界 / 重建边界 / 零线');
    xlabel('x (m)');
    ylabel('y (m)');
    axis equal tight;
    legend('Location', 'best');

    sgtitle('边界偏移纤维路径优化结果');
    save_png_figure_local(fig, plot_paths.summary_png, plot_paths.resolution);

    save_full_domain_plot(x_full, y_full, lsf, material_mask_core, dx, dy, ...
        plot_paths.full_domain_png, plot_paths.resolution);
    save_boundary_compare_plot(x_full, y_full, lsf, material_mask_core, material_mask_full, ...
        dx, dy, diag_info, plot_paths.boundary_compare_png, plot_paths.resolution);
    save_raw_quality_plot(raw_history, raw_metrics, plot_paths.raw_quality_png, plot_paths.resolution);
end

function [history, metrics] = resolve_raw_quality_info(diag_info)
    history = [];
    metrics = [];

    if isfield(diag_info, 'raw_history') && ~isempty(diag_info.raw_history)
        history = diag_info.raw_history;
    elseif isfield(diag_info, 'path_quality_history') && ~isempty(diag_info.path_quality_history)
        history = diag_info.path_quality_history;
    end

    if isfield(diag_info, 'path_quality_raw') && ~isempty(diag_info.path_quality_raw)
        metrics = diag_info.path_quality_raw;
    end
end

function plot_paths = resolve_plot_paths(diag_info)
    plot_paths = struct();
    plot_paths.resolution = 150;
    if isfield(diag_info, 'plot_resolution') && ~isempty(diag_info.plot_resolution)
        plot_paths.resolution = diag_info.plot_resolution;
    end

    base_dir = '';
    if isfield(diag_info, 'plot_dir') && ~isempty(diag_info.plot_dir)
        base_dir = diag_info.plot_dir;
    elseif isfield(diag_info, 'out_dir') && ~isempty(diag_info.out_dir)
        base_dir = diag_info.out_dir;
    elseif isfield(diag_info, 'artifact_dir') && ~isempty(diag_info.artifact_dir)
        base_dir = diag_info.artifact_dir;
    elseif isfield(diag_info, 'save_dir') && ~isempty(diag_info.save_dir)
        base_dir = diag_info.save_dir;
    elseif isfield(diag_info, 'runtime_paths') && isstruct(diag_info.runtime_paths) && ...
            isfield(diag_info.runtime_paths, 'visualization_dir') && ...
            ~isempty(diag_info.runtime_paths.visualization_dir)
        base_dir = diag_info.runtime_paths.visualization_dir;
    else
        base_dir = resolve_runtime_paths().visualization_dir;
    end
    if exist(base_dir, 'dir') ~= 7
        mkdir(base_dir);
    end

    plot_paths.summary_png = fullfile(base_dir, 'optimization_summary.png');
    plot_paths.full_domain_png = fullfile(base_dir, 'full_domain_diagnostic.png');
    plot_paths.boundary_compare_png = fullfile(base_dir, 'boundary_vs_zero_contour.png');
    plot_paths.raw_quality_png = fullfile(base_dir, 'raw_path_quality_history.png');

    if isfield(diag_info, 'plot_paths') && isstruct(diag_info.plot_paths)
        user_paths = diag_info.plot_paths;
        if isfield(user_paths, 'summary_png') && ~isempty(user_paths.summary_png)
            plot_paths.summary_png = user_paths.summary_png;
        end
        if isfield(user_paths, 'full_domain_png') && ~isempty(user_paths.full_domain_png)
            plot_paths.full_domain_png = user_paths.full_domain_png;
        end
        if isfield(user_paths, 'boundary_compare_png') && ~isempty(user_paths.boundary_compare_png)
            plot_paths.boundary_compare_png = user_paths.boundary_compare_png;
        end
        if isfield(user_paths, 'raw_quality_png') && ~isempty(user_paths.raw_quality_png)
            plot_paths.raw_quality_png = user_paths.raw_quality_png;
        end
    end
end

function save_full_domain_plot(x_full, y_full, lsf, material_mask_core, dx, dy, out_path, resolution)
    fig = figure('Visible', 'off', 'Position', [140, 140, 980, 700]);
    contour(x_full, y_full, lsf, 30, 'LineWidth', 0.8, 'Color', [0.35, 0.35, 0.60]);
    hold on;
    contour(x_full, y_full, lsf, [0, 0], 'r', 'LineWidth', 2);
    plot_material_boundary(material_mask_core, dx, dy, '--', [0.10, 0.10, 0.10], 1.2);
    title('全域未裁剪 LSF 诊断');
    xlabel('x (m)');
    ylabel('y (m)');
    axis equal tight;
    grid on;
    save_png_figure_local(fig, out_path, resolution);
    close(fig);
end

function save_boundary_compare_plot(~, ~, lsf, material_mask_core, material_mask_full, ...
    dx, dy, diag_info, out_path, resolution)
    fig = figure('Visible', 'off', 'Position', [160, 160, 980, 700]);
    hold on;

    if isfield(diag_info, 'init_boundary_geometry') && ~isempty(diag_info.init_boundary_geometry)
        geom = diag_info.init_boundary_geometry;
        if isfield(geom, 'pixel_boundary') && ~isempty(geom.pixel_boundary)
            plot(geom.pixel_boundary(:, 1), geom.pixel_boundary(:, 2), '.', ...
                'Color', [0.35, 0.35, 0.35], 'MarkerSize', 8, 'DisplayName', 'pixel boundary');
        end
        boundary_segments = collect_boundary_segments(geom);
        for k = 1:numel(boundary_segments)
            seg = boundary_segments{k};
            plot(seg(:, 1), seg(:, 2), '-', 'Color', [0.10, 0.45, 0.20], ...
                'LineWidth', 1.2, 'DisplayName', ternary_label(k == 1, 'subcell boundary', ''));
        end
    else
        plot_material_boundary(material_mask_core, dx, dy, '-', [0.25, 0.25, 0.25], 1.0);
    end

    zero_geom = extract_levelset_geometry(lsf, dx, dy, 0, material_mask_full);
    for k = 1:numel(zero_geom.segments)
        seg = zero_geom.segments{k};
        plot(seg(:, 1), seg(:, 2), 'r-', 'LineWidth', 1.6, ...
            'DisplayName', ternary_label(k == 1, 'zero contour', ''));
    end
    axis equal tight;
    xlabel('x (m)');
    ylabel('y (m)');
    title('像素边界 / 重建边界 / 零等值线');
    legend('Location', 'best');
    grid on;
    save_png_figure_local(fig, out_path, resolution);
    close(fig);
end

function save_raw_quality_plot(raw_history, raw_metrics, out_path, resolution)
    fig = figure('Visible', 'off', 'Position', [180, 180, 980, 700]);
    if isempty(raw_history)
        axis off;
        text(0.5, 0.55, 'No raw path-quality history', 'HorizontalAlignment', 'center');
        if ~isempty(raw_metrics)
            text(0.5, 0.40, sprintf('mean|turn|=%.3f deg, max|kappa|=%.3e', ...
                raw_metrics.mean_abs_turn_deg, raw_metrics.max_abs_kappa), ...
                'HorizontalAlignment', 'center');
        end
    else
        yyaxis left;
        plot(raw_history.mean_abs_turn_deg, 'LineWidth', 1.6, 'DisplayName', 'mean |turn| (deg)');
        hold on;
        if isfield(raw_history, 'max_abs_kappa')
            plot(raw_history.max_abs_kappa, 'LineWidth', 1.3, 'DisplayName', 'max |kappa|');
        end
        ylabel('turn / curvature');

        yyaxis right;
        plot(raw_history.parallel_spacing_error_percent, 'LineWidth', 1.3, 'DisplayName', 'spacing err (%)');
        hold on;
        plot(raw_history.grad_dev_mean, 'LineWidth', 1.3, 'DisplayName', 'mean ||grad|-1|');
        ylabel('spacing / grad dev');

        xlabel('Iteration');
        title('Raw 路径质量历史');
        grid on;
        legend('Location', 'best');
    end
    save_png_figure_local(fig, out_path, resolution);
    close(fig);
end

function save_png_figure_local(fig, out_path, resolution)
    [parent_dir, ~, ~] = fileparts(out_path);
    if ~isempty(parent_dir) && exist(parent_dir, 'dir') ~= 7
        mkdir(parent_dir);
    end
    try
        exportgraphics(fig, out_path, 'Resolution', resolution);
    catch exportErr
        warning('visualize_results_article:exportgraphicsFailed', ...
            'exportgraphics failed for %s (%s); falling back to print.', ...
            out_path, exportErr.message);
        print(fig, out_path, '-dpng', sprintf('-r%d', resolution));
    end
end

function mask_full = expand_material_mask_to_full_local(mask_core)
    mask_full = false(size(mask_core, 1) + 2, size(mask_core, 2) + 2);
    mask_full(2:end-1, 2:end-1) = logical(mask_core);
    mask_full(1, :) = mask_full(2, :);
    mask_full(end, :) = mask_full(end-1, :);
    mask_full(:, 1) = mask_full(:, 2);
    mask_full(:, end) = mask_full(:, end-1);
end

function plot_material_boundary(material_mask_core, dx, dy, line_style, color_value, line_width)
    hold on;
    geom = reconstruct_material_boundary_subpixel(material_mask_core, dx, dy);
    if ~isfield(geom, 'segments') || isempty(geom.segments)
        return;
    end
    for k = 1:numel(geom.segments)
        seg = geom.segments{k};
        if isempty(seg)
            continue;
        end
        plot(seg(:, 1), seg(:, 2), 'LineStyle', line_style, ...
            'Color', color_value, 'LineWidth', line_width);
    end
end

function segments = collect_boundary_segments(geom)
    segments = {};
    if isfield(geom, 'segments') && ~isempty(geom.segments)
        segments = geom.segments;
    elseif isfield(geom, 'components') && ~isempty(geom.components)
        segments = geom.components;
    end
    if isempty(segments)
        return;
    end
    valid = false(size(segments));
    for k = 1:numel(segments)
        seg = segments{k};
        valid(k) = ~isempty(seg) && size(seg, 1) >= 2 && size(seg, 2) >= 2;
    end
    segments = segments(valid);
end

function label = ternary_label(condition, true_label, false_label)
    if condition
        label = true_label;
    else
        label = false_label;
    end
end
