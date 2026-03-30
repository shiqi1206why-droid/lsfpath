function metrics = capture_baseline_validation(config_name)
%CAPTURE_BASELINE_VALIDATION Run a config and save minimal validation artifacts.

    if nargin < 1 || isempty(config_name)
        config_name = 'fast';
    end

    root_dir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(fullfile(root_dir, 'utilities'), '-begin');
    project_root = get_project_root(root_dir);
    cleanup_path = ensure_project_on_path(project_root); %#ok<NASGU>
    paths = build_project_paths(project_root);
    set(0, 'DefaultFigureVisible', 'off');

    output_root = paths.baseline_dir;
    if ~exist(output_root, 'dir')
        mkdir(output_root);
    end

    run_stamp = datestr(now, 'yyyymmdd_HHMMSS');
    run_dir = fullfile(output_root, sprintf('%s_run_%s', lower(config_name), run_stamp));
    if ~exist(run_dir, 'dir')
        mkdir(run_dir);
    end

    timer_handle = tic;
    results = fiber_levelset(config_name);
    runtime_seconds = toc(timer_handle);
    params = results.params;

    save(fullfile(run_dir, 'results.mat'), 'results', '-v7.3');

    [x_full, y_full] = get_lsf_grid_coordinates(size(results.lsf), params.grid.dx, params.grid.dy);

    fig1 = figure('Visible', 'off', 'Position', [80, 80, 1100, 650]);
    yyaxis left;
    plot(results.compliance_history, 'LineWidth', 1.6);
    ylabel('Compliance');
    yyaxis right;
    plot(results.FCS_history * 100, 'LineWidth', 1.2);
    ylabel('FCS (%)');
    xlabel('Iteration');
    grid on;
    title(sprintf('%s convergence', config_name));
    metrics_convergence_png = fullfile(run_dir, 'convergence.png');
    save_png(fig1, metrics_convergence_png, 130);
    close(fig1);

    fig2 = figure('Visible', 'off', 'Position', [90, 90, 900, 650]);
    contour(x_full, y_full, results.lsf, 28, 'LineWidth', 0.7);
    hold on;
    contour(x_full, y_full, results.lsf, [0 0], 'r', 'LineWidth', 2);
    draw_material_boundary(results.material_mask_core, params.grid.dx, params.grid.dy);
    axis equal tight;
    xlabel('x (m)');
    ylabel('y (m)');
    title(sprintf('%s full-domain diagnostic', config_name));
    metrics_full_domain_png = fullfile(run_dir, 'full_domain_diagnostic.png');
    save_png(fig2, metrics_full_domain_png, 130);
    close(fig2);

    fig3 = figure('Visible', 'off', 'Position', [100, 100, 900, 650]);
    lsf_masked = results.lsf;
    lsf_masked(~results.material_mask_full) = NaN;
    contour(x_full, y_full, lsf_masked, 28, 'LineWidth', 0.7);
    hold on;
    contour(x_full, y_full, lsf_masked, [0 0], 'r', 'LineWidth', 2);
    axis equal tight;
    xlabel('x (m)');
    ylabel('y (m)');
    title(sprintf('%s masked LSF', config_name));
    metrics_final_png = fullfile(run_dir, 'final_lsf.png');
    save_png(fig3, metrics_final_png, 130);
    close(fig3);

    fig4 = figure('Visible', 'off', 'Position', [120, 120, 1000, 650]);
    if isfield(results, 'path_quality_history') && ~isempty(results.path_quality_history)
        hst = results.path_quality_history;
        yyaxis left;
        plot(hst.mean_abs_turn_deg, 'LineWidth', 1.5);
        hold on;
        plot(hst.max_abs_kappa, 'LineWidth', 1.2);
        ylabel('turn / kappa');
        yyaxis right;
        plot(hst.parallel_spacing_error_percent, 'LineWidth', 1.2);
        hold on;
        plot(hst.grad_dev_mean, 'LineWidth', 1.2);
        ylabel('spacing err / grad dev');
        legend('mean|turn| (deg)', 'max|kappa|', 'spacing err (%)', 'mean||grad|-1|', ...
            'Location', 'best');
    else
        text(0.5, 0.5, 'No raw path quality history', 'HorizontalAlignment', 'center');
        axis off;
    end
    xlabel('Iteration');
    grid on;
    title(sprintf('%s raw path quality history', config_name));
    metrics_raw_png = fullfile(run_dir, 'raw_path_quality_history.png');
    save_png(fig4, metrics_raw_png, 130);
    close(fig4);

    metrics = struct();
    metrics.run_stamp = run_stamp;
    metrics.config = config_name;
    metrics.runtime_seconds = runtime_seconds;
    metrics.final_compliance = results.final_compliance;
    metrics.final_FCS = results.final_FCS;
    metrics.best_compliance = results.best_compliance;
    metrics.best_iter = results.best_iter;
    metrics.path_quality_raw = results.path_quality_raw;
    metrics.interface_diagnostics = struct( ...
        'outside_phi_nonpositive_count', results.interface_diagnostics.outside_phi_nonpositive_count, ...
        'outside_velocity_nonzero_count', results.interface_diagnostics.outside_velocity_nonzero_count, ...
        'last_reinit_method', get_last_reinit_method(results));
    metrics.artifacts = struct( ...
        'run_dir', run_dir, ...
        'results_mat', fullfile(run_dir, 'results.mat'), ...
        'convergence_png', metrics_convergence_png, ...
        'full_domain_png', metrics_full_domain_png, ...
        'final_lsf_png', metrics_final_png, ...
        'raw_quality_png', metrics_raw_png);

    metrics_json = fullfile(run_dir, 'validation_metrics.json');
    fid = fopen(metrics_json, 'w');
    if fid < 0
        error('无法写入 validation_metrics.json');
    end
    cleanup_fid = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, '%s', jsonencode(metrics));
end

function save_png(fig, fig_path, resolution)
    try
        exportgraphics(fig, fig_path, 'Resolution', resolution);
    catch
        print(fig, fig_path, '-dpng', sprintf('-r%d', resolution));
    end
end

function draw_material_boundary(material_mask_core, dx, dy)
    boundary_mask = bwperim(material_mask_core);
    [by, bx] = find(boundary_mask);
    plot((bx - 0.5) * dx, (by - 0.5) * dy, 'k--', 'LineWidth', 1.1);
end

function method_name = get_last_reinit_method(results)
    method_name = '';
    if isfield(results, 'interface_diagnostics') && ...
            isfield(results.interface_diagnostics, 'last_reinit_info') && ...
            isstruct(results.interface_diagnostics.last_reinit_info) && ...
            isfield(results.interface_diagnostics.last_reinit_info, 'method_used')
        method_name = results.interface_diagnostics.last_reinit_info.method_used;
    end
end
