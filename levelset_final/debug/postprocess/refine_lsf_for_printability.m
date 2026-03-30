function refinement = refine_lsf_for_printability(base_lsf, params, out_dir, material_mask, opts)
% 在最终优化后的lsf上做制造友好refinement，并用FE重新校验柔度。
%
% 目标：
% 1) 优先降低原始等值线路径的折转与曲率
% 2) 最终导出的打印路径仍保持平滑
% 3) 对每个候选重新计算theta/FE/compliance/FCS，避免只看图像

    if nargin < 4 || (nargin == 4 && isstruct(material_mask))
        opts = material_mask;
        material_mask = [];
    end
    if nargin < 5 || isempty(opts)
        opts = struct();
    end

    opts = apply_defaults(opts, params);
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    if isempty(material_mask)
        material_mask = load_material_mask_from_topology(params);
    end
    export_opts = struct( ...
        'target_ds', opts.target_ds, ...
        'smooth_window', opts.path_smooth_window, ...
        'smooth_iters', opts.path_smooth_iters, ...
        'min_points', opts.min_points, ...
        'raw_primary_only', true);
    if isfield(opts, 'init_boundary_geometry') && ~isempty(opts.init_boundary_geometry)
        export_opts.init_boundary_geometry = opts.init_boundary_geometry;
    end

    baseline_dir = fullfile(out_dir, 'baseline');
    base_state = evaluate_refinement_state(base_lsf, params, material_mask, opts, true);
    export_opts.raw_path_quality = compute_raw_path_quality_metrics(base_lsf, params.grid.dx, params.grid.dy, material_mask);
    base_summary = export_printable_paths_from_lsf(base_lsf, params.grid.dx, params.grid.dy, ...
        opts.levels, baseline_dir, material_mask, export_opts);
    base_entry = build_candidate_entry('baseline', NaN, NaN, NaN, base_lsf, ...
        base_state, base_summary, opts.initial_compliance);
    base_entry.valid = true;
    base_entry.selected = false;

    candidates = repmat(base_entry, 0, 1);
    candidates(end+1) = base_entry; %#ok<AGROW>

    cand_root = fullfile(out_dir, 'candidates');
    if ~exist(cand_root, 'dir')
        mkdir(cand_root);
    end

    cand_id = 0;
    for win = opts.zero_smooth_windows
        for niter = opts.zero_smooth_iters
            cand_id = cand_id + 1;
            cand_name = sprintf('cand_%02d_zero_w%d_i%d', cand_id, win, niter);
            cand_dir = fullfile(cand_root, cand_name);

            lsf_candidate = build_candidate_from_smoothed_zero_contours(base_lsf, ...
                params.grid.dx, params.grid.dy, opts.zero_target_ds, win, niter, material_mask);
            cand_state = evaluate_refinement_state(lsf_candidate, params, material_mask, opts, false);
            export_opts.raw_path_quality = compute_raw_path_quality_metrics(lsf_candidate, params.grid.dx, params.grid.dy, material_mask);
            cand_summary = export_printable_paths_from_lsf(lsf_candidate, params.grid.dx, params.grid.dy, ...
                opts.levels, cand_dir, material_mask, export_opts);

            entry = build_candidate_entry(cand_name, NaN, NaN, NaN, lsf_candidate, ...
                cand_state, cand_summary, opts.initial_compliance);
            entry.mode = 'zero_contour';
            entry.param1 = win;
            entry.param2 = niter;
            entry.param3 = opts.zero_target_ds;
            entry.compliance_rel_loss = (entry.final_compliance - base_entry.final_compliance) / ...
                max(base_entry.final_compliance, eps);
            entry.valid = isfinite(entry.final_compliance) && ...
                entry.compliance_rel_loss <= opts.max_rel_compliance_loss;
            entry.selected = false;
            candidates(end+1) = entry; %#ok<AGROW>
        end
    end

    for band_factor = opts.bandwidth_factors
        for blend = opts.blend_values
            for iter_count = opts.lsf_smooth_iters
                cand_id = cand_id + 1;
                cand_name = sprintf('cand_%02d_b%.2f_l%.2f_i%d', ...
                    cand_id, band_factor, blend, iter_count);
                cand_dir = fullfile(cand_root, cand_name);

                lsf_candidate = smooth_lsf_in_narrow_band(base_lsf, blend, iter_count, ...
                    band_factor * params.grid.h, material_mask);
                zero_mask = compute_zero_mask_from_lsf(lsf_candidate, params.grid.h);
                zero_mask = zero_mask & expand_mask_to_full(material_mask);
                lsf_candidate = fmm_reinitialize(lsf_candidate, params.grid.dx, params.grid.dy, ...
                    zero_mask, material_mask);

                cand_state = evaluate_refinement_state(lsf_candidate, params, material_mask, opts, false);
                export_opts.raw_path_quality = compute_raw_path_quality_metrics(lsf_candidate, params.grid.dx, params.grid.dy, material_mask);
                cand_summary = export_printable_paths_from_lsf(lsf_candidate, params.grid.dx, params.grid.dy, ...
                    opts.levels, cand_dir, material_mask, export_opts);

                entry = build_candidate_entry(cand_name, band_factor, blend, iter_count, ...
                    lsf_candidate, cand_state, cand_summary, opts.initial_compliance);
                entry.mode = 'band_lsf';
                entry.param1 = band_factor;
                entry.param2 = blend;
                entry.param3 = iter_count;
                entry.compliance_rel_loss = (entry.final_compliance - base_entry.final_compliance) / ...
                    max(base_entry.final_compliance, eps);
                entry.valid = isfinite(entry.final_compliance) && ...
                    entry.compliance_rel_loss <= opts.max_rel_compliance_loss;
                entry.selected = false;
                candidates(end+1) = entry; %#ok<AGROW>
            end
        end
    end

    valid_mask = [candidates.valid];
    valid_candidates = candidates(valid_mask);
    best_idx_in_valid = select_best_candidate(valid_candidates);
    best_candidate = valid_candidates(best_idx_in_valid);
    best_candidate.selected = true;

    selected_idx = find(strcmp({candidates.name}, best_candidate.name), 1, 'first');
    candidates(selected_idx).selected = true;

    compare_fig = fullfile(out_dir, 'refinement_compare.png');
    save_refinement_compare_figure(base_lsf, best_candidate.lsf, params, material_mask, compare_fig, ...
        base_entry.name, best_candidate.name);
    best_lsf_png = fullfile(out_dir, 'manufacturing_lsf.png');
    save_single_lsf_figure(best_candidate.lsf, params, material_mask, best_lsf_png, ...
        sprintf('Manufacturing LSF (%s)', best_candidate.name));

    table_path = fullfile(out_dir, 'refinement_candidates.csv');
    write_candidate_table(table_path, candidates);

    summary_path = fullfile(out_dir, 'refinement_summary.txt');
    write_refinement_summary(summary_path, base_entry, best_candidate, opts, table_path, compare_fig);

    refinement = struct();
    refinement.base = strip_large_fields(base_entry);
    refinement.best = strip_large_fields(best_candidate);
    refinement.max_rel_compliance_loss = opts.max_rel_compliance_loss;
    refinement.levels = opts.levels(:)';
    refinement.candidate_csv = table_path;
    refinement.summary_txt = summary_path;
    refinement.compare_png = compare_fig;
    refinement.best_lsf_png = best_lsf_png;
    refinement.best_overlay_png = best_candidate.overlay_png;
    refinement.best_metrics_csv = best_candidate.metrics_csv;
    refinement.best_name = best_candidate.name;
    refinement.baseline_name = base_entry.name;
    refinement.refined_selected = ~strcmp(best_candidate.name, base_entry.name);
    refinement.candidate_count = numel(candidates);
    refinement.path_smooth_window = opts.path_smooth_window;
    refinement.path_smooth_iters = opts.path_smooth_iters;
    refinement.primary_selection_basis = ...
        'raw_max_abs_kappa_then_raw_mean_abs_turn_deg_then_final_compliance';
    refinement.smoothing_role = 'auxiliary_printability_only';
end

function opts = apply_defaults(opts, params)
    if ~isfield(opts, 'levels') || isempty(opts.levels)
        opts.levels = (-2:2) * params.grid.h;
    end
    if ~isfield(opts, 'target_ds') || isempty(opts.target_ds)
        opts.target_ds = params.grid.h / 4;
    end
    if ~isfield(opts, 'path_smooth_window') || isempty(opts.path_smooth_window)
        opts.path_smooth_window = 9;
    end
    if ~isfield(opts, 'path_smooth_iters') || isempty(opts.path_smooth_iters)
        opts.path_smooth_iters = 3;
    end
    if ~isfield(opts, 'min_points') || isempty(opts.min_points)
        opts.min_points = 20;
    end
    if ~isfield(opts, 'blend_values') || isempty(opts.blend_values)
        opts.blend_values = [0.05, 0.10, 0.15];
    end
    if ~isfield(opts, 'lsf_smooth_iters') || isempty(opts.lsf_smooth_iters)
        opts.lsf_smooth_iters = [1, 2];
    end
    if ~isfield(opts, 'bandwidth_factors') || isempty(opts.bandwidth_factors)
        opts.bandwidth_factors = [2.0, 3.0];
    end
    if ~isfield(opts, 'max_rel_compliance_loss') || isempty(opts.max_rel_compliance_loss)
        opts.max_rel_compliance_loss = 0.005;
    end
    if ~isfield(opts, 'initial_compliance') || isempty(opts.initial_compliance)
        opts.initial_compliance = NaN;
    end
    if ~isfield(opts, 'theta_adjust_limit_deg') || isempty(opts.theta_adjust_limit_deg)
        opts.theta_adjust_limit_deg = 0.20;
    end
    if ~isfield(opts, 'zero_target_ds') || isempty(opts.zero_target_ds)
        opts.zero_target_ds = params.grid.h / 4;
    end
    if ~isfield(opts, 'zero_smooth_windows') || isempty(opts.zero_smooth_windows)
        opts.zero_smooth_windows = [9, 11, 13];
    end
    if ~isfield(opts, 'zero_smooth_iters') || isempty(opts.zero_smooth_iters)
        opts.zero_smooth_iters = [2, 3, 4];
    end
    opts.theta_adjust_limit = opts.theta_adjust_limit_deg * pi / 180;
end

function material_mask = load_material_mask_from_topology(params)
    paths = resolve_runtime_paths(params);
    topo = load(paths.topology_file);
    if ~isfield(topo, 'struc')
        error('topo_result.mat 缺少 struc 字段。');
    end

    struc = topo.struc;
    if size(struc, 1) ~= params.grid.nely || size(struc, 2) ~= params.grid.nelx
        struc = imresize(struc, [params.grid.nely, params.grid.nelx], 'nearest');
    end
    [material_mask, ~] = clean_material_mask(struc, ...
        params.init.min_component_size, params.init.morph_radius);
end

function state = evaluate_refinement_state(lsf, params, material_mask, opts, use_baseline_state)
    if nargin >= 5 && use_baseline_state && ...
            isfield(opts, 'baseline_theta') && ~isempty(opts.baseline_theta) && ...
            isfield(opts, 'baseline_compliance') && ~isempty(opts.baseline_compliance)
        state = evaluate_state_with_theta(lsf, opts.baseline_theta, ...
            params.grid.nelx, params.grid.nely, material_mask, ...
            params.material.E_L, params.material.E_T, params.material.nu_LT, ...
            params.material.G_LT, params.material.thickness, params.load.F_mag, ...
            params.grid.dx, params.grid.dy);
        state.compliance = opts.baseline_compliance;
        if isfield(opts, 'baseline_FCS') && ~isempty(opts.baseline_FCS)
            state.FCS = opts.baseline_FCS;
        end
        state.theta_target = opts.baseline_theta;
        return;
    end

    if isfield(opts, 'theta_reference') && ~isempty(opts.theta_reference)
        state = evaluate_candidate_state(lsf, opts.theta_reference, opts.theta_adjust_limit, ...
            params.grid.dx, params.grid.dy, params.grid.nelx, params.grid.nely, ...
            material_mask, params.material.E_L, params.material.E_T, params.material.nu_LT, ...
            params.material.G_LT, params.material.thickness, params.load.F_mag, ...
            params.smooth.eta, params.smooth.iterations);
    else
        [theta_direct, theta_target] = advance_theta_state(lsf, [], params.opt.delta_theta_max, ...
            params.grid.dx, params.grid.dy, material_mask, params.smooth.eta, params.smooth.iterations);
        state = evaluate_state_with_theta(lsf, theta_direct, params.grid.nelx, params.grid.nely, ...
            material_mask, params.material.E_L, params.material.E_T, params.material.nu_LT, ...
            params.material.G_LT, params.material.thickness, params.load.F_mag, ...
            params.grid.dx, params.grid.dy);
        state.theta_target = theta_target;
    end
end

function lsf_new = smooth_lsf_in_narrow_band(lsf, blend, num_iters, bandwidth, material_mask)
    lsf_new = lsf;
    kernel = [0, 1, 0; 1, 0, 1; 0, 1, 0] / 4;
    active_mask = abs(lsf_new) <= bandwidth;
    if nargin >= 5 && ~isempty(material_mask)
        active_mask = active_mask & expand_mask_to_full(material_mask);
    end

    for k = 1:num_iters
        avg_field = conv2(lsf_new, kernel, 'same');
        lsf_new(active_mask) = (1 - blend) * lsf_new(active_mask) + blend * avg_field(active_mask);
        lsf_new(1, :) = lsf_new(2, :);
        lsf_new(end, :) = lsf_new(end-1, :);
        lsf_new(:, 1) = lsf_new(:, 2);
        lsf_new(:, end) = lsf_new(:, end-1);
    end
end

function lsf_candidate = build_candidate_from_smoothed_zero_contours(base_lsf, dx, dy, target_ds, win, niter, material_mask)
    [ny, nx] = size(base_lsf);
    x_coords = linspace(0, dx * (nx - 2), nx);
    y_coords = linspace(0, dy * (ny - 2), ny);

    C = contourc(x_coords, y_coords, base_lsf, [0, 0]);
    segments = parse_contourc_segments_local(C);
    zero_mask = false(size(base_lsf));

    for i = 1:numel(segments)
        pts = segments{i};
        if size(pts, 1) < 5
            continue;
        end
        clipped_segments = clip_segment_to_material_local(pts, material_mask, dx, dy, target_ds);
        for j = 1:numel(clipped_segments)
            clipped_pts = clipped_segments{j};
            if size(clipped_pts, 1) < 5
                continue;
            end
            x = clipped_pts(:, 1);
            y = clipped_pts(:, 2);
            closed = hypot(x(1) - x(end), y(1) - y(end)) <= 1.5 * target_ds;
            if closed
                x = x(1:end-1);
                y = y(1:end-1);
            end
            [x_res, y_res] = resample_polyline_local(x, y, target_ds, closed);
            [x_sm, y_sm] = smooth_polyline_local(x_res, y_res, win, niter, closed);
            [x_fin, y_fin] = resample_polyline_local(x_sm, y_sm, target_ds, closed);
            zero_mask = rasterize_polyline_to_mask(zero_mask, x_fin, y_fin, x_coords, y_coords);
        end
    end

    if any(zero_mask(:))
        zero_mask = bwmorph(zero_mask, 'clean');
        zero_mask = bwmorph(zero_mask, 'thin', Inf);
    else
        zero_mask = compute_zero_mask_from_lsf(base_lsf, min(dx, dy));
    end

    zero_mask = zero_mask & expand_mask_to_full(material_mask);
    lsf_candidate = fmm_reinitialize(base_lsf, dx, dy, zero_mask, material_mask);
end

function clipped_segments = clip_segment_to_material_local(raw_xy, material_mask, dx, dy, target_ds)
    xr = raw_xy(:, 1);
    yr = raw_xy(:, 2);
    closed = hypot(xr(1) - xr(end), yr(1) - yr(end)) <= 1.5 * target_ds;
    clip_ds = min(target_ds, min(dx, dy) / 4);
    [xc, yc] = resample_polyline_local(xr, yr, clip_ds, closed);
    inside = is_points_in_material_local(xc, yc, material_mask, dx, dy);

    clipped_segments = {};
    start_idx = [];
    for k = 1:numel(inside)
        if inside(k) && isempty(start_idx)
            start_idx = k;
        elseif ~inside(k) && ~isempty(start_idx)
            seg = compact_segment_local(xc(start_idx:k-1), yc(start_idx:k-1));
            if size(seg, 1) >= 2
                clipped_segments{end+1} = seg; %#ok<AGROW>
            end
            start_idx = [];
        end
    end

    if ~isempty(start_idx)
        seg = compact_segment_local(xc(start_idx:end), yc(start_idx:end));
        if size(seg, 1) >= 2
            clipped_segments{end+1} = seg; %#ok<AGROW>
        end
    end
end

function inside = is_points_in_material_local(x, y, material_mask, dx, dy)
    [nely, nelx] = size(material_mask);
    x = x(:);
    y = y(:);
    inside = false(size(x));

    finite_mask = isfinite(x) & isfinite(y);
    if ~any(finite_mask)
        return;
    end

    x_f = x(finite_mask);
    y_f = y(finite_mask);
    col = floor(x_f ./ dx) + 1;
    row = floor(y_f ./ dy) + 1;
    col = min(max(col, 1), nelx);
    row = min(max(row, 1), nely);
    linear_idx = sub2ind([nely, nelx], row, col);
    inside(finite_mask) = material_mask(linear_idx);
end

function seg = compact_segment_local(x, y)
    seg = [x(:), y(:)];
    keep = [true; hypot(diff(seg(:, 1)), diff(seg(:, 2))) > 1e-12];
    seg = seg(keep, :);
end

function entry = build_candidate_entry(name, band_factor, blend, iter_count, lsf, state, summary, initial_compliance)
    entry = struct();
    entry.name = name;
    entry.mode = 'baseline';
    entry.band_factor = band_factor;
    entry.blend = blend;
    entry.lsf_smooth_iters = iter_count;
    entry.param1 = NaN;
    entry.param2 = NaN;
    entry.param3 = NaN;
    entry.lsf = lsf;
    entry.final_compliance = state.compliance;
    entry.final_FCS = state.FCS;
    entry.raw_mean_abs_turn_deg = summary.raw_mean_abs_turn_deg;
    entry.raw_max_abs_kappa = summary.raw_max_abs_kappa;
    entry.raw_min_turn_radius = summary.raw_min_turn_radius;
    entry.smooth_mean_abs_turn_deg = summary.smooth_mean_abs_turn_deg;
    entry.smooth_max_abs_kappa = summary.smooth_max_abs_kappa;
    entry.smooth_min_turn_radius = summary.smooth_min_turn_radius;
    entry.max_deviation_from_raw = summary.max_deviation_from_raw;
    entry.segment_count = summary.segment_count;
    entry.overlay_png = summary.overlay_png;
    entry.metrics_csv = summary.metrics_csv;
    entry.summary_txt = summary.summary_txt;
    entry.compliance_rel_loss = NaN;
    entry.valid = false;
    entry.selected = false;
    if isfinite(initial_compliance) && initial_compliance > 0
        entry.improvement_ratio = (initial_compliance - state.compliance) / initial_compliance * 100;
    else
        entry.improvement_ratio = NaN;
    end
end

function idx = select_best_candidate(candidates)
    raw_kappa = [candidates.raw_max_abs_kappa]';
    raw_turn = [candidates.raw_mean_abs_turn_deg]';
    compliance = [candidates.final_compliance]';
    score = [raw_kappa, raw_turn, compliance];
    [~, idx] = sortrows(score, [1, 2, 3]);
    idx = idx(1);
end

function save_refinement_compare_figure(base_lsf, best_lsf, params, material_mask, fig_path, base_name, best_name)
    x = linspace(0, params.grid.Lx, size(base_lsf, 2));
    y = linspace(0, params.grid.Ly, size(base_lsf, 1));
    levels = (-2:2) * params.grid.h;
    mask_full = expand_mask_to_full(material_mask);
    base_plot = base_lsf;
    best_plot = best_lsf;
    base_plot(~mask_full) = NaN;
    best_plot(~mask_full) = NaN;

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [120, 120, 1200, 520]);
    tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    contour(x, y, base_plot, levels, 'LineWidth', 1.2);
    axis equal; grid on;
    xlabel('x (m)'); ylabel('y (m)');
    title(sprintf('Base LSF Contours (%s)', base_name), 'Interpreter', 'none');

    nexttile;
    contour(x, y, best_plot, levels, 'LineWidth', 1.2);
    axis equal; grid on;
    xlabel('x (m)'); ylabel('y (m)');
    title(sprintf('Refined LSF Contours (%s)', best_name), 'Interpreter', 'none');

    save_png_figure(fig, fig_path, 180);
    close(fig);
end

function save_single_lsf_figure(lsf, params, material_mask, fig_path, title_str)
    x = linspace(0, params.grid.Lx, size(lsf, 2));
    y = linspace(0, params.grid.Ly, size(lsf, 1));
    levels = (-2:2) * params.grid.h;
    mask_full = expand_mask_to_full(material_mask);
    lsf_plot = lsf;
    lsf_plot(~mask_full) = NaN;
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [120, 120, 760, 560]);
    cmap = lines(numel(levels));
    hold on;
    for i = 1:numel(levels)
        lw = 1.3;
        if abs(levels(i)) < 1e-12
            lw = 2.4;
        end
        contour(x, y, lsf_plot, [levels(i), levels(i)], 'Color', cmap(i, :), 'LineWidth', lw);
    end
    axis equal tight;
    grid on;
    xlabel('x (m)');
    ylabel('y (m)');
    title(title_str, 'Interpreter', 'none');
    legend(compose('%.3f', levels), 'Location', 'bestoutside');
    save_png_figure(fig, fig_path, 180);
    close(fig);
end

function save_png_figure(fig, fig_path, resolution)
    if nargin < 3 || isempty(resolution)
        resolution = 180;
    end

    try
        exportgraphics(fig, fig_path, 'Resolution', resolution);
    catch exportErr
        warning('refine_lsf_for_printability:exportgraphicsFailed', ...
            'exportgraphics failed for %s (%s); falling back to print.', ...
            fig_path, exportErr.message);
        print(fig, fig_path, '-dpng', sprintf('-r%d', resolution));
    end
end

function mask_full = expand_mask_to_full(mask_core)
    mask_full = false(size(mask_core, 1) + 2, size(mask_core, 2) + 2);
    mask_full(2:end-1, 2:end-1) = logical(mask_core);
    mask_full(1, :) = mask_full(2, :);
    mask_full(end, :) = mask_full(end-1, :);
    mask_full(:, 1) = mask_full(:, 2);
    mask_full(:, end) = mask_full(:, end-1);
end

function write_candidate_table(csv_path, candidates)
    T = table( ...
        string({candidates.name})', ...
        string({candidates.mode})', ...
        [candidates.valid]', ...
        [candidates.selected]', ...
        [candidates.band_factor]', ...
        [candidates.blend]', ...
        [candidates.lsf_smooth_iters]', ...
        [candidates.param1]', ...
        [candidates.param2]', ...
        [candidates.param3]', ...
        [candidates.final_compliance]', ...
        [candidates.improvement_ratio]', ...
        [candidates.compliance_rel_loss]', ...
        [candidates.final_FCS]', ...
        [candidates.raw_mean_abs_turn_deg]', ...
        [candidates.raw_max_abs_kappa]', ...
        [candidates.raw_min_turn_radius]', ...
        [candidates.smooth_mean_abs_turn_deg]', ...
        [candidates.smooth_max_abs_kappa]', ...
        [candidates.smooth_min_turn_radius]', ...
        [candidates.max_deviation_from_raw]', ...
        [candidates.segment_count]', ...
        'VariableNames', { ...
            'name', 'mode', 'valid', 'selected', 'band_factor', 'blend', 'lsf_smooth_iters', ...
            'param1', 'param2', 'param3', ...
            'final_compliance', 'improvement_ratio', 'compliance_rel_loss', 'final_FCS', ...
            'raw_mean_abs_turn_deg', 'raw_max_abs_kappa', 'raw_min_turn_radius', ...
            'smooth_mean_abs_turn_deg', 'smooth_max_abs_kappa', 'smooth_min_turn_radius', ...
            'max_deviation_from_raw', 'segment_count'});
    writetable(T, csv_path);
end

function write_refinement_summary(path_str, base_entry, best_entry, opts, table_path, compare_fig)
    fid = fopen(path_str, 'w');
    if fid < 0
        error('无法写入摘要文件: %s', path_str);
    end
    cleanup_obj = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, 'baseline_name: %s\n', base_entry.name);
    fprintf(fid, 'best_name: %s\n', best_entry.name);
    fprintf(fid, 'best_mode: %s\n', best_entry.mode);
    fprintf(fid, 'primary_selection_basis: raw_max_abs_kappa_then_raw_mean_abs_turn_deg_then_final_compliance\n');
    fprintf(fid, 'smoothing_role: auxiliary_printability_only\n');
    fprintf(fid, 'max_rel_compliance_loss: %.6f\n', opts.max_rel_compliance_loss);
    fprintf(fid, 'base_final_compliance: %.10e\n', base_entry.final_compliance);
    fprintf(fid, 'best_final_compliance: %.10e\n', best_entry.final_compliance);
    fprintf(fid, 'base_improvement_ratio: %.6f\n', base_entry.improvement_ratio);
    fprintf(fid, 'best_improvement_ratio: %.6f\n', best_entry.improvement_ratio);
    fprintf(fid, 'best_compliance_rel_loss: %.6f\n', best_entry.compliance_rel_loss);
    fprintf(fid, 'base_raw_mean_abs_turn_deg: %.6f\n', base_entry.raw_mean_abs_turn_deg);
    fprintf(fid, 'best_raw_mean_abs_turn_deg: %.6f\n', best_entry.raw_mean_abs_turn_deg);
    fprintf(fid, 'base_raw_max_abs_kappa: %.6e\n', base_entry.raw_max_abs_kappa);
    fprintf(fid, 'best_raw_max_abs_kappa: %.6e\n', best_entry.raw_max_abs_kappa);
    fprintf(fid, 'base_smooth_mean_abs_turn_deg: %.6f\n', base_entry.smooth_mean_abs_turn_deg);
    fprintf(fid, 'best_smooth_mean_abs_turn_deg: %.6f\n', best_entry.smooth_mean_abs_turn_deg);
    fprintf(fid, 'base_smooth_max_abs_kappa: %.6e\n', base_entry.smooth_max_abs_kappa);
    fprintf(fid, 'best_smooth_max_abs_kappa: %.6e\n', best_entry.smooth_max_abs_kappa);
    fprintf(fid, 'candidate_table: %s\n', table_path);
    fprintf(fid, 'compare_figure: %s\n', compare_fig);
    fprintf(fid, 'best_lsf_png: %s\n', strrep(compare_fig, 'refinement_compare.png', 'manufacturing_lsf.png'));
    fprintf(fid, 'best_overlay_png: %s\n', best_entry.overlay_png);
end

function s = strip_large_fields(entry)
    s = rmfield(entry, {'lsf'});
end

function segments = parse_contourc_segments_local(C)
    segments = {};
    if isempty(C)
        return;
    end
    idx = 1;
    ncol = size(C, 2);
    while idx <= ncol
        if idx + 1 > ncol
            break;
        end
        npts = round(C(2, idx));
        j0 = idx + 1;
        j1 = idx + npts;
        if npts <= 1 || j1 > ncol
            break;
        end
        segments{end+1} = C(:, j0:j1)'; %#ok<AGROW>
        idx = j1 + 1;
    end
end

function [xo, yo] = resample_polyline_local(x, y, target_ds, closed)
    x = x(:); y = y(:);
    keep = [true; hypot(diff(x), diff(y)) > 1e-12];
    x = x(keep);
    y = y(keep);

    if numel(x) < 2
        xo = x;
        yo = y;
        return;
    end

    if closed
        if hypot(x(1) - x(end), y(1) - y(end)) <= 1e-12
            x = x(1:end-1);
            y = y(1:end-1);
        end
        if numel(x) < 2
            xo = x;
            yo = y;
            return;
        end
        xw = [x; x(1)];
        yw = [y; y(1)];
    else
        xw = x;
        yw = y;
    end
    ds = hypot(diff(xw), diff(yw));
    keep_step = [true; ds > 1e-12];
    xw = xw(keep_step);
    yw = yw(keep_step);

    if numel(xw) < 2
        xo = xw;
        yo = yw;
        return;
    end

    ds = hypot(diff(xw), diff(yw));
    s = [0; cumsum(ds)];
    [s, unique_idx] = unique(s, 'stable');
    xw = xw(unique_idx);
    yw = yw(unique_idx);
    total_len = s(end);
    if numel(s) < 2 || total_len <= 0
        xo = xw;
        yo = yw;
        return;
    end
    n_new = max(20, ceil(total_len / max(target_ds, 1e-6)));
    if closed
        s_new = linspace(0, total_len, n_new + 1)';
        s_new(end) = [];
    else
        s_new = linspace(0, total_len, n_new)';
    end
    xo = interp1(s, xw, s_new, 'linear');
    yo = interp1(s, yw, s_new, 'linear');
end

function [xs, ys] = smooth_polyline_local(x, y, win, niter, closed)
    xs = x(:);
    ys = y(:);
    for k = 1:niter
        if closed
            xs = periodic_moving_average_local(xs, win);
            ys = periodic_moving_average_local(ys, win);
        else
            x_prev = xs;
            y_prev = ys;
            xs = movmean(xs, win, 'Endpoints', 'shrink');
            ys = movmean(ys, win, 'Endpoints', 'shrink');
            xs(1) = x_prev(1); xs(end) = x_prev(end);
            ys(1) = y_prev(1); ys(end) = y_prev(end);
        end
    end
end

function ys = periodic_moving_average_local(x, win)
    x = x(:);
    n = numel(x);
    hw = floor(win / 2);
    if n <= win + 2
        ys = x;
        return;
    end
    xpad = [x(end-hw+1:end); x; x(1:hw)];
    kernel = ones(win, 1) / win;
    ypad = conv(xpad, kernel, 'same');
    ys = ypad(hw+1:hw+n);
end

function mask = rasterize_polyline_to_mask(mask, x, y, x_coords, y_coords)
    col = interp1(x_coords, 1:numel(x_coords), x(:), 'nearest', 'extrap');
    row = interp1(y_coords, 1:numel(y_coords), y(:), 'nearest', 'extrap');
    col = max(1, min(numel(x_coords), round(col)));
    row = max(1, min(numel(y_coords), round(row)));
    idx = sub2ind(size(mask), row, col);
    mask(idx) = true;
end
