function summary = export_printable_paths_from_lsf(lsf, dx, dy, levels, out_dir, material_mask_core, opts)
% 将水平集等值线导出为打印友好路径（平滑后）
%
% 输入:
%   lsf    - 水平集函数 (含ghost cells)
%   dx,dy  - 网格间距
%   levels - 等值线水平数组，例如 [-2h -h 0 h 2h]
%   out_dir- 输出目录
%   material_mask_core - 材料域掩膜（核心网格）
%   opts   - 参数结构体
%
% 输出:
%   summary - 汇总统计结构体

    if nargin < 6 || (nargin == 6 && isstruct(material_mask_core))
        opts = material_mask_core;
        material_mask_core = [];
    end
    if nargin < 7 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'target_ds') || isempty(opts.target_ds)
        opts.target_ds = min(dx, dy) / 4;
    end
    if ~isfield(opts, 'smooth_window') || isempty(opts.smooth_window)
        opts.smooth_window = 9;
    end
    if ~isfield(opts, 'smooth_iters') || isempty(opts.smooth_iters)
        opts.smooth_iters = 3;
    end
    if ~isfield(opts, 'min_points') || isempty(opts.min_points)
        opts.min_points = 20;
    end
    if ~isfield(opts, 'auto_optimize') || isempty(opts.auto_optimize)
        opts.auto_optimize = true;
    end
    if ~isfield(opts, 'max_deviation_limit') || isempty(opts.max_deviation_limit)
        opts.max_deviation_limit = 0.75 * min(dx, dy);
    end
    if ~isfield(opts, 'candidate_windows') || isempty(opts.candidate_windows)
        opts.candidate_windows = [7, 9, 11, 13];
    end
    if ~isfield(opts, 'candidate_iters') || isempty(opts.candidate_iters)
        opts.candidate_iters = [2, 3, 4];
    end
    if ~isfield(opts, 'candidate_methods') || isempty(opts.candidate_methods)
        opts.candidate_methods = {'moving_average', 'chaikin'};
    end
    if ~isfield(opts, 'chaikin_iters') || isempty(opts.chaikin_iters)
        opts.chaikin_iters = [1, 2, 3];
    end
    if mod(opts.smooth_window, 2) == 0
        opts.smooth_window = opts.smooth_window + 1;
    end
    if ~isfield(opts, 'raw_path_quality') || isempty(opts.raw_path_quality)
        opts.raw_path_quality = compute_raw_path_quality_metrics(lsf, dx, dy, material_mask_core);
    end
    if ~isfield(opts, 'raw_primary_only') || isempty(opts.raw_primary_only)
        opts.raw_primary_only = true;
    end
    if ~isfield(opts, 'init_boundary_geometry')
        opts.init_boundary_geometry = [];
    end

    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    [ny, nx] = size(lsf);
    if isempty(material_mask_core)
        material_mask_core = true(ny - 2, nx - 2);
    elseif ~isequal(size(material_mask_core), [ny - 2, nx - 2])
        error('material_mask_core尺寸错误：期望[%d,%d]，实际[%d,%d]。', ...
            ny - 2, nx - 2, size(material_mask_core, 1), size(material_mask_core, 2));
    else
        material_mask_core = logical(material_mask_core);
    end

    [x_coords, y_coords] = get_lsf_grid_coordinates(size(lsf), dx, dy);

    fig = figure('Visible', 'off', 'Position', [120, 120, 1300, 820]);
    tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile(1);
    hold on; axis equal; grid on;
    title('Raw Paths');
    xlabel('x (m)'); ylabel('y (m)');

    nexttile(2);
    hold on; axis equal; grid on;
    title('Smoothed Paths (auxiliary)');
    xlabel('x (m)'); ylabel('y (m)');

    metrics_rows = {};
    seg_global_id = 0;
    sampled_point_count = 0;
    sampled_point_violations = 0;

    contour_entries = extract_path_contours(lsf, levels, x_coords, y_coords, ...
        material_mask_core, dx, dy, opts);
    for ci = 1:numel(contour_entries)
        lv = contour_entries{ci}.level;
        raw_xy = contour_entries{ci}.xy;
        if size(raw_xy, 1) < opts.min_points
            continue;
        end

        seg_global_id = seg_global_id + 1;
        smooth_out = smooth_path_candidates(raw_xy, material_mask_core, dx, dy, opts);
        xr = smooth_out.raw_x;
        yr = smooth_out.raw_y;
        xs = smooth_out.smooth_x;
        ys = smooth_out.smooth_y;
        mraw = smooth_out.raw_metrics;
        msmooth = smooth_out.smooth_metrics;
        meta = smooth_out.meta;

        inside_smooth = is_points_in_material(xs, ys, material_mask_core, dx, dy);
        sampled_point_count = sampled_point_count + numel(inside_smooth);
        sampled_point_violations = sampled_point_violations + nnz(~inside_smooth);

        % 可视化
        nexttile(1);
        plot(xr, yr, '-', 'LineWidth', 1.1);
        nexttile(2);
        plot(xs, ys, '-', 'LineWidth', 1.5);

        % 导出CSV
        csv_name = sprintf('path_level_%+.4f_seg_%03d.csv', lv, seg_global_id);
        csv_path = fullfile(out_dir, csv_name);
        write_path_to_file(csv_path, xr, yr, xs, ys);

        metrics_rows(end+1, :) = { ... %#ok<AGROW>
            seg_global_id, lv, numel(xr), numel(xs), ...
            mraw.length, msmooth.length, ...
            mraw.mean_abs_turn_deg, msmooth.mean_abs_turn_deg, ...
            mraw.max_abs_turn_deg, msmooth.max_abs_turn_deg, ...
            mraw.max_abs_kappa, msmooth.max_abs_kappa, ...
            mraw.min_turn_radius, msmooth.min_turn_radius, ...
            msmooth.max_deviation_from_raw, ...
            string(meta.method), meta.param1, meta.param2};
    end

    % 汇总表
    metrics_header = { ...
        'seg_id','level','raw_n','smooth_n', ...
        'raw_len','smooth_len', ...
        'raw_mean_abs_turn_deg','smooth_mean_abs_turn_deg', ...
        'raw_max_abs_turn_deg','smooth_max_abs_turn_deg', ...
        'raw_max_abs_kappa','smooth_max_abs_kappa', ...
        'raw_min_turn_radius','smooth_min_turn_radius', ...
        'smooth_max_deviation_from_raw', ...
        'smooth_method','smooth_param1','smooth_param2'};
    if isempty(metrics_rows)
        warning('export_printable_paths_from_lsf: 未提取到可用路径段。');
        metrics_rows = cell(0, numel(metrics_header));
    end
    metrics_table = cell2table(metrics_rows, 'VariableNames', metrics_header);

    metrics_csv = fullfile(out_dir, 'path_smoothing_metrics.csv');
    writetable(metrics_table, metrics_csv);

    % 汇总统计
    summary = struct();
    summary.segment_count = height(metrics_table);
    summary.levels = levels(:)';
    summary.target_ds = opts.target_ds;
    summary.smooth_window = opts.smooth_window;
    summary.smooth_iters = opts.smooth_iters;
    summary.auto_optimize = logical(opts.auto_optimize);
    summary.max_deviation_limit = opts.max_deviation_limit;
    if height(metrics_table) == 0
        summary.raw_mean_abs_turn_deg = NaN;
        summary.smooth_mean_abs_turn_deg = NaN;
        summary.raw_max_abs_turn_deg = NaN;
        summary.smooth_max_abs_turn_deg = NaN;
        summary.raw_max_abs_kappa = NaN;
        summary.smooth_max_abs_kappa = NaN;
        summary.raw_min_turn_radius = NaN;
        summary.smooth_min_turn_radius = NaN;
        summary.max_deviation_from_raw = NaN;
        summary.selected_methods = strings(0, 1);
    else
        summary.raw_mean_abs_turn_deg = mean(metrics_table.raw_mean_abs_turn_deg);
        summary.smooth_mean_abs_turn_deg = mean(metrics_table.smooth_mean_abs_turn_deg);
        summary.raw_max_abs_turn_deg = max(metrics_table.raw_max_abs_turn_deg);
        summary.smooth_max_abs_turn_deg = max(metrics_table.smooth_max_abs_turn_deg);
        summary.raw_max_abs_kappa = max(metrics_table.raw_max_abs_kappa);
        summary.smooth_max_abs_kappa = max(metrics_table.smooth_max_abs_kappa);
        summary.raw_min_turn_radius = min(metrics_table.raw_min_turn_radius);
        summary.smooth_min_turn_radius = min(metrics_table.smooth_min_turn_radius);
        summary.max_deviation_from_raw = max(metrics_table.smooth_max_deviation_from_raw);
        summary.selected_methods = unique(metrics_table.smooth_method)';
    end
    summary.metrics_csv = metrics_csv;
    summary.material_mask_size = size(material_mask_core);
    summary.sampled_point_count = sampled_point_count;
    summary.sampled_point_violations = sampled_point_violations;
    summary.raw_path_quality = opts.raw_path_quality;
    summary.raw_primary_only = logical(opts.raw_primary_only);
    summary.primary_metric_source = 'raw_lsf_contours';
    summary.smoothing_role = 'auxiliary_printability_only';
    if ~isempty(opts.init_boundary_geometry)
        summary.init_boundary_geometry = opts.init_boundary_geometry;
    end

    fig_path = fullfile(out_dir, 'printable_paths_overlay.png');
    save_png_figure(fig, fig_path, 140);
    close(fig);
    summary.overlay_png = fig_path;

    summary_txt = fullfile(out_dir, 'path_smoothing_summary.txt');
    summary.summary_txt = summary_txt;
    write_summary_text(summary_txt, summary);
end

function segments = parse_contourc_segments(C)
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
        pts = C(:, j0:j1)';
        segments{end+1} = pts; %#ok<AGROW>
        idx = j1 + 1;
    end
end

function save_png_figure(fig, fig_path, resolution)
    if nargin < 3 || isempty(resolution)
        resolution = 140;
    end

    try
        exportgraphics(fig, fig_path, 'Resolution', resolution);
    catch exportErr
        warning('export_printable_paths_from_lsf:exportgraphicsFailed', ...
            'exportgraphics failed for %s (%s); falling back to print.', ...
            fig_path, exportErr.message);
        print(fig, fig_path, '-dpng', sprintf('-r%d', resolution));
    end
end

function clipped_segments = clip_segment_to_material(raw_xy, material_mask_core, dx, dy, opts)
    xr = raw_xy(:, 1);
    yr = raw_xy(:, 2);
    closed = hypot(xr(1) - xr(end), yr(1) - yr(end)) <= 1.5 * opts.target_ds;
    clip_ds = min(opts.target_ds, min(dx, dy) / 4);
    [xc, yc] = resample_polyline(xr, yr, clip_ds, closed);
    inside = is_points_in_material(xc, yc, material_mask_core, dx, dy);

    clipped_segments = {};
    start_idx = [];
    for k = 1:numel(inside)
        if inside(k) && isempty(start_idx)
            start_idx = k;
        elseif ~inside(k) && ~isempty(start_idx)
            seg = compact_segment(xc(start_idx:k-1), yc(start_idx:k-1));
            if size(seg, 1) >= 2
                clipped_segments{end+1} = seg; %#ok<AGROW>
            end
            start_idx = [];
        end
    end
    if ~isempty(start_idx)
        seg = compact_segment(xc(start_idx:end), yc(start_idx:end));
        if size(seg, 1) >= 2
            clipped_segments{end+1} = seg; %#ok<AGROW>
        end
    end
end

function inside = is_points_in_material(x, y, material_mask_core, dx, dy)
    [nely, nelx] = size(material_mask_core);
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
    inside(finite_mask) = material_mask_core(linear_idx);
end

function seg = compact_segment(x, y)
    seg = [x(:), y(:)];
    keep = [true; hypot(diff(seg(:, 1)), diff(seg(:, 2))) > 1e-12];
    seg = seg(keep, :);
end

function [xs, ys, mraw, msmooth, meta] = smooth_segment_for_print(xr, yr, material_mask_core, dx, dy, opts)
    xr = xr(:);
    yr = yr(:);

    % 去重复点
    keep = [true; hypot(diff(xr), diff(yr)) > 1e-12];
    xr = xr(keep);
    yr = yr(keep);

    if numel(xr) < 5
        xs = xr;
        ys = yr;
        mraw = polyline_metrics(xr, yr);
        msmooth = mraw;
        msmooth.max_deviation_from_raw = 0;
        meta = struct('method', 'identity', 'param1', 0, 'param2', 0);
        return;
    end

    closed = hypot(xr(1)-xr(end), yr(1)-yr(end)) <= 1.5 * opts.target_ds;
    if closed
        xr = xr(1:end-1);
        yr = yr(1:end-1);
    end

    [x_res, y_res] = resample_polyline(xr, yr, opts.target_ds, closed);
    if ~is_segment_inside_material(x_res, y_res, material_mask_core, dx, dy)
        x_res = xr;
        y_res = yr;
    end
    mraw = polyline_metrics(x_res, y_res);
    [xs, ys, msmooth, meta] = choose_best_smoothed_segment( ...
        x_res, y_res, mraw, material_mask_core, dx, dy, opts, closed);
end

function [xo, yo] = resample_polyline(x, y, target_ds, closed)
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

function [xs, ys] = smooth_polyline(x, y, win, niter, closed)
    xs = x(:); ys = y(:);
    for k = 1:niter
        if closed
            xs = periodic_moving_average(xs, win);
            ys = periodic_moving_average(ys, win);
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

function [xs_best, ys_best, m_best, meta_best] = choose_best_smoothed_segment(x_res, y_res, mraw, material_mask_core, dx, dy, opts, closed)
    [x_default, y_default] = smooth_polyline(x_res, y_res, opts.smooth_window, opts.smooth_iters, closed);
    [xs_default, ys_default] = resample_polyline(x_default, y_default, opts.target_ds, closed);
    if is_segment_inside_material(xs_default, ys_default, material_mask_core, dx, dy)
        m_default = polyline_metrics(xs_default, ys_default);
        m_default.max_deviation_from_raw = max_min_distance(xs_default, ys_default, x_res, y_res);
        xs_best = xs_default;
        ys_best = ys_default;
        m_best = m_default;
        meta_best = struct('method', 'moving_average', 'param1', opts.smooth_window, 'param2', opts.smooth_iters);
    else
        xs_best = x_res;
        ys_best = y_res;
        m_best = mraw;
        m_best.max_deviation_from_raw = 0;
        meta_best = struct('method', 'identity', 'param1', 0, 'param2', 0);
    end

    if ~opts.auto_optimize
        return;
    end

    candidates = {};
    for method_idx = 1:numel(opts.candidate_methods)
        method = opts.candidate_methods{method_idx};
        switch lower(method)
            case 'moving_average'
                for win = opts.candidate_windows
                    win_use = win + mod(win + 1, 2);
                    for niter = opts.candidate_iters
                        [x_sm, y_sm] = smooth_polyline(x_res, y_res, win_use, niter, closed);
                        [xs, ys] = resample_polyline(x_sm, y_sm, opts.target_ds, closed);
                        if ~is_segment_inside_material(xs, ys, material_mask_core, dx, dy)
                            continue;
                        end
                        m = polyline_metrics(xs, ys);
                        m.max_deviation_from_raw = max_min_distance(xs, ys, x_res, y_res);
                        candidates{end+1} = {xs, ys, m, struct('method', 'moving_average', 'param1', win_use, 'param2', niter)}; %#ok<AGROW>
                    end
                end
            case 'chaikin'
                for niter = opts.chaikin_iters
                    [x_sm, y_sm] = chaikin_smooth_polyline(x_res, y_res, niter, closed);
                    [xs, ys] = resample_polyline(x_sm, y_sm, opts.target_ds, closed);
                    if ~is_segment_inside_material(xs, ys, material_mask_core, dx, dy)
                        continue;
                    end
                    m = polyline_metrics(xs, ys);
                    m.max_deviation_from_raw = max_min_distance(xs, ys, x_res, y_res);
                    candidates{end+1} = {xs, ys, m, struct('method', 'chaikin', 'param1', niter, 'param2', 0)}; %#ok<AGROW>
                end
        end
    end

    for i = 1:numel(candidates)
        cand = candidates{i};
        xs = cand{1};
        ys = cand{2};
        m = cand{3};
        meta = cand{4};
        if m.max_deviation_from_raw > opts.max_deviation_limit
            continue;
        end
        if is_better_segment(m, meta, m_best, meta_best)
            xs_best = xs;
            ys_best = ys;
            m_best = m;
            meta_best = meta;
        end
    end
end

function tf = is_segment_inside_material(x, y, material_mask_core, dx, dy)
    inside = is_points_in_material(x, y, material_mask_core, dx, dy);
    tf = all(inside);
end

function tf = is_better_segment(m_new, meta_new, m_old, meta_old)
    score_new = [m_new.max_abs_kappa, m_new.mean_abs_turn_deg, m_new.max_deviation_from_raw];
    score_old = [m_old.max_abs_kappa, m_old.mean_abs_turn_deg, m_old.max_deviation_from_raw];

    if any(score_new < score_old) && ~any(score_new > score_old)
        tf = true;
        return;
    end
    if any(score_old < score_new) && ~any(score_old > score_new)
        tf = false;
        return;
    end

    tf = false;
    if score_new(1) < score_old(1) * (1 - 1e-6)
        tf = true;
    elseif abs(score_new(1) - score_old(1)) <= max(1e-9, 1e-6 * score_old(1))
        if score_new(2) < score_old(2) * (1 - 1e-6)
            tf = true;
        elseif abs(score_new(2) - score_old(2)) <= max(1e-9, 1e-6 * score_old(2))
            if score_new(3) < score_old(3) * (1 - 1e-6)
                tf = true;
            elseif abs(score_new(3) - score_old(3)) <= max(1e-9, 1e-6 * score_old(3))
                tf = prefer_method(meta_new, meta_old);
            end
        end
    end
end

function tf = prefer_method(meta_new, meta_old)
    order = struct('chaikin', 1, 'moving_average', 2, 'identity', 3);
    key_new = lower(meta_new.method);
    key_old = lower(meta_old.method);
    if isfield(order, key_new) && isfield(order, key_old)
        tf = order.(key_new) < order.(key_old);
    else
        tf = false;
    end
end

function [xs, ys] = chaikin_smooth_polyline(x, y, niter, closed)
    xs = x(:);
    ys = y(:);
    for k = 1:niter
        if closed
            xw = [xs; xs(1)];
            yw = [ys; ys(1)];
            nseg = numel(xs);
            x_new = zeros(2 * nseg, 1);
            y_new = zeros(2 * nseg, 1);
            idx = 1;
            for i = 1:nseg
                p0 = [xw(i), yw(i)];
                p1 = [xw(i+1), yw(i+1)];
                q = 0.75 * p0 + 0.25 * p1;
                r = 0.25 * p0 + 0.75 * p1;
                x_new(idx:idx+1) = [q(1); r(1)];
                y_new(idx:idx+1) = [q(2); r(2)];
                idx = idx + 2;
            end
            xs = x_new;
            ys = y_new;
        else
            n = numel(xs);
            x_new = zeros(2 * n, 1);
            y_new = zeros(2 * n, 1);
            x_new(1) = xs(1);
            y_new(1) = ys(1);
            pos = 2;
            for i = 1:(n-1)
                p0 = [xs(i), ys(i)];
                p1 = [xs(i+1), ys(i+1)];
                q = 0.75 * p0 + 0.25 * p1;
                r = 0.25 * p0 + 0.75 * p1;
                x_new(pos:pos+1) = [q(1); r(1)];
                y_new(pos:pos+1) = [q(2); r(2)];
                pos = pos + 2;
            end
            x_new(pos) = xs(end);
            y_new(pos) = ys(end);
            xs = x_new(1:pos);
            ys = y_new(1:pos);
        end
    end
end

function ys = periodic_moving_average(x, win)
    x = x(:);
    n = numel(x);
    hw = floor(win/2);
    if n <= win + 2
        ys = x;
        return;
    end
    xpad = [x(end-hw+1:end); x; x(1:hw)];
    kernel = ones(win, 1) / win;
    ypad = conv(xpad, kernel, 'same');
    ys = ypad(hw+1:hw+n);
end

function m = polyline_metrics(x, y)
    x = x(:); y = y(:);
    n = numel(x);
    if n < 4
        m = struct('length', 0, 'mean_abs_turn_deg', 0, 'max_abs_turn_deg', 0, ...
            'max_abs_kappa', 0, 'min_turn_radius', inf);
        return;
    end

    seg = hypot(diff(x), diff(y));
    L = sum(seg);
    theta = atan2(diff(y), diff(x));
    dtheta = wrap_to_pi_local(diff(theta));
    ds_mid = max((seg(1:end-1) + seg(2:end)) / 2, 1e-9);
    kappa = abs(dtheta) ./ ds_mid;

    m = struct();
    m.length = L;
    m.mean_abs_turn_deg = mean(abs(dtheta)) * 180 / pi;
    m.max_abs_turn_deg = max(abs(dtheta)) * 180 / pi;
    m.max_abs_kappa = max(kappa);
    if m.max_abs_kappa > 0
        m.min_turn_radius = 1 / m.max_abs_kappa;
    else
        m.min_turn_radius = inf;
    end
end

function w = wrap_to_pi_local(a)
    w = mod(a + pi, 2*pi) - pi;
end

function dmax = max_min_distance(xa, ya, xb, yb)
    % 计算A点集到B点集的最大最近邻距离（简单O(NM)，当前规模足够）
    xa = xa(:); ya = ya(:);
    xb = xb(:); yb = yb(:);
    dmax = 0;
    for i = 1:numel(xa)
        d2 = (xb - xa(i)).^2 + (yb - ya(i)).^2;
        dmin = sqrt(min(d2));
        if dmin > dmax
            dmax = dmin;
        end
    end
end

function write_path_csv(csv_path, xr, yr, xs, ys)
    n = max(numel(xr), numel(xs));
    xr_pad = nan(n, 1); yr_pad = nan(n, 1);
    xs_pad = nan(n, 1); ys_pad = nan(n, 1);
    xr_pad(1:numel(xr)) = xr;
    yr_pad(1:numel(yr)) = yr;
    xs_pad(1:numel(xs)) = xs;
    ys_pad(1:numel(ys)) = ys;
    T = table(xr_pad, yr_pad, xs_pad, ys_pad, ...
        'VariableNames', {'x_raw','y_raw','x_smooth','y_smooth'});
    writetable(T, csv_path);
end

function write_summary_text(path_str, s)
    fid = fopen(path_str, 'w');
    if fid < 0
        error('无法写入摘要文件: %s', path_str);
    end
    c = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, 'Printable Path Smoothing Summary\n');
    fprintf(fid, 'segment_count: %d\n', s.segment_count);
    fprintf(fid, 'target_ds: %.6e\n', s.target_ds);
    fprintf(fid, 'smooth_window: %d\n', s.smooth_window);
    fprintf(fid, 'smooth_iters: %d\n', s.smooth_iters);
    fprintf(fid, 'auto_optimize: %d\n', s.auto_optimize);
    fprintf(fid, 'max_deviation_limit: %.6e\n', s.max_deviation_limit);
    fprintf(fid, 'raw_mean_abs_turn_deg: %.6f\n', s.raw_mean_abs_turn_deg);
    fprintf(fid, 'smooth_mean_abs_turn_deg: %.6f\n', s.smooth_mean_abs_turn_deg);
    fprintf(fid, 'raw_max_abs_turn_deg: %.6f\n', s.raw_max_abs_turn_deg);
    fprintf(fid, 'smooth_max_abs_turn_deg: %.6f\n', s.smooth_max_abs_turn_deg);
    fprintf(fid, 'raw_max_abs_kappa: %.6e\n', s.raw_max_abs_kappa);
    fprintf(fid, 'smooth_max_abs_kappa: %.6e\n', s.smooth_max_abs_kappa);
    fprintf(fid, 'raw_min_turn_radius: %.6e\n', s.raw_min_turn_radius);
    fprintf(fid, 'smooth_min_turn_radius: %.6e\n', s.smooth_min_turn_radius);
    fprintf(fid, 'max_deviation_from_raw: %.6e\n', s.max_deviation_from_raw);
    fprintf(fid, 'sampled_point_count: %d\n', s.sampled_point_count);
    fprintf(fid, 'sampled_point_violations: %d\n', s.sampled_point_violations);
    if isfield(s, 'primary_metric_source')
        fprintf(fid, 'primary_metric_source: %s\n', s.primary_metric_source);
    end
    if isfield(s, 'smoothing_role')
        fprintf(fid, 'smoothing_role: %s\n', s.smoothing_role);
    end
    if isfield(s, 'raw_primary_only')
        fprintf(fid, 'raw_primary_only: %d\n', s.raw_primary_only);
    end
    if isfield(s, 'raw_path_quality') && ~isempty(s.raw_path_quality)
        fprintf(fid, 'raw_zero_mean_abs_turn_deg: %.6f\n', s.raw_path_quality.mean_abs_turn_deg);
        fprintf(fid, 'raw_zero_max_abs_kappa: %.6e\n', s.raw_path_quality.max_abs_kappa);
        fprintf(fid, 'raw_parallel_spacing_error_percent: %.6f\n', s.raw_path_quality.parallel_spacing_error_percent);
        fprintf(fid, 'raw_grad_dev_mean: %.6e\n', s.raw_path_quality.grad_dev_mean);
    end
    fprintf(fid, 'selected_methods: %s\n', strjoin(cellstr(s.selected_methods), ', '));
    fprintf(fid, 'metrics_csv: %s\n', s.metrics_csv);
    fprintf(fid, 'overlay_png: %s\n', s.overlay_png);
end
