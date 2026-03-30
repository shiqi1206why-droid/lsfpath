function metrics = compute_lsf_path_quality(lsf, dx, dy, material_mask_core, opts)
%COMPUTE_LSF_PATH_QUALITY Raw geometric quality metrics for the current LSF.

    if nargin < 5 || isempty(opts)
        opts = struct();
    end
    h = min(dx, dy);
    if ~isfield(opts, 'resample_ds') || isempty(opts.resample_ds)
        opts.resample_ds = h / 4;
    end
    if ~isfield(opts, 'grad_bandwidth') || isempty(opts.grad_bandwidth)
        opts.grad_bandwidth = 0.75 * h;
    end
    if ~isfield(opts, 'parallel_levels') || isempty(opts.parallel_levels)
        opts.parallel_levels = h * (1:3);
    end
    if ~isfield(opts, 'grad_outlier_bounds') || isempty(opts.grad_outlier_bounds)
        opts.grad_outlier_bounds = [0.5, 1.5];
    end
    if ~isfield(opts, 'boundary_overlap_bandwidth') || isempty(opts.boundary_overlap_bandwidth)
        opts.boundary_overlap_bandwidth = 2.0 * h;
    end

    material_mask_full = normalize_mask_to_lsf_grid_local(material_mask_core, size(lsf));
    interior_mask_full = false(size(lsf));
    interior_mask_full(2:end-1, 2:end-1) = true;
    interior_mask_full = interior_mask_full & material_mask_full;

    zero_contour_info = extract_subcell_isocontours(lsf, 0, dx, dy, material_mask_full);
    zero_geometry = struct();
    zero_geometry.level = 0;
    zero_geometry.segments = zero_contour_info.segments;
    zero_geometry.segment_count = zero_contour_info.segment_count;
    zero_geometry.component_count = zero_contour_info.segment_count;
    zero_geometry.segment_lengths = zero_contour_info.segment_lengths;
    zero_geometry.total_length = sum(zero_contour_info.segment_lengths);
    zero_geometry.success = zero_contour_info.success;

    segment_metrics = repmat(empty_segment_metrics(), 0, 1);
    total_length = 0;
    weighted_turn = 0;
    max_kappa = 0;
    max_turn = 0;

    for idx = 1:numel(zero_geometry.segments)
        seg = zero_geometry.segments{idx};
        if size(seg, 1) < 3
            continue;
        end
        closed = hypot(seg(1, 1) - seg(end, 1), seg(1, 2) - seg(end, 2)) <= 1.5 * opts.resample_ds;
        [xs, ys] = resample_polyline(seg(:, 1), seg(:, 2), opts.resample_ds, closed);
        sm = compute_segment_metrics(xs, ys);
        if sm.length <= 0
            continue;
        end
        total_length = total_length + sm.length;
        weighted_turn = weighted_turn + sm.mean_abs_turn_deg * sm.length;
        max_kappa = max(max_kappa, sm.max_abs_kappa);
        max_turn = max(max_turn, sm.max_abs_turn_deg);
        sm.segment_id = idx;
        segment_metrics(end + 1) = sm; %#ok<AGROW>
    end

    [grad_y, grad_x] = gradient(lsf, dy, dx);
    grad_mag = hypot(grad_x, grad_y);
    near_zero_band = abs(lsf) <= opts.grad_bandwidth;
    boundary_distance = compute_material_boundary_distance(material_mask_full, h);
    boundary_band = boundary_distance <= opts.boundary_overlap_bandwidth;

    grad_scopes = struct();
    grad_scopes.all_material = compute_gradient_scope_stats( ...
        grad_mag, near_zero_band, boundary_band, material_mask_full, opts.grad_outlier_bounds);
    grad_scopes.interior_material_excluding_ghost = compute_gradient_scope_stats( ...
        grad_mag, near_zero_band, boundary_band, interior_mask_full, opts.grad_outlier_bounds);

    primary_scope_name = 'interior_material_excluding_ghost';
    primary_scope = grad_scopes.(primary_scope_name);
    if primary_scope.valid_count == 0
        primary_scope_name = 'all_material';
        primary_scope = grad_scopes.(primary_scope_name);
    end

    spacing_levels = opts.parallel_levels(:)';
    material_phi = abs(lsf(material_mask_full));
    if isempty(material_phi)
        max_inside_distance = 0;
    else
        max_inside_distance = max(material_phi);
    end
    spacing_levels = spacing_levels(spacing_levels <= max_inside_distance + 1e-12 & spacing_levels > 0);
    spacing_measured = nan(size(spacing_levels));
    spacing_rel_error = nan(size(spacing_levels));
    for idx = 1:numel(spacing_levels)
        lv = spacing_levels(idx);
        spacing_measured(idx) = verify_path_spacing(lsf, 0, lv, dx, dy);
        if isfinite(spacing_measured(idx)) && lv > 0
            spacing_rel_error(idx) = abs(spacing_measured(idx) - lv) / lv;
        end
    end

    metrics = struct();
    metrics.zero_contour = zero_geometry;
    metrics.zero_segment_metrics = segment_metrics;
    metrics.zero_segment_count = numel(segment_metrics);
    metrics.zero_total_length = total_length;
    if total_length > 0
        metrics.zero_mean_abs_turn_deg = weighted_turn / total_length;
    else
        metrics.zero_mean_abs_turn_deg = NaN;
    end
    metrics.zero_max_abs_turn_deg = max_turn;
    metrics.zero_max_curvature = max_kappa;
    if max_kappa > 0
        metrics.zero_min_turn_radius = 1 / max_kappa;
    else
        metrics.zero_min_turn_radius = inf;
    end
    metrics.parallel_levels = spacing_levels;
    metrics.parallel_spacing_measured = spacing_measured;
    metrics.parallel_spacing_relative_error = spacing_rel_error;
    finite_spacing_error = spacing_rel_error(isfinite(spacing_rel_error));
    if isempty(finite_spacing_error)
        metrics.parallel_spacing_mean_error = NaN;
        metrics.parallel_spacing_max_error = NaN;
    else
        metrics.parallel_spacing_mean_error = mean(finite_spacing_error);
        metrics.parallel_spacing_max_error = max(finite_spacing_error);
    end

    metrics.grad_bandwidth = opts.grad_bandwidth;
    metrics.gradient_scopes = grad_scopes;
    metrics.gradient_primary_scope = primary_scope_name;
    metrics.grad_near_zero_count = primary_scope.near_zero_count;
    metrics.grad_abs_deviation_mean = primary_scope.grad_abs_deviation_mean;
    metrics.grad_abs_deviation_p95 = primary_scope.grad_abs_deviation_p95;
    metrics.grad_abs_deviation_max = primary_scope.grad_abs_deviation_max;
    metrics.grad_outlier_ratio_0p5_1p5 = primary_scope.outlier_ratio_0p5_1p5;
    metrics.near_zero_grad_outlier_ratio = primary_scope.near_zero_outlier_ratio;
    metrics.high_grad_outlier_ratio = primary_scope.high_outlier_ratio;
    metrics.low_grad_outlier_ratio = primary_scope.low_outlier_ratio;
    metrics.high_grad_boundary_overlap_ratio = primary_scope.high_grad_boundary_overlap_ratio;
    metrics.near_zero_grad_median = primary_scope.near_zero_grad_median;
    metrics.near_zero_grad_p95 = primary_scope.near_zero_grad_p95;
    metrics.grad_outlier_ratio_0p5_1p5_all_material = grad_scopes.all_material.outlier_ratio_0p5_1p5;
    metrics.near_zero_grad_outlier_ratio_all_material = grad_scopes.all_material.near_zero_outlier_ratio;
    metrics.high_grad_boundary_overlap_ratio_all_material = grad_scopes.all_material.high_grad_boundary_overlap_ratio;
end

function material_mask_full = normalize_mask_to_lsf_grid_local(material_mask, lsf_size)
    if isequal(size(material_mask), lsf_size)
        material_mask_full = logical(material_mask);
        return;
    end

    core_size = [lsf_size(1) - 2, lsf_size(2) - 2];
    if all(size(material_mask) == core_size)
        material_mask_full = false(lsf_size);
        material_mask_full(2:end-1, 2:end-1) = logical(material_mask);
        material_mask_full(1, :) = material_mask_full(2, :);
        material_mask_full(end, :) = material_mask_full(end-1, :);
        material_mask_full(:, 1) = material_mask_full(:, 2);
        material_mask_full(:, end) = material_mask_full(:, end-1);
        return;
    end

    error('compute_lsf_path_quality: material_mask size mismatch.');
end

function boundary_distance = compute_material_boundary_distance(material_mask_full, h)
    core_mask = material_mask_full(2:end-1, 2:end-1);
    core_distance = inf(size(core_mask));
    if any(core_mask(:))
        distance_map = bwdist(~core_mask) * h;
        core_distance(core_mask) = distance_map(core_mask);
    end
    boundary_distance = expand_core_field_to_full(core_distance, size(material_mask_full));
    boundary_distance(~material_mask_full) = inf;
end

function field_full = expand_core_field_to_full(field_core, lsf_size)
    core_size = [lsf_size(1) - 2, lsf_size(2) - 2];
    if ~all(size(field_core) == core_size)
        error('compute_lsf_path_quality: core field size mismatch.');
    end
    field_full = inf(lsf_size);
    field_full(2:end-1, 2:end-1) = field_core;
    field_full(1, :) = field_full(2, :);
    field_full(end, :) = field_full(end-1, :);
    field_full(:, 1) = field_full(:, 2);
    field_full(:, end) = field_full(:, end-1);
end

function scope = compute_gradient_scope_stats(grad_mag, near_zero_band, boundary_band, scope_mask, outlier_bounds)
    valid_mask = scope_mask & isfinite(grad_mag);
    near_zero_mask = valid_mask & near_zero_band;
    high_mask = valid_mask & grad_mag > outlier_bounds(2);
    low_mask = valid_mask & grad_mag < outlier_bounds(1);
    outlier_mask = high_mask | low_mask;
    near_zero_outlier_mask = near_zero_mask & (grad_mag > outlier_bounds(2) | grad_mag < outlier_bounds(1));

    near_zero_grad = grad_mag(near_zero_mask);
    near_zero_grad = near_zero_grad(isfinite(near_zero_grad));
    near_zero_dev = abs(near_zero_grad - 1);

    scope = struct();
    scope.valid_count = nnz(valid_mask);
    scope.near_zero_count = nnz(near_zero_mask);
    scope.high_count = nnz(high_mask);
    scope.low_count = nnz(low_mask);
    scope.outlier_count = nnz(outlier_mask);
    scope.high_boundary_overlap_count = nnz(high_mask & boundary_band);
    scope.outlier_ratio_0p5_1p5 = safe_ratio(scope.outlier_count, scope.valid_count);
    scope.near_zero_outlier_ratio = safe_ratio(nnz(near_zero_outlier_mask), scope.near_zero_count);
    scope.high_outlier_ratio = safe_ratio(scope.high_count, scope.valid_count);
    scope.low_outlier_ratio = safe_ratio(scope.low_count, scope.valid_count);
    scope.high_grad_boundary_overlap_ratio = safe_ratio(scope.high_boundary_overlap_count, scope.high_count);
    scope.near_zero_grad_median = median_or_nan(near_zero_grad);
    scope.near_zero_grad_p95 = prctile_or_nan(near_zero_grad, 95);
    scope.grad_abs_deviation_mean = mean_or_nan(near_zero_dev);
    scope.grad_abs_deviation_p95 = prctile_or_nan(near_zero_dev, 95);
    scope.grad_abs_deviation_max = max_or_nan(near_zero_dev);
end

function sm = empty_segment_metrics()
    sm = struct('segment_id', NaN, 'length', 0, 'mean_abs_turn_deg', 0, ...
        'max_abs_turn_deg', 0, 'max_abs_kappa', 0, 'min_turn_radius', inf);
end

function sm = compute_segment_metrics(x, y)
    x = x(:);
    y = y(:);
    keep = [true; hypot(diff(x), diff(y)) > 1e-12];
    x = x(keep);
    y = y(keep);

    sm = empty_segment_metrics();
    if numel(x) < 4
        return;
    end

    seg_len = hypot(diff(x), diff(y));
    theta = atan2(diff(y), diff(x));
    dtheta = wrap_to_pi_local(diff(theta));
    ds_mid = max((seg_len(1:end-1) + seg_len(2:end)) / 2, 1e-12);
    kappa = abs(dtheta) ./ ds_mid;

    sm.length = sum(seg_len);
    sm.mean_abs_turn_deg = mean(abs(dtheta)) * 180 / pi;
    sm.max_abs_turn_deg = max(abs(dtheta)) * 180 / pi;
    sm.max_abs_kappa = max(kappa);
    if sm.max_abs_kappa > 0
        sm.min_turn_radius = 1 / sm.max_abs_kappa;
    end
end

function [xo, yo] = resample_polyline(x, y, target_ds, closed)
    x = x(:);
    y = y(:);
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
        xw = [x; x(1)];
        yw = [y; y(1)];
    else
        xw = x;
        yw = y;
    end

    seg_len = hypot(diff(xw), diff(yw));
    keep_seg = [true; seg_len > 1e-12];
    xw = xw(keep_seg);
    yw = yw(keep_seg);
    if numel(xw) < 2
        xo = xw;
        yo = yw;
        return;
    end

    seg_len = hypot(diff(xw), diff(yw));
    s = [0; cumsum(seg_len)];
    total_len = s(end);
    if total_len <= 0
        xo = xw;
        yo = yw;
        return;
    end

    n_new = max(20, ceil(total_len / max(target_ds, 1e-12)));
    if closed
        s_query = linspace(0, total_len, n_new + 1)';
        s_query(end) = [];
    else
        s_query = linspace(0, total_len, n_new)';
    end
    xo = interp1(s, xw, s_query, 'linear');
    yo = interp1(s, yw, s_query, 'linear');
end

function value = safe_ratio(num, den)
    if den <= 0
        value = NaN;
    else
        value = num / den;
    end
end

function val = mean_or_nan(data)
    data = data(isfinite(data));
    if isempty(data)
        val = NaN;
    else
        val = mean(data);
    end
end

function val = median_or_nan(data)
    data = data(isfinite(data));
    if isempty(data)
        val = NaN;
    else
        val = median(data);
    end
end

function val = max_or_nan(data)
    data = data(isfinite(data));
    if isempty(data)
        val = NaN;
    else
        val = max(data);
    end
end

function val = prctile_or_nan(data, p)
    data = data(isfinite(data));
    if isempty(data)
        val = NaN;
    else
        val = prctile(data, p);
    end
end

function wrapped = wrap_to_pi_local(angle_value)
    wrapped = mod(angle_value + pi, 2 * pi) - pi;
end
