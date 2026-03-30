function init_info = compute_boundary_offset_stats(lsf, material_mask, dx, dy, ...
    delta_phi_target, delta_phi_used, max_inner_distance, boundary_geometry, phi_boundary_core)
    % 统计边界偏移初始化的各项诊断信息

    mask = logical(material_mask);
    if nargin < 8 || isempty(boundary_geometry)
        boundary_geometry = reconstruct_material_boundary_subpixel(mask, dx, dy);
    end
    if nargin < 9 || isempty(phi_boundary_core)
        sign_reference = ones(size(mask));
        sign_reference(mask) = -1;
        [phi_boundary_core, dist_info] = build_signed_distance_from_segments( ...
            boundary_geometry.x_centers, boundary_geometry.y_centers, ...
            boundary_geometry.segments, sign_reference, true(size(mask)), 0);
        if ~dist_info.success
            error('无法根据重建边界恢复符号距离场。');
        end
    end

    contour_info = extract_subcell_isocontours(lsf, 0, dx, dy, expand_material_mask_to_full(mask));
    [contour_x, contour_y] = flatten_segments(contour_info.segments);

    distances = contour_to_boundary_distance(contour_info.segments, boundary_geometry.segments);
    if isempty(distances)
        mean_offset = 0;
        std_offset = 0;
        mean_error = 0;
        max_error = 0;
        distance_samples = [];
        sample_indices = [];
    else
        mean_offset = mean(distances);
        std_offset = std(distances);
        mean_error = mean(distances - delta_phi_used);
        max_error = max(abs(distances - delta_phi_used));
        samples_limit = 1000;
        if numel(distances) > samples_limit
            sample_indices = round(linspace(1, numel(distances), samples_limit));
        else
            sample_indices = 1:numel(distances);
        end
        distance_samples = distances(sample_indices);
    end

    inward_distance = -phi_boundary_core;
    thin_mask = mask & (inward_distance < delta_phi_used);
    total_cells = max(1, nnz(mask));
    thin_ratio = nnz(thin_mask) / total_cells;

    init_info = struct();
    init_info.mean_offset = mean_offset;
    init_info.std_offset = std_offset;
    init_info.mean_error = mean_error;
    init_info.max_error = max_error;
    init_info.num_samples = numel(distances);
    init_info.distance_samples = distance_samples;
    init_info.distance_sample_indices = sample_indices;
    init_info.distance_sample_errors = distance_samples - delta_phi_used;
    init_info.contour = struct('x', contour_x(:)', 'y', contour_y(:)');
    init_info.boundary_coords = struct( ...
        'x', boundary_geometry.pixel_boundary(:, 1)', ...
        'y', boundary_geometry.pixel_boundary(:, 2)');
    init_info.thin_ratio = thin_ratio;
    init_info.num_thin_cells = nnz(thin_mask);
    init_info.total_material_cells = total_cells;
    if total_cells > 0
        init_info.min_thickness = 2 * min(inward_distance(mask));
    else
        init_info.min_thickness = 0;
    end
    init_info.max_inner_distance = max_inner_distance;
    init_info.delta_phi_target = delta_phi_target;
    init_info.delta_phi_used = delta_phi_used;
    init_info.boundary_method = boundary_geometry.method;
    init_info.boundary_geometry = boundary_geometry;
    init_info.phi_boundary_core = phi_boundary_core;
    init_info.phi_boundary_full = expand_phi_to_full_grid(phi_boundary_core);
    init_info.grad_phi_near_zero = compute_zero_band_gradient_stats(lsf, dx, dy, mask);
end

function [x_points, y_points] = flatten_segments(segments)
    x_points = [];
    y_points = [];
    for k = 1:numel(segments)
        seg = segments{k};
        x_points = [x_points, seg(:, 1)']; %#ok<AGROW>
        y_points = [y_points, seg(:, 2)']; %#ok<AGROW>
        if k < numel(segments)
            x_points(end+1) = NaN; %#ok<AGROW>
            y_points(end+1) = NaN; %#ok<AGROW>
        end
    end
end

function distances = contour_to_boundary_distance(contour_segments, boundary_segments)
    distances = [];
    if isempty(contour_segments) || isempty(boundary_segments)
        return;
    end

    for k = 1:numel(contour_segments)
        pts = contour_segments{k};
        if size(pts, 1) < 2
            continue;
        end
        seg_dist = inf(size(pts, 1), 1);
        for j = 1:numel(boundary_segments)
            target = boundary_segments{j};
            if size(target, 1) < 2
                continue;
            end
            seg_dist = min(seg_dist, point_to_polyline_distance(pts, target));
        end
        distances = [distances; seg_dist(isfinite(seg_dist))]; %#ok<AGROW>
    end
end

function min_distance = point_to_polyline_distance(points, polyline)
    min_distance = inf(size(points, 1), 1);
    for s = 1:size(polyline, 1) - 1
        a = polyline(s, :);
        b = polyline(s + 1, :);
        ab = b - a;
        denom = dot(ab, ab);
        if denom <= eps
            proj = repmat(a, size(points, 1), 1);
        else
            t = ((points(:, 1) - a(1)) * ab(1) + (points(:, 2) - a(2)) * ab(2)) / denom;
            t = min(max(t, 0), 1);
            proj = a + t .* ab;
        end
        dist = hypot(points(:, 1) - proj(:, 1), points(:, 2) - proj(:, 2));
        min_distance = min(min_distance, dist);
    end
end

function phi_full = expand_phi_to_full_grid(phi_core)
    phi_full = zeros(size(phi_core, 1) + 2, size(phi_core, 2) + 2);
    phi_full(2:end-1, 2:end-1) = phi_core;
    phi_full(1, :) = phi_full(2, :);
    phi_full(end, :) = phi_full(end-1, :);
    phi_full(:, 1) = phi_full(:, 2);
    phi_full(:, end) = phi_full(:, end-1);
end

function grad_stats = compute_zero_band_gradient_stats(lsf, dx, dy, material_mask)
    [grad_y, grad_x] = gradient(lsf, dy, dx);
    grad_mag = hypot(grad_x, grad_y);

    material_mask_full = expand_material_mask_to_full(material_mask);
    band_mask = (abs(lsf) <= min(dx, dy)) & material_mask_full;
    band_values = grad_mag(band_mask);
    if isempty(band_values)
        grad_stats = struct('mean', NaN, 'std', NaN, 'max_abs_deviation', NaN, ...
            'p95_abs_deviation', NaN, 'sample_count', 0);
        return;
    end

    abs_deviation = abs(band_values - 1);
    grad_stats = struct();
    grad_stats.mean = mean(band_values);
    grad_stats.std = std(band_values);
    grad_stats.max_abs_deviation = max(abs_deviation);
    grad_stats.p95_abs_deviation = prctile(abs_deviation, 95);
    grad_stats.sample_count = numel(band_values);
end
