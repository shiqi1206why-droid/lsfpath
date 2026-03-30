function [lsf_new, diagnostics] = build_signed_distance_from_geometry(lsf_old, geometry, dx, dy, material_mask)
    % 基于子像素折线几何重建有符号距离场

    material_mask_full = normalize_mask_to_lsf_grid_local(material_mask, size(lsf_old));
    sign_field = sign(lsf_old);
    sign_field(sign_field == 0) = 1;
    sign_field(~material_mask_full) = 1;

    diagnostics = struct('method', 'subcell_signed_distance', ...
        'fallback_used', false, 'component_count', geometry.component_count, ...
        'total_length', geometry.total_length);

    if geometry.component_count == 0 || isempty(geometry.segment_start)
        error('build_signed_distance_from_geometry: 几何折线为空。');
    end

    [x_full, y_full] = get_lsf_grid_coordinates(size(lsf_old), dx, dy);
    [X, Y] = meshgrid(x_full, y_full);

    T = inf(size(lsf_old));
    query_mask = material_mask_full;
    query_points_x = X(query_mask);
    query_points_y = Y(query_mask);
    T(query_mask) = compute_point_to_segments_distance_local(query_points_x, query_points_y, ...
        geometry.segment_start, geometry.segment_end);

    finite_inside = T(query_mask & isfinite(T));
    if isempty(finite_inside)
        Tmax = 0;
    else
        Tmax = max(finite_inside);
    end
    pad = max(10 * max(dx, dy), min(dx, dy));
    outside_distance = Tmax + pad;
    if outside_distance <= 0
        outside_distance = pad;
    end

    T(~isfinite(T) & material_mask_full) = outside_distance;
    T(~material_mask_full) = outside_distance;

    lsf_new = sign_field .* T;
    lsf_new(~material_mask_full) = outside_distance;
    lsf_new(1, :) = lsf_new(2, :);
    lsf_new(end, :) = lsf_new(end-1, :);
    lsf_new(:, 1) = lsf_new(:, 2);
    lsf_new(:, end) = lsf_new(:, end-1);

    diagnostics.outside_distance = outside_distance;
    diagnostics.min_inside_distance = min(T(material_mask_full));
    diagnostics.max_inside_distance = max(T(material_mask_full));
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

    error('build_signed_distance_from_geometry: material_mask尺寸不匹配。');
end

function min_distance = compute_point_to_segments_distance_local(xq, yq, seg_start, seg_end)
    px = xq(:);
    py = yq(:);
    point_count = numel(px);
    segment_count = size(seg_start, 1);
    min_distance_sq = inf(point_count, 1);

    sx0 = seg_start(:, 1).';
    sy0 = seg_start(:, 2).';
    sx1 = seg_end(:, 1).';
    sy1 = seg_end(:, 2).';
    vx = sx1 - sx0;
    vy = sy1 - sy0;
    vv = max(vx.^2 + vy.^2, eps);
    batch_size = max(1, floor(2e6 / max(segment_count, 1)));

    for start_idx = 1:batch_size:point_count
        end_idx = min(start_idx + batch_size - 1, point_count);
        px_batch = px(start_idx:end_idx);
        py_batch = py(start_idx:end_idx);
        wx = px_batch - sx0;
        wy = py_batch - sy0;
        t = (wx .* vx + wy .* vy) ./ vv;
        t = min(max(t, 0), 1);
        proj_x = sx0 + t .* vx;
        proj_y = sy0 + t .* vy;
        dist_sq = (px_batch - proj_x).^2 + (py_batch - proj_y).^2;
        min_distance_sq(start_idx:end_idx) = min(dist_sq, [], 2);
    end

    min_distance = sqrt(max(min_distance_sq, 0));
end
