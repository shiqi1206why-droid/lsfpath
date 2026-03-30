function phi_boundary = compute_signed_distance_to_boundary(material_mask, boundary_geometry, dx, dy)
    % 基于重建边界折线计算单元中心到边界的有符号距离场

    mask = logical(material_mask);
    if ~isfield(boundary_geometry, 'segment_start') || isempty(boundary_geometry.segment_start)
        error('compute_signed_distance_to_boundary: 边界几何为空。');
    end

    [nely, nelx] = size(mask);
    x_centers = ((1:nelx) - 0.5) * dx;
    y_centers = ((1:nely) - 0.5) * dy;
    [X, Y] = meshgrid(x_centers, y_centers);

    distances = compute_point_to_segments_distance(X(:), Y(:), ...
        boundary_geometry.segment_start, boundary_geometry.segment_end);
    distances = reshape(distances, size(mask));

    phi_boundary = distances;
    phi_boundary(mask) = -phi_boundary(mask);
end

function min_distance = compute_point_to_segments_distance(xq, yq, seg_start, seg_end)
    px = xq(:);
    py = yq(:);
    point_count = numel(px);
    segment_count = size(seg_start, 1);

    min_distance_sq = inf(point_count, 1);
    batch_size = max(1, floor(1e6 / max(segment_count, 1)));

    sx0 = seg_start(:, 1).';
    sy0 = seg_start(:, 2).';
    sx1 = seg_end(:, 1).';
    sy1 = seg_end(:, 2).';
    vx = sx1 - sx0;
    vy = sy1 - sy0;
    vv = max(vx.^2 + vy.^2, eps);

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
