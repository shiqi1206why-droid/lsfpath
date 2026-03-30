function distance_field = compute_distance_to_polylines(x_coords, y_coords, segments)
%COMPUTE_DISTANCE_TO_POLYLINES Minimum Euclidean distance from grid points to polylines.

    [X, Y] = meshgrid(x_coords, y_coords);
    points = [X(:), Y(:)];
    min_dist2 = inf(size(points, 1), 1);

    for seg_idx = 1:numel(segments)
        seg = segments{seg_idx};
        if size(seg, 1) < 2
            continue;
        end
        for k = 1:(size(seg, 1) - 1)
            p0 = seg(k, :);
            p1 = seg(k + 1, :);
            edge = p1 - p0;
            edge_len2 = sum(edge .^ 2);
            if edge_len2 <= eps
                diff_vec = points - p0;
                dist2 = sum(diff_vec .^ 2, 2);
            else
                rel_vec = points - p0;
                proj = (rel_vec * edge') / edge_len2;
                proj = min(max(proj, 0), 1);
                closest = p0 + proj .* edge;
                diff_vec = points - closest;
                dist2 = sum(diff_vec .^ 2, 2);
            end
            min_dist2 = min(min_dist2, dist2);
        end
    end

    distance_field = reshape(sqrt(min_dist2), size(X));
end
