function [phi, info] = rebuild_signed_distance_from_segments(segments, x_coords, y_coords, sign_reference, material_mask_full)
%REBUILD_SIGNED_DISTANCE_FROM_SEGMENTS Build a signed distance field from polylines.

    if nargin < 5 || isempty(material_mask_full)
        material_mask_full = true(numel(y_coords), numel(x_coords));
    end

    [X, Y] = meshgrid(x_coords, y_coords);
    unsigned_distance = inf(size(X));
    total_segment_edges = 0;

    for si = 1:numel(segments)
        seg = segments{si};
        if size(seg, 1) < 2
            continue;
        end
        for k = 1:size(seg, 1) - 1
            x1 = seg(k, 1);
            y1 = seg(k, 2);
            x2 = seg(k + 1, 1);
            y2 = seg(k + 1, 2);
            edge_len_sq = (x2 - x1)^2 + (y2 - y1)^2;
            if edge_len_sq <= eps
                continue;
            end
            total_segment_edges = total_segment_edges + 1;
            t = ((X - x1) * (x2 - x1) + (Y - y1) * (y2 - y1)) / edge_len_sq;
            t = max(0, min(1, t));
            proj_x = x1 + t * (x2 - x1);
            proj_y = y1 + t * (y2 - y1);
            dist_seg = hypot(X - proj_x, Y - proj_y);
            unsigned_distance = min(unsigned_distance, dist_seg);
        end
    end

    if all(~isfinite(unsigned_distance(:)))
        error('rebuild_signed_distance_from_segments: 未能从折线构造有效距离场。');
    end

    sign_field = sign(double(sign_reference));
    sign_field(sign_field == 0) = 1;
    sign_field(~material_mask_full) = 1;
    phi = unsigned_distance .* sign_field;

    finite_inside = unsigned_distance(isfinite(unsigned_distance) & material_mask_full);
    if isempty(finite_inside)
        outside_distance = max(diff_or_one(x_coords), diff_or_one(y_coords));
    else
        outside_distance = max(finite_inside) + max(diff_or_one(x_coords), diff_or_one(y_coords));
    end
    phi(~material_mask_full) = outside_distance;

    phi(1, :) = phi(2, :);
    phi(end, :) = phi(end-1, :);
    phi(:, 1) = phi(:, 2);
    phi(:, end) = phi(:, end-1);

    info = struct();
    info.segment_count = numel(segments);
    info.segment_edge_count = total_segment_edges;
    info.outside_distance = outside_distance;
    info.min_unsigned_distance = min(unsigned_distance(isfinite(unsigned_distance)), [], 'omitnan');
    info.max_unsigned_distance = max(unsigned_distance(isfinite(unsigned_distance)), [], 'omitnan');
end

function d = diff_or_one(coords)
    if numel(coords) >= 2
        d = mean(diff(coords));
    else
        d = 1;
    end
end
