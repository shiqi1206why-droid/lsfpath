function [signed_distance, info] = build_signed_distance_from_segments(x_coords, y_coords, ...
    segments, sign_reference, valid_mask, outside_value)
%BUILD_SIGNED_DISTANCE_FROM_SEGMENTS Rebuild a signed distance field from polylines.

    if nargin < 5 || isempty(valid_mask)
        valid_mask = true(size(sign_reference));
    end
    if nargin < 6 || isempty(outside_value) || ~isfinite(outside_value)
        outside_value = max(diff_or_one(x_coords), diff_or_one(y_coords));
    end
    outside_value = abs(outside_value);
    if outside_value <= 0
        outside_value = max(diff_or_one(x_coords), diff_or_one(y_coords));
    end

    if ~isequal(size(sign_reference), size(valid_mask))
        error('sign_reference and valid_mask must have identical sizes.');
    end
    if size(sign_reference, 2) ~= numel(x_coords) || size(sign_reference, 1) ~= numel(y_coords)
        error('Coordinate vectors are inconsistent with sign_reference size.');
    end

    signed_distance = outside_value * ones(size(sign_reference), class(sign_reference));
    cleaned_segments = sanitize_segments(segments);
    info = struct( ...
        'success', false, ...
        'segment_count', numel(segments), ...
        'used_segment_count', numel(cleaned_segments), ...
        'min_distance', NaN, ...
        'max_distance', NaN, ...
        'outside_value', outside_value);

    if isempty(cleaned_segments)
        return;
    end

    valid_mask = logical(valid_mask);
    [X, Y] = meshgrid(x_coords, y_coords);
    points = [X(valid_mask), Y(valid_mask)];

    if isempty(points)
        return;
    end

    min_dist = inf(size(points, 1), 1);
    for k = 1:numel(cleaned_segments)
        seg = cleaned_segments{k};
        dist_k = point_to_polyline_distance(points, seg);
        min_dist = min(min_dist, dist_k);
    end

    if ~any(isfinite(min_dist))
        return;
    end

    sign_field = sign(sign_reference);
    sign_field(sign_field == 0) = 1;
    sign_field(~valid_mask) = 1;

    signed_distance(valid_mask) = min_dist .* sign_field(valid_mask);
    signed_distance(~valid_mask) = outside_value;

    info.success = true;
    info.min_distance = min(min_dist);
    info.max_distance = max(min_dist);
end

function cleaned_segments = sanitize_segments(segments)
    cleaned_segments = {};
    for k = 1:numel(segments)
        seg = segments{k};
        if size(seg, 2) ~= 2
            continue;
        end
        finite_rows = all(isfinite(seg), 2);
        seg = seg(finite_rows, :);
        if size(seg, 1) < 2
            continue;
        end
        keep = [true; hypot(diff(seg(:, 1)), diff(seg(:, 2))) > 1e-12];
        seg = seg(keep, :);
        if size(seg, 1) >= 2
            cleaned_segments{end + 1, 1} = seg; %#ok<AGROW>
        end
    end
end

function dist = point_to_polyline_distance(points, polyline)
    dist = inf(size(points, 1), 1);

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
        d = hypot(points(:, 1) - proj(:, 1), points(:, 2) - proj(:, 2));
        dist = min(dist, d);
    end
end

function step = diff_or_one(coords)
    if numel(coords) <= 1
        step = 1;
        return;
    end
    d = diff(coords);
    d = d(isfinite(d) & d > 0);
    if isempty(d)
        step = 1;
    else
        step = min(d);
    end
end
