function contour_info = extract_subcell_isocontours(field, level, dx, dy, valid_mask)
%EXTRACT_SUBCELL_ISOCONTOURS Extract subcell contours using marching squares.
% contourc performs a marching-squares style linear interpolation on the
% scalar field, which is the desired subcell contour reconstruction here.

    if nargin < 5 || isempty(valid_mask)
        valid_mask = true(size(field));
    end

    if ~isequal(size(field), size(valid_mask))
        error('field and valid_mask must have identical sizes.');
    end

    valid_mask = logical(valid_mask);
    [x_coords, y_coords] = get_lsf_grid_coordinates(size(field), dx, dy);

    masked_field = field;
    finite_valid_mask = valid_mask & isfinite(masked_field);
    finite_valid = masked_field(finite_valid_mask);
    if isempty(finite_valid)
        contour_info = struct();
        contour_info.level = level;
        contour_info.segments = {};
        contour_info.segment_lengths = zeros(0, 1);
        contour_info.segment_count = 0;
        contour_info.raw_segment_count = 0;
        contour_info.success = false;
        contour_info.x_coords = x_coords;
        contour_info.y_coords = y_coords;
        return;
    end

    % Fill invalid/outside points by nearest valid samples to avoid creating
    % artificial contour components along the valid-mask boundary.
    if any(~finite_valid_mask(:))
        [~, nearest_idx] = bwdist(finite_valid_mask);
        replace_idx = ~finite_valid_mask;
        masked_field(replace_idx) = masked_field(nearest_idx(replace_idx));
    end

    C = contourc(x_coords, y_coords, masked_field, [level, level]);
    raw_segments = parse_contourc_segments(C);
    segments = clip_segments_to_valid_mask(raw_segments, x_coords, y_coords, valid_mask);

    lengths = zeros(numel(segments), 1);
    for k = 1:numel(segments)
        dxy = diff(segments{k}, 1, 1);
        lengths(k) = sum(hypot(dxy(:, 1), dxy(:, 2)));
    end

    contour_info = struct();
    contour_info.level = level;
    contour_info.segments = segments;
    contour_info.segment_lengths = lengths;
    contour_info.segment_count = numel(segments);
    contour_info.raw_segment_count = numel(raw_segments);
    contour_info.success = ~isempty(segments);
    contour_info.x_coords = x_coords;
    contour_info.y_coords = y_coords;
end

function segments = clip_segments_to_valid_mask(raw_segments, x_coords, y_coords, valid_mask)
    mask_values = double(valid_mask);
    segments = {};
    for k = 1:numel(raw_segments)
        seg = raw_segments{k};
        if size(seg, 1) < 2
            continue;
        end
        clipped_parts = split_segment_by_mask(seg, x_coords, y_coords, mask_values);
        for p = 1:numel(clipped_parts)
            part = clipped_parts{p};
            if size(part, 1) < 2
                continue;
            end
            keep = [true; hypot(diff(part(:, 1)), diff(part(:, 2))) > 1e-12];
            part = part(keep, :);
            if size(part, 1) >= 2
                segments{end + 1, 1} = part; %#ok<AGROW>
            end
        end
    end
end

function parts = split_segment_by_mask(seg, x_coords, y_coords, mask_values)
    vals = interp2(x_coords, y_coords, mask_values, seg(:, 1), seg(:, 2), 'linear', 0);
    vals(~isfinite(vals)) = 0;
    inside = vals >= 0.5;

    parts = {};
    if ~any(inside)
        return;
    end

    current = zeros(0, 2);
    npt = size(seg, 1);
    for i = 1:(npt - 1)
        p0 = seg(i, :);
        p1 = seg(i + 1, :);
        v0 = vals(i);
        v1 = vals(i + 1);
        in0 = inside(i);
        in1 = inside(i + 1);

        if in0 && isempty(current)
            current = p0;
        end

        if in0 && in1
            current = [current; p1]; %#ok<AGROW>
            continue;
        end

        if in0 && ~in1
            p_cross = interpolate_to_mask_half(p0, p1, v0, v1);
            current = [current; p_cross]; %#ok<AGROW>
            if size(current, 1) >= 2
                parts{end + 1, 1} = current; %#ok<AGROW>
            end
            current = zeros(0, 2);
            continue;
        end

        if ~in0 && in1
            p_cross = interpolate_to_mask_half(p0, p1, v0, v1);
            current = [p_cross; p1];
            continue;
        end
    end

    if inside(end)
        if isempty(current)
            current = seg(end, :);
        end
        if size(current, 1) >= 2
            parts{end + 1, 1} = current; %#ok<AGROW>
        end
    end
end

function p = interpolate_to_mask_half(p0, p1, v0, v1)
    if ~isfinite(v0) || ~isfinite(v1) || abs(v1 - v0) < 1e-12
        t = 0.5;
    else
        t = (0.5 - v0) / (v1 - v0);
        t = min(max(t, 0), 1);
    end
    p = p0 + t * (p1 - p0);
end
