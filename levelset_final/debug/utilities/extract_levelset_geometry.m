function geometry = extract_levelset_geometry(lsf, arg2, arg3, arg4, material_mask)
    % 使用 contourc 提取子像素等值线几何

    if nargin < 5
        material_mask = [];
    end

    if arg2 > 0 && arg3 > 0
        dx = arg2;
        dy = arg3;
        level = arg4;
    else
        level = arg2;
        dx = arg3;
        dy = arg4;
    end

    [x_full, y_full] = get_lsf_grid_coordinates(size(lsf), dx, dy);
    contour_matrix = contourc(x_full, y_full, lsf, [level, level]);
    segments = parse_contourc_segments_local(contour_matrix);

    if ~isempty(material_mask)
        material_mask = normalize_mask_to_lsf_grid_local(material_mask, size(lsf));
        segments = filter_segments_by_material_local(segments, material_mask, x_full, y_full);
    end

    segments = segments(cellfun(@(s) size(s, 1) >= 2, segments));
    geometry = build_geometry_struct_local(segments, level);
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
        if npts < 2 || j1 > ncol
            break;
        end
        segments{end + 1, 1} = C(:, j0:j1).'; %#ok<AGROW>
        idx = j1 + 1;
    end
end

function geometry = build_geometry_struct_local(segments, level)
    geometry = struct();
    geometry.level = level;
    geometry.segments = segments;
    geometry.component_count = numel(segments);
    geometry.segment_count = numel(segments);
    geometry.segment_start = zeros(0, 2);
    geometry.segment_end = zeros(0, 2);
    geometry.component_lengths = zeros(1, numel(segments));
    geometry.total_length = 0;
    geometry.x = [];
    geometry.y = [];
    geometry.point_count = 0;

    all_x = [];
    all_y = [];
    seg_start = zeros(0, 2);
    seg_end = zeros(0, 2);
    lengths = zeros(1, numel(segments));
    total_length = 0;

    for k = 1:numel(segments)
        polyline = segments{k};
        if isempty(polyline)
            continue;
        end
        diffs = diff(polyline, 1, 1);
        seg_len = sum(hypot(diffs(:, 1), diffs(:, 2)));
        lengths(k) = seg_len;
        total_length = total_length + seg_len;
        seg_start = [seg_start; polyline(1:end-1, :)]; %#ok<AGROW>
        seg_end = [seg_end; polyline(2:end, :)]; %#ok<AGROW>
        if ~isempty(all_x)
            all_x(end + 1, 1) = NaN; %#ok<AGROW>
            all_y(end + 1, 1) = NaN; %#ok<AGROW>
        end
        all_x = [all_x; polyline(:, 1)]; %#ok<AGROW>
        all_y = [all_y; polyline(:, 2)]; %#ok<AGROW>
    end

    geometry.segment_start = seg_start;
    geometry.segment_end = seg_end;
    geometry.component_lengths = lengths;
    geometry.total_length = total_length;
    geometry.x = all_x(:)';
    geometry.y = all_y(:)';
    geometry.point_count = sum(cellfun(@(s) size(s, 1), segments));
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

    error('extract_levelset_geometry: material_mask尺寸不匹配。');
end

function segments = filter_segments_by_material_local(segments, material_mask_full, x_full, y_full)
    filtered = {};
    mask_values = double(material_mask_full);
    for k = 1:numel(segments)
        polyline = segments{k};
        if size(polyline, 1) < 2
            continue;
        end
        sample_x = polyline(:, 1);
        sample_y = polyline(:, 2);
        inside = interp2(x_full, y_full, mask_values, sample_x, sample_y, 'nearest', 0) > 0.5;
        if any(inside)
            filtered{end + 1, 1} = polyline; %#ok<AGROW>
        end
    end
    segments = filtered;
end
