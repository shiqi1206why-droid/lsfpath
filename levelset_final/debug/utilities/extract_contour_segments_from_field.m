function contour_info = extract_contour_segments_from_field(field, x_coords, y_coords, level)
%EXTRACT_CONTOUR_SEGMENTS_FROM_FIELD Extract contour polylines at a level.
%   contour_info.segments is a cell array of [x,y] polylines.

    if nargin < 4 || isempty(level)
        level = 0;
    end

    C = contourc(x_coords, y_coords, double(field), [level, level]);
    segments = parse_contourc_segments_local(C);

    all_x = [];
    all_y = [];
    valid_segments = {};
    for i = 1:numel(segments)
        seg = segments{i};
        if size(seg, 1) < 2 || any(~isfinite(seg(:)))
            continue;
        end
        valid_segments{end+1} = seg; %#ok<AGROW>
        all_x = [all_x; seg(:, 1)]; %#ok<AGROW>
        all_y = [all_y; seg(:, 2)]; %#ok<AGROW>
    end

    contour_info = struct();
    contour_info.level = level;
    contour_info.segments = valid_segments;
    contour_info.count = numel(valid_segments);
    contour_info.x = all_x(:)';
    contour_info.y = all_y(:)';
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
        if npts <= 1 || j1 > ncol
            break;
        end
        segments{end+1} = C(:, j0:j1)'; %#ok<AGROW>
        idx = j1 + 1;
    end
end
