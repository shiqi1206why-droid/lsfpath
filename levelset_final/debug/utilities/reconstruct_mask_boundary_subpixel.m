function geometry = reconstruct_mask_boundary_subpixel(material_mask, dx, dy)
%RECONSTRUCT_MASK_BOUNDARY_SUBPIXEL Extract subcell mask boundary using contouring.
%
% This is equivalent to marching squares with linear interpolation on the
% binary material indicator sampled at cell centers.

    mask = logical(material_mask);
    [nely, nelx] = size(mask);
    x_core = (0.5:1:(nelx - 0.5)) * dx;
    y_core = (0.5:1:(nely - 0.5)) * dy;

    C = contourc(x_core, y_core, double(mask), [0.5, 0.5]);
    [segments, levels] = parse_contourc_segments(C); %#ok<ASGLU>

    boundary_mask = bwperim(mask);
    [py, px] = find(boundary_mask);

    geometry = struct();
    geometry.method = 'marching_squares_linear';
    geometry.segments = segments;
    geometry.segment_count = numel(segments);
    geometry.point_count = sum(cellfun(@(seg) size(seg, 1), segments));
    geometry.pixel_boundary_x = (px - 0.5) * dx;
    geometry.pixel_boundary_y = (py - 0.5) * dy;
    geometry.x_core = x_core;
    geometry.y_core = y_core;

    [geom_x, geom_y] = concatenate_segments(segments);
    geometry.boundary_x = geom_x;
    geometry.boundary_y = geom_y;
end

function [x_all, y_all] = concatenate_segments(segments)
    x_all = zeros(0, 1);
    y_all = zeros(0, 1);
    for idx = 1:numel(segments)
        seg = segments{idx};
        if isempty(seg)
            continue;
        end
        if ~isempty(x_all)
            x_all(end + 1, 1) = NaN; %#ok<AGROW>
            y_all(end + 1, 1) = NaN; %#ok<AGROW>
        end
        x_all = [x_all; seg(:, 1)]; %#ok<AGROW>
        y_all = [y_all; seg(:, 2)]; %#ok<AGROW>
    end
end
