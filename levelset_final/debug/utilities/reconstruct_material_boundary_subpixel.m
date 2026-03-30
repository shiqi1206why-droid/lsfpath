function geometry = reconstruct_material_boundary_subpixel(material_mask, dx, dy)
%RECONSTRUCT_MATERIAL_BOUNDARY_SUBPIXEL Reconstruct the mask boundary with
% marching squares and linear interpolation on the cell-centered mask field.

    mask = logical(material_mask);
    [nely, nelx] = size(mask);

    x_centers = ((1:nelx) - 0.5) * dx;
    y_centers = ((1:nely) - 0.5) * dy;

    % Add one-cell false padding so contour extraction is valid even when
    % material touches the design-domain boundary (including all-ones masks).
    mask_padded = false(nely + 2, nelx + 2);
    mask_padded(2:end-1, 2:end-1) = mask;
    x_padded = ((1:(nelx + 2)) - 1.5) * dx;
    y_padded = ((1:(nely + 2)) - 1.5) * dy;
    C = contourc(x_padded, y_padded, double(mask_padded), [0.5, 0.5]);
    raw_segments = parse_contourc_segments(C);
    segments = sanitize_boundary_segments(raw_segments, nelx, nely, dx, dy);

    boundary_mask = bwperim(mask);
    [pixel_y, pixel_x] = find(boundary_mask);

    geometry = struct();
    geometry.method = 'marching_squares_linear';
    geometry.mask_size = size(mask);
    geometry.segments = segments;
    geometry.segment_count = numel(segments);
    geometry.pixel_boundary = [(pixel_x - 0.5) * dx, (pixel_y - 0.5) * dy];
    geometry.x_centers = x_centers;
    geometry.y_centers = y_centers;
    geometry.domain_bounds = [0, nelx * dx, 0, nely * dy];
    geometry.success = ~isempty(segments);
end

function segments = sanitize_boundary_segments(raw_segments, nelx, nely, dx, dy)
    xmin = 0;
    xmax = nelx * dx;
    ymin = 0;
    ymax = nely * dy;
    tol = 1e-12;

    segments = {};
    for k = 1:numel(raw_segments)
        seg = raw_segments{k};
        if size(seg, 1) < 2
            continue;
        end

        finite_rows = all(isfinite(seg), 2);
        seg = seg(finite_rows, :);
        if size(seg, 1) < 2
            continue;
        end

        seg(:, 1) = min(max(seg(:, 1), xmin), xmax);
        seg(:, 2) = min(max(seg(:, 2), ymin), ymax);

        % Remove duplicated consecutive points caused by clamping.
        keep = [true; hypot(diff(seg(:, 1)), diff(seg(:, 2))) > tol];
        seg = seg(keep, :);
        if size(seg, 1) < 2
            continue;
        end

        % Remove tiny degenerate components.
        dxy = diff(seg, 1, 1);
        if sum(hypot(dxy(:, 1), dxy(:, 2))) <= max(dx, dy) * 1e-9
            continue;
        end

        segments{end + 1, 1} = seg; %#ok<AGROW>
    end
end
