function boundary_geometry = reconstruct_subcell_boundary_geometry(material_mask, dx, dy)
    % 使用 marching squares + 线性插值重建材料域子像素边界

    mask = logical(material_mask);
    if ~any(mask(:))
        error('reconstruct_subcell_boundary_geometry: 材料掩膜为空。');
    end

    [nely, nelx] = size(mask);
    x_centers = ((1:nelx) - 0.5) * dx;
    y_centers = ((1:nely) - 0.5) * dy;

    contour_matrix = contourc(x_centers, y_centers, double(mask), [0.5, 0.5]);
    components = parse_contour_components(contour_matrix);
    if isempty(components)
        if all(mask(:))
            rect = [0, 0; nelx * dx, 0; nelx * dx, nely * dy; 0, nely * dy; 0, 0];
            components = {rect};
        else
            error('reconstruct_subcell_boundary_geometry: 未能提取材料边界。');
        end
    end

    all_x = [];
    all_y = [];
    component_lengths = zeros(numel(components), 1);
    valid_components = cell(0, 1);
    segment_start = zeros(0, 2);
    segment_end = zeros(0, 2);

    for k = 1:numel(components)
        polyline = components{k};
        if size(polyline, 1) < 2
            continue;
        end
        if norm(polyline(1, :) - polyline(end, :)) > 1e-12
            polyline(end + 1, :) = polyline(1, :);
        end
        valid_components{end + 1, 1} = polyline; %#ok<AGROW>

        diffs = diff(polyline, 1, 1);
        component_lengths(numel(valid_components)) = sum(hypot(diffs(:, 1), diffs(:, 2)));
        segment_start = [segment_start; polyline(1:end-1, :)]; %#ok<AGROW>
        segment_end = [segment_end; polyline(2:end, :)]; %#ok<AGROW>

        if ~isempty(all_x)
            all_x(end + 1, 1) = NaN; %#ok<AGROW>
            all_y(end + 1, 1) = NaN; %#ok<AGROW>
        end
        all_x = [all_x; polyline(:, 1)]; %#ok<AGROW>
        all_y = [all_y; polyline(:, 2)]; %#ok<AGROW>
    end

    component_lengths = component_lengths(1:numel(valid_components));
    if isempty(valid_components) || isempty(segment_start)
        error('reconstruct_subcell_boundary_geometry: 边界折线退化，无法构建有效线段。');
    end

    boundary_geometry = struct();
    boundary_geometry.method = 'marching_squares_linear';
    boundary_geometry.level = 0.5;
    boundary_geometry.component_count = numel(valid_components);
    boundary_geometry.components = valid_components;
    boundary_geometry.x = all_x(:)';
    boundary_geometry.y = all_y(:)';
    boundary_geometry.segment_start = segment_start;
    boundary_geometry.segment_end = segment_end;
    boundary_geometry.component_lengths = component_lengths(:)';
    boundary_geometry.total_length = sum(component_lengths);
    boundary_geometry.x_centers = x_centers(:)';
    boundary_geometry.y_centers = y_centers(:)';
    boundary_geometry.mask_size = [nely, nelx];
end

function components = parse_contour_components(C)
    components = {};
    if isempty(C)
        return;
    end

    idx = 1;
    while idx < size(C, 2)
        point_count = C(2, idx);
        polyline = C(:, idx + 1:idx + point_count).';
        if size(polyline, 1) >= 2
            components{end + 1, 1} = polyline; %#ok<AGROW>
        end
        idx = idx + point_count + 1;
    end
end
