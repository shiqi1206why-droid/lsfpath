function init_info = compute_boundary_offset_stats(lsf, material_mask, dx, dy, delta_phi_target, delta_phi_used, max_inner_distance)
    % 统计边界偏移初始化的各项诊断信息

    mask = logical(material_mask);
    h = min(dx, dy);

    d_in = bwdist(~mask) * h;
    d_out = bwdist(mask) * h;
    phi_boundary_core = d_out - d_in;  % Inner negative, outer positive

    phi_boundary = zeros(size(lsf));
    phi_boundary(2:end-1, 2:end-1) = phi_boundary_core;
    phi_boundary(1, :) = phi_boundary(2, :);
    phi_boundary(end, :) = phi_boundary(end-1, :);
    phi_boundary(:, 1) = phi_boundary(:, 2);
    phi_boundary(:, end) = phi_boundary(:, end-1);

    Lx = dx * (size(lsf, 2) - 2);
    Ly = dy * (size(lsf, 1) - 2);
    x_coords = linspace(0, Lx, size(lsf, 2));
    y_coords = linspace(0, Ly, size(lsf, 1));
    C = contourc(x_coords, y_coords, lsf, [0 0]);
    [contour_x, contour_y] = contourc_to_points(C);

    [Xgrid, Ygrid] = meshgrid(x_coords, y_coords);
    distances = [];
    if ~isempty(contour_x)
        % 使用nearest插值避免边界NaN，保留所有采样点
        phi_interp = interp2(Xgrid, Ygrid, phi_boundary, contour_x, contour_y, 'nearest', 0);
        distances = abs(phi_interp);
        
        % === 统计量NaN保护：过滤非有限值 ===
        % 参考：代码修改交流.md (Codex方案 2025-10-17)
        distances = distances(isfinite(distances));
    end

    if isempty(distances)
        mean_offset = 0;
        std_offset = 0;
        mean_error = 0;
        max_error = 0;
        distance_samples = [];
        sample_indices = [];
    else
        mean_offset = mean(distances);
        std_offset = std(distances);
        mean_error = mean(distances - delta_phi_used);
        max_error = max(abs(distances - delta_phi_used));
        samples_limit = 1000;
        if numel(distances) > samples_limit
            sample_indices = round(linspace(1, numel(distances), samples_limit));
        else
            sample_indices = 1:numel(distances);
        end
        distance_samples = distances(sample_indices);
    end
    
    % 保存完整的等值线坐标用于可视化（不过滤）

    thin_mask = mask & (d_in < delta_phi_used);
    total_cells = max(1, nnz(mask));
    thin_ratio = nnz(thin_mask) / total_cells;

    boundary_mask = bwperim(mask);
    [by, bx] = find(boundary_mask);
    boundary_coords_x = (bx - 0.5) * dx;
    boundary_coords_y = (by - 0.5) * dy;

    init_info = struct();
    init_info.mean_offset = mean_offset;
    init_info.std_offset = std_offset;
    init_info.mean_error = mean_error;
    init_info.max_error = max_error;
    init_info.num_samples = numel(distances);
    init_info.distance_samples = distance_samples;
    init_info.distance_sample_indices = sample_indices;
    init_info.distance_sample_errors = distance_samples - delta_phi_used;
    init_info.contour = struct('x', contour_x(:)', 'y', contour_y(:)');
    init_info.boundary_coords = struct('x', boundary_coords_x(:)', 'y', boundary_coords_y(:)');
    init_info.thin_ratio = thin_ratio;
    init_info.num_thin_cells = nnz(thin_mask);
    init_info.total_material_cells = total_cells;
    if total_cells > 0
        init_info.min_thickness = 2 * min(d_in(mask));
    else
        init_info.min_thickness = 0;
    end
    init_info.max_inner_distance = max_inner_distance;
    init_info.delta_phi_target = delta_phi_target;
    init_info.delta_phi_used = delta_phi_used;
end

