function [lsf, parallel_paths, init_info] = construct_boundary_offset_levelset_with_parallel(material_mask, nelx, nely, dx, dy, delta_phi, smooth_opts)
    % 构造零等值线位于材料边界内偏移 Δφ 的水平集函数

    if nargin < 7 || isempty(smooth_opts)
        smooth_opts = struct();
    end

    mask = logical(material_mask);
    if isfield(smooth_opts, 'morph_radius') && smooth_opts.morph_radius > 0
        se = strel('disk', smooth_opts.morph_radius);
        mask = imopen(mask, se);
    end

    if ~any(mask(:))
        error('材料掩膜为空，无法构造水平集函数。');
    end

    h = min(dx, dy);
    boundary_reconstruction = lower(get_option_or_default( ...
        smooth_opts, 'boundary_reconstruction', 'marching_squares_linear'));
    switch boundary_reconstruction
        case 'marching_squares_linear'
            boundary_geometry = reconstruct_material_boundary_subpixel(mask, dx, dy);
        otherwise
            error('不支持的 boundary_reconstruction: %s', boundary_reconstruction);
    end

    if ~isfield(boundary_geometry, 'success') || ~boundary_geometry.success
        error('材料边界重建失败，无法初始化水平集。');
    end

    sign_reference = ones(size(mask));
    sign_reference(mask) = -1;
    [phi_boundary, boundary_distance_info] = build_signed_distance_from_segments( ...
        boundary_geometry.x_centers, boundary_geometry.y_centers, ...
        boundary_geometry.segments, sign_reference, true(size(mask)), 0);
    if ~boundary_distance_info.success
        error('无法根据子像素边界重建初始化符号距离场。');
    end

    max_inner_distance = max(-phi_boundary(mask));
    if ~isfinite(max_inner_distance) || max_inner_distance <= 0
        error('材料域的有效内距非正，无法构造边界内偏移路径。');
    end
    delta_phi_target = delta_phi;
    delta_phi_used = delta_phi_target;
    if delta_phi_used >= max_inner_distance
        delta_phi_used = 0.9 * max_inner_distance;
        warning('delta_phi %.4f 超出最大内距 %.4f，改用 %.4f。', ...
            delta_phi_target, max_inner_distance, delta_phi_used);
    end

    phi_main = phi_boundary + delta_phi_used;

    lsf = zeros(nely+2, nelx+2);
    lsf(2:end-1, 2:end-1) = phi_main;
    lsf(1, :) = lsf(2, :);
    lsf(end, :) = lsf(end-1, :);
    lsf(:, 1) = lsf(:, 2);
    lsf(:, end) = lsf(:, end-1);

    material_mask_full = expand_mask_to_full_grid(mask, size(lsf));

    zero_mask = false(size(lsf));
    zero_mask(2:end-1, 2:end-1) = (abs(phi_main) <= (0.5 * h)) & mask;
    zero_mask = zero_mask & material_mask_full;
    zero_mask = apply_neumann_extension(zero_mask);
    if ~any(zero_mask(:))
        [~, min_idx] = min(abs(phi_main(mask)));
        mask_linear = find(mask);
        selected_core = mask_linear(min_idx);
        [seed_i, seed_j] = ind2sub(size(mask), selected_core);
        zero_mask(seed_i + 1, seed_j + 1) = true;
        zero_mask = apply_neumann_extension(zero_mask);
    end

    reinit_opts = struct();
    reinit_opts.method = get_option_or_default(smooth_opts, 'reinit_method', 'subcell_signed_distance');
    if isfield(smooth_opts, 'zero_geometry_min_points') && ~isempty(smooth_opts.zero_geometry_min_points)
        reinit_opts.zero_geometry_min_points = smooth_opts.zero_geometry_min_points;
    end
    if isfield(smooth_opts, 'zero_geometry_min_length') && ~isempty(smooth_opts.zero_geometry_min_length)
        reinit_opts.zero_geometry_min_length = smooth_opts.zero_geometry_min_length;
    end
    [lsf, reinit_info] = fmm_reinitialize(lsf, dx, dy, zero_mask, mask, reinit_opts);

    spacing = h;
    max_distance = max(abs(lsf(material_mask_full)));
    num_levels = floor(max_distance / spacing);
    base_levels = spacing * (1:min(num_levels, 3));

    parallel_paths = struct();
    parallel_paths.positive_levels = base_levels;
    parallel_paths.negative_levels = -base_levels;
    parallel_paths.spacing = spacing;
    parallel_paths.zero_mask = zero_mask;
    parallel_paths.boundary_geometry = boundary_geometry;

    init_info = compute_boundary_offset_stats(lsf, mask, dx, dy, delta_phi_target, delta_phi_used, ...
        max_inner_distance, boundary_geometry, phi_boundary);
    init_info.zero_mask = zero_mask;
    init_info.spacing = spacing;
    init_info.delta_phi_target = delta_phi_target;
    init_info.delta_phi_used = delta_phi_used;
    init_info.phi_boundary_core = phi_boundary;
    init_info.phi_boundary_full = zeros(size(lsf));
    init_info.phi_boundary_full(2:end-1, 2:end-1) = phi_boundary;
    init_info.phi_boundary_full = apply_neumann_extension(init_info.phi_boundary_full);
    init_info.boundary_geometry = boundary_geometry;
    init_info.boundary_reconstruction_method = boundary_reconstruction;
    init_info.boundary_distance_info = boundary_distance_info;
    init_info.material_mask_full = material_mask_full;
    init_info.reinit_info = reinit_info;
end

function mask_full = expand_mask_to_full_grid(mask_core, target_size)
    mask_full = false(target_size);
    mask_full(2:end-1, 2:end-1) = logical(mask_core);
    mask_full = apply_neumann_extension(mask_full);
end

function field = apply_neumann_extension(field)
    field(1, :) = field(2, :);
    field(end, :) = field(end-1, :);
    field(:, 1) = field(:, 2);
    field(:, end) = field(:, end-1);
end

function value = get_option_or_default(opts, name, default_value)
    value = default_value;
    if isfield(opts, name) && ~isempty(opts.(name))
        value = opts.(name);
    end
end

