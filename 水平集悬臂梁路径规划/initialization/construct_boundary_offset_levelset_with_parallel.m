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
    d_in = bwdist(~mask) * h;
    d_out = bwdist(mask) * h;
    phi_boundary = d_out - d_in;  % Inner negative, outer positive

    max_inner_distance = max(d_in(:));
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

    zero_mask = false(size(lsf));
    zero_mask(2:end-1, 2:end-1) = abs(phi_main) <= (0.5 * h);
    zero_mask(1, :) = zero_mask(2, :);
    zero_mask(end, :) = zero_mask(end-1, :);
    zero_mask(:, 1) = zero_mask(:, 2);
    zero_mask(:, end) = zero_mask(:, end-1);
    % 在全域传播，避免外圈被掩膜截断
    lsf = fmm_reinitialize(lsf, dx, dy, zero_mask, []);

    spacing = h;
    max_distance = max(abs(lsf(:)));
    num_levels = floor(max_distance / spacing);
    base_levels = spacing * (1:min(num_levels, 3));

    parallel_paths = struct();
    parallel_paths.positive_levels = base_levels;
    parallel_paths.negative_levels = -base_levels;
    parallel_paths.spacing = spacing;
    parallel_paths.zero_mask = zero_mask;

    init_info = compute_boundary_offset_stats(lsf, mask, dx, dy, delta_phi_target, delta_phi_used, max_inner_distance);
    init_info.zero_mask = zero_mask;
    init_info.spacing = spacing;
    init_info.delta_phi_target = delta_phi_target;
    init_info.delta_phi_used = delta_phi_used;
end

