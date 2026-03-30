function [lsf, diagnostics] = fmm_reinitialize(lsf, dx, dy, zero_mask, material_mask, opts)
%FMM_REINITIALIZE Reinitialize the level set with a geometry-first strategy.
% Primary path:
%   1) extract subcell zero contours from the current lsf
%   2) rebuild the signed distance field from those contours
% Fallback:
%   masked fast marching with the provided zero seeds

    if nargin < 5 || isempty(material_mask)
        material_mask = true(size(lsf));
    end
    if nargin < 6 || isempty(opts)
        opts = struct();
    end
    opts = apply_reinit_defaults(opts);

    if nargin < 4 || isempty(zero_mask)
        tol_zero = 0.5 * min(dx, dy);
        zero_mask = abs(lsf) <= tol_zero;
    else
        zero_mask = logical(zero_mask);
    end

    material_mask = normalize_mask_to_lsf_grid(material_mask, size(lsf), 'material_mask');
    zero_mask = normalize_mask_to_lsf_grid(zero_mask, size(lsf), 'zero_mask');
    zero_mask = zero_mask & material_mask;
    local_shell_mask = [];
    if ~isempty(opts.local_shell_mask)
        local_shell_mask = normalize_mask_to_lsf_grid(opts.local_shell_mask, size(lsf), 'local_shell_mask');
        local_shell_mask = local_shell_mask & material_mask;
    end

    diagnostics = struct();
    diagnostics.requested_method = opts.method;
    diagnostics.method_used = '';
    diagnostics.used_method = '';
    diagnostics.fallback_used = false;
    diagnostics.fallback_reason = '';
    diagnostics.zero_seed_count = nnz(zero_mask);
    diagnostics.contour_segment_count = 0;
    diagnostics.outside_distance = NaN;
    diagnostics.local_shell_requested = ~isempty(local_shell_mask) && any(local_shell_mask(:));
    diagnostics.local_shell_size = nnz(local_shell_mask);
    diagnostics.local_shell_applied = false;
    diagnostics.preserve_outside_shell = logical(opts.preserve_outside_shell);

    method_name = lower(strtrim(string(opts.method)));
    if method_name == "subcell_signed_distance"
        [lsf_geom, geom_info] = reinitialize_from_zero_geometry(lsf, dx, dy, material_mask, opts);
        diagnostics.geometry = geom_info;
        diagnostics.contour_segment_count = geom_info.segment_count;
        if geom_info.success
            diagnostics.method_used = 'subcell_signed_distance';
            diagnostics.used_method = diagnostics.method_used;
            diagnostics.outside_distance = geom_info.outside_distance;
            lsf = merge_local_shell_update(lsf, lsf_geom, material_mask, local_shell_mask, ...
                logical(opts.preserve_outside_shell), geom_info.outside_distance);
            diagnostics.local_shell_applied = diagnostics.local_shell_requested;
            return;
        end
        diagnostics.fallback_used = true;
        diagnostics.fallback_reason = geom_info.failure_reason;
    elseif method_name == "masked_fmm" || method_name == "masked_fmm_fallback"
        [lsf_candidate, fmm_info] = masked_fmm_reinitialize(lsf, dx, dy, zero_mask, material_mask);
        diagnostics.method_used = 'masked_fmm';
        diagnostics.used_method = diagnostics.method_used;
        diagnostics.fmm = fmm_info;
        diagnostics.outside_distance = fmm_info.outside_distance;
        lsf = merge_local_shell_update(lsf, lsf_candidate, material_mask, local_shell_mask, ...
            logical(opts.preserve_outside_shell), fmm_info.outside_distance);
        diagnostics.local_shell_applied = diagnostics.local_shell_requested;
        return;
    else
        diagnostics.fallback_used = true;
        diagnostics.fallback_reason = 'unsupported_reinit_method';
    end

    [lsf_candidate, fmm_info] = masked_fmm_reinitialize(lsf, dx, dy, zero_mask, material_mask);
    diagnostics.method_used = 'masked_fmm_fallback';
    diagnostics.used_method = diagnostics.method_used;
    diagnostics.fmm = fmm_info;
    diagnostics.outside_distance = fmm_info.outside_distance;
    lsf = merge_local_shell_update(lsf, lsf_candidate, material_mask, local_shell_mask, ...
        logical(opts.preserve_outside_shell), fmm_info.outside_distance);
    diagnostics.local_shell_applied = diagnostics.local_shell_requested;
end

function opts = apply_reinit_defaults(opts)
    if ~isfield(opts, 'method') || isempty(opts.method)
        opts.method = 'subcell_signed_distance';
    end
    if ~isfield(opts, 'zero_geometry_min_points') || isempty(opts.zero_geometry_min_points)
        opts.zero_geometry_min_points = 8;
    end
    if ~isfield(opts, 'zero_geometry_min_length') || isempty(opts.zero_geometry_min_length)
        opts.zero_geometry_min_length = 0.5;
    end
    if ~isfield(opts, 'local_shell_mask') || isempty(opts.local_shell_mask)
        opts.local_shell_mask = [];
    end
    if ~isfield(opts, 'preserve_outside_shell') || isempty(opts.preserve_outside_shell)
        opts.preserve_outside_shell = false;
    end
end

function [lsf_new, info] = reinitialize_from_zero_geometry(lsf, dx, dy, material_mask, opts)
    contour_info = extract_subcell_isocontours(lsf, 0, dx, dy, material_mask);
    info = struct();
    info.success = false;
    info.method = 'subcell_signed_distance';
    info.segment_count = contour_info.segment_count;
    info.raw_segment_count = get_field_or_default(contour_info, 'raw_segment_count', contour_info.segment_count);
    info.failure_reason = '';
    info.outside_distance = NaN;
    info.total_points = 0;
    info.total_length = 0;

    if ~contour_info.success || contour_info.segment_count == 0
        info.failure_reason = 'empty_zero_contour';
        lsf_new = lsf;
        return;
    end
    total_points = 0;
    for si = 1:numel(contour_info.segments)
        total_points = total_points + size(contour_info.segments{si}, 1);
    end
    info.total_points = total_points;
    if total_points < opts.zero_geometry_min_points
        info.failure_reason = 'insufficient_zero_points';
        lsf_new = lsf;
        return;
    end
    if isfield(contour_info, 'segment_lengths')
        total_length = sum(contour_info.segment_lengths);
    else
        total_length = NaN;
    end
    info.total_length = total_length;
    if isfinite(total_length) && total_length < opts.zero_geometry_min_length * min(dx, dy)
        info.failure_reason = 'short_zero_contour';
        lsf_new = lsf;
        return;
    end

    sign_reference = sign(lsf);
    sign_reference(sign_reference == 0) = 1;
    sign_reference(~material_mask) = 1;

    finite_material = abs(lsf(material_mask & isfinite(lsf)));
    if isempty(finite_material)
        outside_distance = 10 * max(dx, dy);
    else
        outside_distance = max(finite_material) + 10 * max(dx, dy);
    end
    info.outside_distance = outside_distance;

    [x_coords, y_coords] = get_lsf_grid_coordinates(size(lsf), dx, dy);
    [lsf_new, dist_info] = build_signed_distance_from_segments( ...
        x_coords, y_coords, contour_info.segments, sign_reference, material_mask, outside_distance);
    if ~dist_info.success
        info.failure_reason = 'distance_rebuild_failed';
        lsf_new = lsf;
        return;
    end

    lsf_new(~material_mask) = outside_distance;
    lsf_new = apply_neumann_boundary(lsf_new);
    info.success = true;
    info.distance_info = dist_info;
end

function [lsf, info] = masked_fmm_reinitialize(lsf, dx, dy, zero_mask, material_mask)
    if ~any(zero_mask(:))
        tol_zero = 0.5 * min(dx, dy);
        zero_mask = (abs(lsf) <= tol_zero) & material_mask;
    end

    seed_fallback_used = false;
    if ~any(zero_mask(:))
        interior_mask = material_mask(2:end-1, 2:end-1);
        if ~any(interior_mask(:))
            error('fmm_reinitialize: 材料域为空，无法执行masked重初始化。');
        end
        core_candidates = find(interior_mask);
        core_lsf = abs(lsf(2:end-1, 2:end-1));
        [~, best_local_idx] = min(core_lsf(core_candidates));
        best_core_linear = core_candidates(best_local_idx);
        [seed_i, seed_j] = ind2sub(size(interior_mask), best_core_linear);
        zero_mask(seed_i + 1, seed_j + 1) = true;
        seed_fallback_used = true;
    end

    sign_field = sign(lsf);
    sign_field(sign_field == 0) = 1;
    sign_field(~material_mask) = 1;

    T = fmm_compute_distance(zero_mask, dx, dy, material_mask);
    finite_material = T(isfinite(T) & material_mask);
    if isempty(finite_material)
        Tmax = 0;
    else
        Tmax = max(finite_material);
    end
    pad = max(10 * max(dx, dy), min(dx, dy));
    outside_distance = max(Tmax + pad, pad);

    unreachable_inside = ~isfinite(T) & material_mask;
    T(unreachable_inside) = outside_distance;
    T(~material_mask) = outside_distance;

    lsf = sign_field .* T;
    lsf(zero_mask) = 0;
    lsf(~material_mask) = outside_distance;
    lsf = apply_neumann_boundary(lsf);

    info = struct();
    info.zero_seed_count = nnz(zero_mask);
    info.seed_fallback_used = seed_fallback_used;
    info.unreachable_inside_count = nnz(unreachable_inside);
    info.outside_distance = outside_distance;
end

function lsf_out = merge_local_shell_update(lsf_old, lsf_candidate, material_mask, local_shell_mask, preserve_outside_shell, outside_distance)
    if nargin < 6 || isempty(outside_distance) || ~isfinite(outside_distance)
        outside_distance = max(10, 10 * max(1, max(abs(lsf_candidate(material_mask)), [], 'omitnan')));
    end

    if preserve_outside_shell && ~isempty(local_shell_mask)
        lsf_out = lsf_old;
        lsf_out(local_shell_mask) = lsf_candidate(local_shell_mask);
    else
        lsf_out = lsf_candidate;
    end

    lsf_out(~material_mask) = outside_distance;
    lsf_out = apply_neumann_boundary(lsf_out);
end

function value = get_field_or_default(data, field_name, default_value)
    value = default_value;
    if isfield(data, field_name) && ~isempty(data.(field_name))
        value = data.(field_name);
    end
end

%% ================== FMM核心算法 ==================

function T = fmm_compute_distance(zero_mask, dx, dy, material_mask)
% FMM_COMPUTE_DISTANCE 使用快速行进法计算到达时间场（距离场）

    [ny, nx] = size(zero_mask);

    T = inf(ny, nx);
    State = zeros(ny, nx, 'uint8'); % 0=Far, 1=Trial, 2=Alive
    F_s = 1.0;

    seed_points = find(zero_mask & material_mask);
    if isempty(seed_points)
        return;
    end
    T(seed_points) = 0;
    State(seed_points) = 2;

    heap = [];
    for idx = 1:length(seed_points)
        sp = seed_points(idx);
        [i, j] = ind2sub([ny, nx], sp);
        neighbors = get_neighbors(i, j, ny, nx);
        for n = 1:size(neighbors, 1)
            ni = neighbors(n, 1);
            nj = neighbors(n, 2);
            nidx = sub2ind([ny, nx], ni, nj);

            if State(ni, nj) == 0 && material_mask(ni, nj)
                T_new = compute_local_update(T, ni, nj, dx, dy, F_s, State);
                T(ni, nj) = T_new;
                State(ni, nj) = 1;
                heap = heap_push(heap, T_new, nidx);
            end
        end
    end

    while ~isempty(heap)
        [heap, ~, idx_min] = heap_pop(heap);
        [i, j] = ind2sub([ny, nx], idx_min);
        State(i, j) = 2;

        neighbors = get_neighbors(i, j, ny, nx);
        for n = 1:size(neighbors, 1)
            ni = neighbors(n, 1);
            nj = neighbors(n, 2);
            nidx = sub2ind([ny, nx], ni, nj);

            if State(ni, nj) ~= 2 && material_mask(ni, nj)
                T_new = compute_local_update(T, ni, nj, dx, dy, F_s, State);

                if State(ni, nj) == 0
                    T(ni, nj) = T_new;
                    State(ni, nj) = 1;
                    heap = heap_push(heap, T_new, nidx);
                elseif T_new < T(ni, nj)
                    T(ni, nj) = T_new;
                    heap = heap_decrease_key(heap, T_new, nidx);
                end
            end
        end
    end
end

function T_new = compute_local_update(T, i, j, dx, dy, F_s, State)
    [ny, nx] = size(T);

    T_L = inf; T_R = inf; T_D = inf; T_U = inf;
    if j > 1 && State(i, j-1) == 2
        T_L = T(i, j-1);
    end
    if j < nx && State(i, j+1) == 2
        T_R = T(i, j+1);
    end
    if i > 1 && State(i-1, j) == 2
        T_D = T(i-1, j);
    end
    if i < ny && State(i+1, j) == 2
        T_U = T(i+1, j);
    end

    a = min(T_L, T_R);
    b = min(T_D, T_U);

    if isinf(a) && isinf(b)
        T_new = inf;
        return;
    elseif isinf(a)
        T_new = b + dy / F_s;
        return;
    elseif isinf(b)
        T_new = a + dx / F_s;
        return;
    end

    A = 1/dx^2 + 1/dy^2;
    B = -2*(a/dx^2 + b/dy^2);
    C = a^2/dx^2 + b^2/dy^2 - 1/F_s^2;
    discriminant = B^2 - 4*A*C;

    if discriminant >= 0
        T_candidate = (-B + sqrt(discriminant)) / (2*A);
        if T_candidate >= max(a, b)
            T_new = T_candidate;
        else
            T_new = min(a + dx/F_s, b + dy/F_s);
        end
    else
        T_new = min(a + dx/F_s, b + dy/F_s);
    end
end

function neighbors = get_neighbors(i, j, ny, nx)
    candidates = [i-1, j; i+1, j; i, j-1; i, j+1];
    valid = candidates(:,1) >= 1 & candidates(:,1) <= ny & ...
            candidates(:,2) >= 1 & candidates(:,2) <= nx;
    neighbors = candidates(valid, :);
end

function heap = heap_push(heap, T_value, idx)
    heap = [heap; T_value, idx]; %#ok<AGROW>
    pos = size(heap, 1);
    while pos > 1
        parent = floor(pos / 2);
        if heap(parent, 1) <= heap(pos, 1)
            break;
        end
        tmp = heap(parent, :);
        heap(parent, :) = heap(pos, :);
        heap(pos, :) = tmp;
        pos = parent;
    end
end

function [heap, T_value, idx] = heap_pop(heap)
    if isempty(heap)
        T_value = [];
        idx = [];
        return;
    end
    T_value = heap(1, 1);
    idx = heap(1, 2);
    heap(1, :) = heap(end, :);
    heap(end, :) = [];
    pos = 1;
    while true
        left = 2 * pos;
        right = left + 1;
        smallest = pos;
        if left <= size(heap, 1) && heap(left, 1) < heap(smallest, 1)
            smallest = left;
        end
        if right <= size(heap, 1) && heap(right, 1) < heap(smallest, 1)
            smallest = right;
        end
        if smallest == pos
            break;
        end
        tmp = heap(pos, :);
        heap(pos, :) = heap(smallest, :);
        heap(smallest, :) = tmp;
        pos = smallest;
    end
end

function heap = heap_decrease_key(heap, T_value, idx)
    pos = find(heap(:, 2) == idx, 1);
    if isempty(pos)
        heap = heap_push(heap, T_value, idx);
        return;
    end

    heap(pos, 1) = T_value;
    while pos > 1
        parent = floor(pos / 2);
        if heap(parent, 1) <= heap(pos, 1)
            break;
        end
        tmp = heap(parent, :);
        heap(parent, :) = heap(pos, :);
        heap(pos, :) = tmp;
        pos = parent;
    end
end
