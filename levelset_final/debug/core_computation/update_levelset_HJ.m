function [lsf_new, diagnostics] = update_levelset_HJ(lsf, velocity, dt, dx, dy, active_mask, hj_opts)
    % 采用二阶 ENO/Godunov + SSPRK2 求解 Hamilton-Jacobi 方程
    % 当局部 stencil 不足、接近活动域边界或梯度退化时，
    % 显式回退到一阶单调迎风格式。

    if nargin < 6 || isempty(active_mask)
        active_mask = true(size(lsf));
    else
        active_mask = normalize_mask_to_lsf_grid(active_mask, size(lsf), 'active_mask');
    end
    if nargin < 7 || isempty(hj_opts)
        hj_opts = struct();
    end

    hj_opts = apply_hj_defaults(hj_opts);
    diagnostics = init_diagnostics(hj_opts);

    if strcmpi(hj_opts.time_integrator, 'ssprk2') && hj_opts.advection_order >= 2
        [rhs_1, diag_1] = compute_hj_rhs(lsf, velocity, dx, dy, active_mask, hj_opts);
        phi_stage_1 = lsf;
        phi_stage_1(active_mask) = lsf(active_mask) + dt * rhs_1(active_mask);
        phi_stage_1 = apply_neumann_boundary(phi_stage_1);

        [rhs_2, diag_2] = compute_hj_rhs(phi_stage_1, velocity, dx, dy, active_mask, hj_opts);
        lsf_new = lsf;
        lsf_new(active_mask) = 0.5 * lsf(active_mask) + ...
            0.5 * (phi_stage_1(active_mask) + dt * rhs_2(active_mask));
        lsf_new = apply_neumann_boundary(lsf_new);

        diagnostics.stage_1 = diag_1;
        diagnostics.stage_2 = diag_2;
        diagnostics.second_order_points = diag_1.second_order_points + diag_2.second_order_points;
        diagnostics.first_order_points = diag_1.first_order_points + diag_2.first_order_points;
        diagnostics.fallback_insufficient_stencil = diag_1.fallback_insufficient_stencil + diag_2.fallback_insufficient_stencil;
        diagnostics.fallback_degenerate_gradient = diag_1.fallback_degenerate_gradient + diag_2.fallback_degenerate_gradient;
        diagnostics.fallback_oscillation_guard = diag_1.fallback_oscillation_guard + diag_2.fallback_oscillation_guard;
        diagnostics.frozen_incomplete_godunov_count = diag_1.frozen_incomplete_godunov_count + diag_2.frozen_incomplete_godunov_count;
        diagnostics.used_first_order_complete_count = diag_1.used_first_order_complete_count + diag_2.used_first_order_complete_count;
        diagnostics.used_second_order_complete_count = diag_1.used_second_order_complete_count + diag_2.used_second_order_complete_count;
    else
        [rhs_1, diag_1] = compute_hj_rhs(lsf, velocity, dx, dy, active_mask, hj_opts);
        lsf_new = lsf;
        lsf_new(active_mask) = lsf(active_mask) + dt * rhs_1(active_mask);
        lsf_new = apply_neumann_boundary(lsf_new);
        diagnostics.stage_1 = diag_1;
        diagnostics.stage_2 = struct([]);
        diagnostics.second_order_points = diag_1.second_order_points;
        diagnostics.first_order_points = diag_1.first_order_points;
        diagnostics.fallback_insufficient_stencil = diag_1.fallback_insufficient_stencil;
        diagnostics.fallback_degenerate_gradient = diag_1.fallback_degenerate_gradient;
        diagnostics.fallback_oscillation_guard = diag_1.fallback_oscillation_guard;
        diagnostics.frozen_incomplete_godunov_count = diag_1.frozen_incomplete_godunov_count;
        diagnostics.used_first_order_complete_count = diag_1.used_first_order_complete_count;
        diagnostics.used_second_order_complete_count = diag_1.used_second_order_complete_count;
    end

    diagnostics.second_order_count = diagnostics.second_order_points;
    diagnostics.first_order_count = diagnostics.first_order_points;
    diagnostics.fallback_count = diagnostics.fallback_insufficient_stencil + ...
        diagnostics.fallback_degenerate_gradient + diagnostics.fallback_oscillation_guard;
end

function hj_opts = apply_hj_defaults(hj_opts)
    if ~isfield(hj_opts, 'advection_order') || isempty(hj_opts.advection_order)
        hj_opts.advection_order = 2;
    end
    if ~isfield(hj_opts, 'time_integrator') || isempty(hj_opts.time_integrator)
        hj_opts.time_integrator = 'ssprk2';
    end
    if ~isfield(hj_opts, 'fallback_first_order') || isempty(hj_opts.fallback_first_order)
        hj_opts.fallback_first_order = true;
    end
    if ~isfield(hj_opts, 'gradient_tol') || isempty(hj_opts.gradient_tol)
        hj_opts.gradient_tol = 1e-12;
    end
    if ~isfield(hj_opts, 'oscillation_ratio_limit') || isempty(hj_opts.oscillation_ratio_limit)
        hj_opts.oscillation_ratio_limit = 4.0;
    end
    if ~isfield(hj_opts, 'eno_smoothness_factor') || isempty(hj_opts.eno_smoothness_factor)
        hj_opts.eno_smoothness_factor = 2.5;
    end
    if ~isfield(hj_opts, 'stencil_mask') || isempty(hj_opts.stencil_mask)
        hj_opts.stencil_mask = [];
    end
    if ~isfield(hj_opts, 'freeze_on_incomplete_godunov') || isempty(hj_opts.freeze_on_incomplete_godunov)
        hj_opts.freeze_on_incomplete_godunov = true;
    end
    if ~isfield(hj_opts, 'stencil_buffer_cells') || isempty(hj_opts.stencil_buffer_cells)
        hj_opts.stencil_buffer_cells = 2;
    end
    if ~isfield(hj_opts, 'rhs_mode') || isempty(hj_opts.rhs_mode)
        hj_opts.rhs_mode = 'legacy';
    end
end

function diagnostics = init_diagnostics(hj_opts)
    diagnostics = struct();
    diagnostics.scheme = 'godunov';
    diagnostics.advection_order = hj_opts.advection_order;
    diagnostics.time_integrator = hj_opts.time_integrator;
    diagnostics.rhs_mode = char(hj_opts.rhs_mode);
    diagnostics.second_order_points = 0;
    diagnostics.first_order_points = 0;
    diagnostics.fallback_insufficient_stencil = 0;
    diagnostics.fallback_degenerate_gradient = 0;
    diagnostics.fallback_oscillation_guard = 0;
    diagnostics.frozen_incomplete_godunov_count = 0;
    diagnostics.used_first_order_complete_count = 0;
    diagnostics.used_second_order_complete_count = 0;
    diagnostics.second_order_count = 0;
    diagnostics.first_order_count = 0;
    diagnostics.fallback_count = 0;
end

function [rhs, diagnostics] = compute_hj_rhs(phi, velocity, dx, dy, active_mask, hj_opts)
    if isempty(hj_opts.stencil_mask)
        stencil_mask = active_mask;
    else
        stencil_mask = normalize_mask_to_lsf_grid(hj_opts.stencil_mask, size(phi), 'stencil_mask');
    end

    rhs_mode = lower(string(hj_opts.rhs_mode));
    switch rhs_mode
        case "legacy"
            [rhs, diagnostics] = compute_hj_rhs_legacy(phi, velocity, dx, dy, active_mask, stencil_mask, hj_opts);
        case "indexed"
            [rhs, diagnostics] = compute_hj_rhs_indexed(phi, velocity, dx, dy, active_mask, stencil_mask, hj_opts);
        case "vectorized_first_order"
            if hj_opts.advection_order < 2
                [rhs, diagnostics] = compute_hj_rhs_vectorized_first_order(phi, velocity, dx, dy, active_mask, stencil_mask, hj_opts);
            else
                [rhs, diagnostics] = compute_hj_rhs_indexed(phi, velocity, dx, dy, active_mask, stencil_mask, hj_opts);
            end
        otherwise
            warning('update_levelset_HJ:UnknownRhsMode', ...
                '未知rhs_mode=%s，回退到legacy模式。', hj_opts.rhs_mode);
            [rhs, diagnostics] = compute_hj_rhs_legacy(phi, velocity, dx, dy, active_mask, stencil_mask, hj_opts);
    end
end

function [rhs, diagnostics] = compute_hj_rhs_legacy(phi, velocity, dx, dy, active_mask, stencil_mask, hj_opts)
    [ny, nx] = size(phi);
    rhs = zeros(size(phi));
    diagnostics = init_diagnostics(hj_opts);

    for i = 2:ny-1
        for j = 2:nx-1
            if ~active_mask(i, j)
                continue;
            end

            v = velocity(i, j);
            if ~isfinite(v) || abs(v) <= hj_opts.gradient_tol
                rhs(i, j) = 0;
                continue;
            end

            [dmx, dpx, dmy, dpy, deriv_diag] = directional_derivatives(phi, i, j, dx, dy, stencil_mask, hj_opts);
            diagnostics.second_order_points = diagnostics.second_order_points + deriv_diag.second_order_used;
            diagnostics.first_order_points = diagnostics.first_order_points + deriv_diag.first_order_used;
            diagnostics.fallback_insufficient_stencil = diagnostics.fallback_insufficient_stencil + deriv_diag.fallback_insufficient_stencil;
            diagnostics.fallback_degenerate_gradient = diagnostics.fallback_degenerate_gradient + deriv_diag.fallback_degenerate_gradient;
            diagnostics.fallback_oscillation_guard = diagnostics.fallback_oscillation_guard + deriv_diag.fallback_oscillation_guard;
            diagnostics.frozen_incomplete_godunov_count = diagnostics.frozen_incomplete_godunov_count + deriv_diag.frozen_incomplete_godunov;
            diagnostics.used_first_order_complete_count = diagnostics.used_first_order_complete_count + deriv_diag.first_order_complete_used;
            diagnostics.used_second_order_complete_count = diagnostics.used_second_order_complete_count + deriv_diag.second_order_complete_used;

            if deriv_diag.freeze_point
                rhs(i, j) = 0;
                continue;
            end

            if v > 0
                grad_godunov = sqrt(max(dmx, 0)^2 + min(dpx, 0)^2 + ...
                    max(dmy, 0)^2 + min(dpy, 0)^2);
            else
                grad_godunov = sqrt(min(dmx, 0)^2 + max(dpx, 0)^2 + ...
                    min(dmy, 0)^2 + max(dpy, 0)^2);
            end
            rhs(i, j) = -v * grad_godunov;
        end
    end

    diagnostics = finalize_rhs_diagnostics(diagnostics);
end

function [rhs, diagnostics] = compute_hj_rhs_indexed(phi, velocity, dx, dy, active_mask, stencil_mask, hj_opts)
    [ny, nx] = size(phi);
    rhs = zeros(size(phi));
    diagnostics = init_diagnostics(hj_opts);

    interior_mask = false(size(phi));
    interior_mask(2:ny-1, 2:nx-1) = true;
    process_mask = interior_mask & active_mask;
    process_mask = process_mask & isfinite(velocity) & (abs(velocity) > hj_opts.gradient_tol);
    linear_idx = find(process_mask);

    for k = 1:numel(linear_idx)
        [i, j] = ind2sub([ny, nx], linear_idx(k));
        v = velocity(i, j);

        [dmx, dpx, dmy, dpy, deriv_diag] = directional_derivatives(phi, i, j, dx, dy, stencil_mask, hj_opts);
        diagnostics.second_order_points = diagnostics.second_order_points + deriv_diag.second_order_used;
        diagnostics.first_order_points = diagnostics.first_order_points + deriv_diag.first_order_used;
        diagnostics.fallback_insufficient_stencil = diagnostics.fallback_insufficient_stencil + deriv_diag.fallback_insufficient_stencil;
        diagnostics.fallback_degenerate_gradient = diagnostics.fallback_degenerate_gradient + deriv_diag.fallback_degenerate_gradient;
        diagnostics.fallback_oscillation_guard = diagnostics.fallback_oscillation_guard + deriv_diag.fallback_oscillation_guard;
        diagnostics.frozen_incomplete_godunov_count = diagnostics.frozen_incomplete_godunov_count + deriv_diag.frozen_incomplete_godunov;
        diagnostics.used_first_order_complete_count = diagnostics.used_first_order_complete_count + deriv_diag.first_order_complete_used;
        diagnostics.used_second_order_complete_count = diagnostics.used_second_order_complete_count + deriv_diag.second_order_complete_used;

        if deriv_diag.freeze_point
            rhs(i, j) = 0;
            continue;
        end

        if v > 0
            grad_godunov = sqrt(max(dmx, 0)^2 + min(dpx, 0)^2 + ...
                max(dmy, 0)^2 + min(dpy, 0)^2);
        else
            grad_godunov = sqrt(min(dmx, 0)^2 + max(dpx, 0)^2 + ...
                min(dmy, 0)^2 + max(dpy, 0)^2);
        end
        rhs(i, j) = -v * grad_godunov;
    end

    diagnostics = finalize_rhs_diagnostics(diagnostics);
end

function [rhs, diagnostics] = compute_hj_rhs_vectorized_first_order(phi, velocity, dx, dy, active_mask, stencil_mask, hj_opts)
    [ny, nx] = size(phi);
    rhs = zeros(size(phi));
    diagnostics = init_diagnostics(hj_opts);

    interior_mask = false(size(phi));
    interior_mask(2:ny-1, 2:nx-1) = true;
    process_mask = interior_mask & active_mask;
    process_mask = process_mask & isfinite(velocity) & (abs(velocity) > hj_opts.gradient_tol);

    left_ok = false(size(phi));
    right_ok = false(size(phi));
    down_ok = false(size(phi));
    up_ok = false(size(phi));
    left_ok(:, 2:end) = stencil_mask(:, 1:end-1);
    right_ok(:, 1:end-1) = stencil_mask(:, 2:end);
    down_ok(2:end, :) = stencil_mask(1:end-1, :);
    up_ok(1:end-1, :) = stencil_mask(2:end, :);
    complete_first_order = left_ok & right_ok & down_ok & up_ok;

    incomplete_mask = process_mask & ~complete_first_order;
    diagnostics.fallback_insufficient_stencil = nnz(incomplete_mask);

    if hj_opts.freeze_on_incomplete_godunov
        diagnostics.frozen_incomplete_godunov_count = nnz(incomplete_mask);
        compute_mask = process_mask & complete_first_order;
    else
        compute_mask = process_mask;
    end

    if any(compute_mask(:))
        dmx = zeros(size(phi));
        dpx = zeros(size(phi));
        dmy = zeros(size(phi));
        dpy = zeros(size(phi));
        dmx(:, 2:end) = (phi(:, 2:end) - phi(:, 1:end-1)) / dx;
        dpx(:, 1:end-1) = (phi(:, 2:end) - phi(:, 1:end-1)) / dx;
        dmy(2:end, :) = (phi(2:end, :) - phi(1:end-1, :)) / dy;
        dpy(1:end-1, :) = (phi(2:end, :) - phi(1:end-1, :)) / dy;

        grad_sq = zeros(size(phi));
        pos_mask = compute_mask & (velocity > 0);
        neg_mask = compute_mask & ~pos_mask;

        grad_sq(pos_mask) = max(dmx(pos_mask), 0).^2 + min(dpx(pos_mask), 0).^2 + ...
            max(dmy(pos_mask), 0).^2 + min(dpy(pos_mask), 0).^2;
        grad_sq(neg_mask) = min(dmx(neg_mask), 0).^2 + max(dpx(neg_mask), 0).^2 + ...
            min(dmy(neg_mask), 0).^2 + max(dpy(neg_mask), 0).^2;
        rhs(compute_mask) = -velocity(compute_mask) .* sqrt(grad_sq(compute_mask));
    end

    first_order_complete_used_mask = compute_mask & complete_first_order;
    diagnostics.first_order_points = nnz(first_order_complete_used_mask);
    diagnostics.used_first_order_complete_count = diagnostics.first_order_points;
    diagnostics = finalize_rhs_diagnostics(diagnostics);
end

function diagnostics = finalize_rhs_diagnostics(diagnostics)
    diagnostics.second_order_count = diagnostics.second_order_points;
    diagnostics.first_order_count = diagnostics.first_order_points;
    diagnostics.fallback_count = diagnostics.fallback_insufficient_stencil + ...
        diagnostics.fallback_degenerate_gradient + diagnostics.fallback_oscillation_guard;
end

function [dmx, dpx, dmy, dpy, diag_info] = directional_derivatives(phi, i, j, dx, dy, stencil_mask, hj_opts)
    left_ok = stencil_mask(i, j-1);
    right_ok = stencil_mask(i, j+1);
    down_ok = stencil_mask(i-1, j);
    up_ok = stencil_mask(i+1, j);
    has_complete_first_order = left_ok && right_ok && down_ok && up_ok;

    dmx = 0;
    dpx = 0;
    dmy = 0;
    dpy = 0;

    diag_info = struct( ...
        'second_order_used', 0, ...
        'first_order_used', 0, ...
        'fallback_insufficient_stencil', 0, ...
        'fallback_degenerate_gradient', 0, ...
        'fallback_oscillation_guard', 0, ...
        'freeze_point', false, ...
        'frozen_incomplete_godunov', 0, ...
        'first_order_complete_used', 0, ...
        'second_order_complete_used', 0);

    if ~has_complete_first_order
        diag_info.fallback_insufficient_stencil = 1;
        if hj_opts.freeze_on_incomplete_godunov
            diag_info.freeze_point = true;
            diag_info.frozen_incomplete_godunov = 1;
            return;
        end
    end

    dmx_1 = (phi(i, j) - phi(i, j-1)) / dx;
    dpx_1 = (phi(i, j+1) - phi(i, j)) / dx;
    dmy_1 = (phi(i, j) - phi(i-1, j)) / dy;
    dpy_1 = (phi(i+1, j) - phi(i, j)) / dy;

    dmx = dmx_1;
    dpx = dpx_1;
    dmy = dmy_1;
    dpy = dpy_1;
    diag_info.first_order_used = 1;
    diag_info.first_order_complete_used = 1;

    if ~has_complete_first_order
        diag_info.first_order_complete_used = 0;
        diag_info.first_order_used = 0;
    end

    if hj_opts.advection_order < 2
        return;
    end

    base_grad = max(abs([dmx_1, dpx_1, dmy_1, dpy_1]));
    if ~isfinite(base_grad) || base_grad < hj_opts.gradient_tol
        diag_info.fallback_degenerate_gradient = 1;
        return;
    end

    [dmx_2, ok_bx, reason_bx] = select_backward_derivative(phi(i, :), j, dx, stencil_mask(i, :), hj_opts);
    [dpx_2, ok_fx, reason_fx] = select_forward_derivative(phi(i, :), j, dx, stencil_mask(i, :), hj_opts);
    [dmy_2, ok_by, reason_by] = select_backward_derivative(phi(:, j).', i, dy, stencil_mask(:, j).', hj_opts);
    [dpy_2, ok_fy, reason_fy] = select_forward_derivative(phi(:, j).', i, dy, stencil_mask(:, j).', hj_opts);

    reasons = {reason_bx, reason_fx, reason_by, reason_fy};
    if ~(ok_bx && ok_fx && ok_by && ok_fy)
        diag_info.fallback_insufficient_stencil = any(strcmp(reasons, 'stencil'));
        diag_info.fallback_oscillation_guard = any(strcmp(reasons, 'oscillation'));
        diag_info.fallback_degenerate_gradient = diag_info.fallback_degenerate_gradient + any(strcmp(reasons, 'degenerate'));
        return;
    end

    dmx = dmx_2;
    dpx = dpx_2;
    dmy = dmy_2;
    dpy = dpy_2;
    diag_info.second_order_used = 1;
    diag_info.first_order_used = 0;
    diag_info.first_order_complete_used = 0;
    diag_info.second_order_complete_used = 1;
end

function [deriv, ok, reason] = select_backward_derivative(line_values, idx, step, line_mask, hj_opts)
    deriv = (line_values(idx) - line_values(idx-1)) / step;
    ok = false;
    reason = 'stencil';

    can_upwind2 = idx >= 3 && line_mask(idx-1) && line_mask(idx-2);
    can_central2 = idx >= 2 && idx <= numel(line_values)-1 && line_mask(idx-1) && line_mask(idx+1);
    if ~(can_upwind2 || can_central2)
        return;
    end

    first_order = deriv;
    second_candidates = [];
    smoothness = [];

    if can_upwind2
        second_candidates(end+1) = (3 * line_values(idx) - 4 * line_values(idx-1) + line_values(idx-2)) / (2 * step); %#ok<AGROW>
        smoothness(end+1) = abs(line_values(idx) - 2 * line_values(idx-1) + line_values(idx-2)); %#ok<AGROW>
    end
    if can_central2
        second_candidates(end+1) = (line_values(idx+1) - line_values(idx-1)) / (2 * step); %#ok<AGROW>
        smoothness(end+1) = abs(line_values(idx+1) - 2 * line_values(idx) + line_values(idx-1)); %#ok<AGROW>
    end

    [~, best_idx] = min(smoothness);
    candidate = second_candidates(best_idx);
    if ~isfinite(candidate)
        reason = 'degenerate';
        return;
    end
    selected_smoothness = smoothness(best_idx);
    smoothness_limit = hj_opts.eno_smoothness_factor * max(abs(first_order) * step, hj_opts.gradient_tol * step);
    if selected_smoothness > smoothness_limit
        reason = 'oscillation';
        return;
    end
    if abs(candidate) > hj_opts.oscillation_ratio_limit * max(abs(first_order), hj_opts.gradient_tol)
        reason = 'oscillation';
        return;
    end

    deriv = candidate;
    ok = true;
    reason = '';
end

function [deriv, ok, reason] = select_forward_derivative(line_values, idx, step, line_mask, hj_opts)
    deriv = (line_values(idx+1) - line_values(idx)) / step;
    ok = false;
    reason = 'stencil';

    can_forward2 = idx <= numel(line_values)-2 && line_mask(idx+1) && line_mask(idx+2);
    can_central2 = idx >= 2 && idx <= numel(line_values)-1 && line_mask(idx-1) && line_mask(idx+1);
    if ~(can_forward2 || can_central2)
        return;
    end

    first_order = deriv;
    second_candidates = [];
    smoothness = [];

    if can_forward2
        second_candidates(end+1) = (-line_values(idx+2) + 4 * line_values(idx+1) - 3 * line_values(idx)) / (2 * step); %#ok<AGROW>
        smoothness(end+1) = abs(line_values(idx+2) - 2 * line_values(idx+1) + line_values(idx)); %#ok<AGROW>
    end
    if can_central2
        second_candidates(end+1) = (line_values(idx+1) - line_values(idx-1)) / (2 * step); %#ok<AGROW>
        smoothness(end+1) = abs(line_values(idx+1) - 2 * line_values(idx) + line_values(idx-1)); %#ok<AGROW>
    end

    [~, best_idx] = min(smoothness);
    candidate = second_candidates(best_idx);
    if ~isfinite(candidate)
        reason = 'degenerate';
        return;
    end
    selected_smoothness = smoothness(best_idx);
    smoothness_limit = hj_opts.eno_smoothness_factor * max(abs(first_order) * step, hj_opts.gradient_tol * step);
    if selected_smoothness > smoothness_limit
        reason = 'oscillation';
        return;
    end
    if abs(candidate) > hj_opts.oscillation_ratio_limit * max(abs(first_order), hj_opts.gradient_tol)
        reason = 'oscillation';
        return;
    end

    deriv = candidate;
    ok = true;
    reason = '';
end
