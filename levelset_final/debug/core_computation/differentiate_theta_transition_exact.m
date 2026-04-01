function [pullback, diagnostics] = differentiate_theta_transition_exact(dC_dtheta_next, transition_cache, diff_opts)
%DIFFERENTIATE_THETA_TRANSITION_EXACT Reverse the shared theta transition discretization.

    if nargin < 3 || isempty(diff_opts)
        diff_opts = struct();
    end
    if ~isfield(diff_opts, 'limiter_mode') || isempty(diff_opts.limiter_mode)
        diff_opts.limiter_mode = 'hard';
    end
    if ~isfield(diff_opts, 'soft_limiter_beta') || isempty(diff_opts.soft_limiter_beta)
        diff_opts.soft_limiter_beta = 20.0;
    end
    if ~isfield(diff_opts, 'grad_floor') || isempty(diff_opts.grad_floor)
        diff_opts.grad_floor = 0.05;
    end

    material_mask_core = transition_cache.material_mask_core;
    dC_dtheta_next = sanitize_theta_adjoint(dC_dtheta_next, material_mask_core);

    limiter_weight = compute_limiter_pullback_weight(transition_cache, diff_opts);
    dtheta_target = dC_dtheta_next .* limiter_weight;
    dtheta_target(transition_cache.wrap_discontinuity_mask) = 0;
    if ~transition_cache.theta_prev_supplied
        dtheta_target = dC_dtheta_next;
    end

    [dz_real, dz_imag] = theta_target_pullback(dtheta_target, transition_cache.z_smooth, material_mask_core);
    [dtheta_raw, smooth_diag] = smooth_chain_pullback(dz_real, dz_imag, transition_cache.smooth_cache);
    dtheta_raw(~material_mask_core) = 0;

    [top_left_contrib, bottom_left_contrib, bottom_right_contrib, top_right_contrib, theta_raw_diag] = ...
        theta_raw_to_phi_pullback(dtheta_raw, transition_cache.angle_cache, diff_opts.grad_floor);

    pullback = struct();
    pullback.mode = 'exact';
    pullback.top_left_contrib = top_left_contrib;
    pullback.bottom_left_contrib = bottom_left_contrib;
    pullback.bottom_right_contrib = bottom_right_contrib;
    pullback.top_right_contrib = top_right_contrib;
    pullback.material_mask_core = material_mask_core;

    diagnostics = struct();
    diagnostics.limiter_mode = char(diff_opts.limiter_mode);
    diagnostics.saturation_ratio = safe_mask_fraction(transition_cache.limiter_saturated_mask, material_mask_core);
    diagnostics.degenerate_ratio = safe_mask_fraction(transition_cache.degenerate_grad_mask, material_mask_core);
    diagnostics.wrap_discontinuity_ratio = safe_mask_fraction(transition_cache.wrap_discontinuity_mask, material_mask_core);
    diagnostics.zero_gradient_due_to_limiter_ratio = safe_mask_fraction( ...
        transition_cache.limiter_saturated_mask, material_mask_core);
    diagnostics.smooth_clamp_ratio = smooth_diag.clamp_ratio;
    diagnostics.theta_raw_grad_floor = diff_opts.grad_floor;
    diagnostics.theta_raw_guard_ratio = theta_raw_diag.guard_ratio;
end

function dtheta = sanitize_theta_adjoint(dtheta, material_mask_core)
    dtheta = real(dtheta);
    dtheta(~isfinite(dtheta)) = 0;
    if ~isempty(material_mask_core)
        dtheta(~material_mask_core) = 0;
    end
end

function limiter_weight = compute_limiter_pullback_weight(transition_cache, diff_opts)
    switch lower(string(diff_opts.limiter_mode))
        case "soft_experiment"
            beta = diff_opts.soft_limiter_beta;
            limiter_weight = 1 ./ (1 + exp(beta * ...
                (abs(transition_cache.theta_diff_wrapped) - transition_cache.delta_theta_max)));
            limiter_weight(transition_cache.limiter_kink_mask) = 0.5;
        otherwise
            limiter_weight = double(transition_cache.limiter_linear_mask);
            limiter_weight(transition_cache.limiter_saturated_mask) = 0;
            limiter_weight(transition_cache.limiter_kink_mask) = 0;
    end
    limiter_weight(~transition_cache.material_mask_core) = 0;
end

function [dz_real, dz_imag] = theta_target_pullback(dtheta_target, z_smooth, material_mask_core)
    x = real(z_smooth);
    y = imag(z_smooth);
    denom = max(x .* x + y .* y, 1e-12);
    dz_real = dtheta_target .* (-0.5 * y ./ denom);
    dz_imag = dtheta_target .* (0.5 * x ./ denom);
    dz_real(~material_mask_core) = 0;
    dz_imag(~material_mask_core) = 0;
end

function [dtheta_raw, diagnostics] = smooth_chain_pullback(dz_real, dz_imag, smooth_cache)
    laplacian_kernel = [0, 1, 0; 1, -4, 1; 0, 1, 0];
    g_real = dz_real;
    g_imag = dz_imag;
    clamp_count = 0;
    total_count = 0;

    for iter = numel(smooth_cache.iterations):-1:1
        iter_cache = smooth_cache.iterations(iter);
        clamp_count = clamp_count + nnz(iter_cache.clamped_mask);
        total_count = total_count + numel(iter_cache.clamped_mask);

        [g_real, g_imag] = normalize_adjoint(g_real, g_imag, iter_cache);
        g_real = apply_neumann_boundary_adjoint(g_real);
        g_imag = apply_neumann_boundary_adjoint(g_imag);
        g_real = g_real + smooth_cache.eta * conv2(g_real, rot90(laplacian_kernel, 2), 'same');
        g_imag = g_imag + smooth_cache.eta * conv2(g_imag, rot90(laplacian_kernel, 2), 'same');
    end

    theta_in = smooth_cache.theta_input;
    dtheta_raw = g_real .* (-2 * sin(2 * theta_in)) + g_imag .* (2 * cos(2 * theta_in));
    dtheta_raw(~isfinite(dtheta_raw)) = 0;

    diagnostics = struct();
    diagnostics.clamp_ratio = safe_fraction(clamp_count, total_count);
end

function [gin_real, gin_imag] = normalize_adjoint(gout_real, gout_imag, iter_cache)
    x = real(iter_cache.z_after_bc);
    y = imag(iter_cache.z_after_bc);
    r = iter_cache.clamped_norm;
    clamp_mask = iter_cache.clamped_mask;

    dot_term = x .* gout_real + y .* gout_imag;
    gin_real = gout_real ./ r - x .* dot_term ./ (r .^ 3);
    gin_imag = gout_imag ./ r - y .* dot_term ./ (r .^ 3);

    if any(clamp_mask(:))
        gin_real(clamp_mask) = gout_real(clamp_mask) ./ r(clamp_mask);
        gin_imag(clamp_mask) = gout_imag(clamp_mask) ./ r(clamp_mask);
    end

    gin_real(~isfinite(gin_real)) = 0;
    gin_imag(~isfinite(gin_imag)) = 0;
end

function field = apply_neumann_boundary_adjoint(field)
    field(:, end-1) = field(:, end-1) + field(:, end);
    field(:, end) = 0;

    field(:, 2) = field(:, 2) + field(:, 1);
    field(:, 1) = 0;

    field(end-1, :) = field(end-1, :) + field(end, :);
    field(end, :) = 0;

    field(2, :) = field(2, :) + field(1, :);
    field(1, :) = 0;
end

function [top_left_contrib, bottom_left_contrib, bottom_right_contrib, top_right_contrib, diagnostics] = ...
    theta_raw_to_phi_pullback(dtheta_raw, angle_cache, grad_floor)
    gx = angle_cache.dphi_dx;
    gy = angle_cache.dphi_dy;
    degenerate_grad_mask = angle_cache.degenerate_grad_mask;
    material_mask_core = angle_cache.material_mask_core;
    dx = angle_cache.dx;
    dy = angle_cache.dy;
    dN_dx = [-0.5 / dx, 0.5 / dx, 0.5 / dx, -0.5 / dx];
    dN_dy = [-0.5 / dy, -0.5 / dy, 0.5 / dy, 0.5 / dy];

    grad_sq = gx .* gx + gy .* gy;
    denom = max(grad_sq, 1e-12);
    gx_coeff = dtheta_raw .* (-gy ./ denom);
    gy_coeff = dtheta_raw .* (gx ./ denom);

    guard_mask = degenerate_grad_mask | (grad_sq < grad_floor * grad_floor);
    invalid_mask = guard_mask | ~material_mask_core | ~isfinite(gx_coeff) | ~isfinite(gy_coeff);
    gx_coeff(invalid_mask) = 0;
    gy_coeff(invalid_mask) = 0;

    top_left_contrib = gx_coeff * dN_dx(1) + gy_coeff * dN_dy(1);
    bottom_left_contrib = gx_coeff * dN_dx(2) + gy_coeff * dN_dy(2);
    bottom_right_contrib = gx_coeff * dN_dx(3) + gy_coeff * dN_dy(3);
    top_right_contrib = gx_coeff * dN_dx(4) + gy_coeff * dN_dy(4);

    diagnostics = struct();
    diagnostics.guard_ratio = safe_mask_fraction(guard_mask, material_mask_core);
end

function ratio = safe_mask_fraction(mask, base_mask)
    mask = logical(mask);
    base_mask = logical(base_mask);
    denom = nnz(base_mask);
    if denom == 0
        ratio = 0;
    else
        ratio = nnz(mask & base_mask) / denom;
    end
end

function value = safe_fraction(num, den)
    if den <= 0
        value = 0;
    else
        value = num / den;
    end
end
