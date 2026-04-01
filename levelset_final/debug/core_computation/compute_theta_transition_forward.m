function transition = compute_theta_transition_forward(lsf, theta_prev, delta_theta_max, dx, dy, material_mask_core, smooth_eta, smooth_iterations)
%COMPUTE_THETA_TRANSITION_FORWARD Shared forward core for phi -> theta transition.

    if nargin < 6 || isempty(material_mask_core)
        material_mask_core = [];
    end
    if nargin < 8 && (isempty(material_mask_core) || isscalar(material_mask_core))
        smooth_iterations = smooth_eta;
        smooth_eta = material_mask_core;
        material_mask_core = [];
    elseif isscalar(material_mask_core) && (nargin < 8 || isempty(smooth_iterations))
        smooth_iterations = smooth_eta;
        smooth_eta = material_mask_core;
        material_mask_core = [];
    elseif isscalar(material_mask_core) && (nargin < 7 || isempty(smooth_eta))
        smooth_eta = material_mask_core;
        material_mask_core = [];
    end
    if nargin < 7 || isempty(smooth_eta)
        smooth_eta = 0.10;
    end
    if nargin < 8 || isempty(smooth_iterations)
        smooth_iterations = 2;
    end

    [theta_raw, angle_cache] = compute_fiber_angles_from_lsf(lsf, dx, dy, material_mask_core);
    theta_for_smooth = theta_raw;
    if ~isempty(material_mask_core)
        theta_for_smooth(~material_mask_core) = 0;
    end

    [z_smooth, smooth_cache] = angle_smooth_vectorized(theta_for_smooth, smooth_eta, smooth_iterations);
    theta_target = mod(0.5 * angle(z_smooth), pi);
    if ~isempty(material_mask_core)
        theta_target(~material_mask_core) = NaN;
    end

    theta_prev_supplied = ~(nargin < 2 || isempty(theta_prev));
    if theta_prev_supplied
        theta_prev_used = theta_prev;
        if ~isempty(material_mask_core)
            theta_prev_used(~material_mask_core) = NaN;
        end
        theta_diff_raw = theta_target - theta_prev_used;
        theta_diff_wrapped = atan2(sin(theta_diff_raw), cos(theta_diff_raw));
        theta_step = sign(theta_diff_wrapped) .* min(abs(theta_diff_wrapped), delta_theta_max);
        theta_next = mod(theta_prev_used + theta_step, pi);
        limiter_linear_mask = abs(theta_diff_wrapped) < (delta_theta_max - 1e-12);
        limiter_saturated_mask = abs(theta_diff_wrapped) > (delta_theta_max + 1e-12);
        limiter_kink_mask = ~(limiter_linear_mask | limiter_saturated_mask);
        wrap_discontinuity_mask = abs(abs(theta_diff_wrapped) - pi) <= 1e-12;
    else
        theta_prev_used = theta_target;
        theta_diff_raw = zeros(size(theta_target));
        theta_diff_wrapped = zeros(size(theta_target));
        theta_step = zeros(size(theta_target));
        theta_next = theta_target;
        limiter_linear_mask = true(size(theta_target));
        limiter_saturated_mask = false(size(theta_target));
        limiter_kink_mask = false(size(theta_target));
        wrap_discontinuity_mask = false(size(theta_target));
    end

    if ~isempty(material_mask_core)
        theta_prev_used(~material_mask_core) = NaN;
        theta_next(~material_mask_core) = NaN;
        limiter_linear_mask(~material_mask_core) = false;
        limiter_saturated_mask(~material_mask_core) = false;
        limiter_kink_mask(~material_mask_core) = false;
        wrap_discontinuity_mask(~material_mask_core) = false;
    end

    transition = struct();
    transition.theta_raw = theta_raw;
    transition.theta_target = theta_target;
    transition.theta_next = theta_next;
    transition.cache = struct( ...
        'phi_base', lsf, ...
        'theta_prev', theta_prev_used, ...
        'theta_prev_supplied', theta_prev_supplied, ...
        'theta_raw', theta_raw, ...
        'theta_for_smooth', theta_for_smooth, ...
        'z_smooth', z_smooth, ...
        'smooth_cache', smooth_cache, ...
        'theta_target', theta_target, ...
        'theta_diff_raw', theta_diff_raw, ...
        'theta_diff_wrapped', theta_diff_wrapped, ...
        'theta_step', theta_step, ...
        'theta_next', theta_next, ...
        'limiter_linear_mask', limiter_linear_mask, ...
        'limiter_saturated_mask', limiter_saturated_mask, ...
        'limiter_kink_mask', limiter_kink_mask, ...
        'wrap_discontinuity_mask', wrap_discontinuity_mask, ...
        'degenerate_grad_mask', angle_cache.degenerate_grad_mask, ...
        'angle_cache', angle_cache, ...
        'material_mask_core', angle_cache.material_mask_core, ...
        'delta_theta_max', delta_theta_max, ...
        'dx', dx, ...
        'dy', dy, ...
        'smooth_eta', smooth_eta, ...
        'smooth_iterations', smooth_iterations);
end
