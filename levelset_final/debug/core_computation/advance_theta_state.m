function [theta_next, theta_target] = advance_theta_state(lsf, theta_prev, delta_theta_max, dx, dy, material_mask_core, smooth_eta, smooth_iterations)
    % 统一角度状态推进算子：
    % 1) theta_target = smooth(theta(lsf))
    % 2) theta_next   = limit(theta_prev -> theta_target, delta_theta_max)
    % 这样可保证当lsf不变时，theta状态不发生额外漂移

    if nargin < 6 || isempty(material_mask_core)
        material_mask_core = [];
    end
    if nargin < 8 && (isempty(material_mask_core) || isscalar(material_mask_core))
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

    theta_raw = compute_fiber_angles_from_lsf(lsf, dx, dy, material_mask_core);
    theta_for_smooth = theta_raw;
    if ~isempty(material_mask_core)
        theta_for_smooth(~material_mask_core) = 0;
    end
    z = angle_smooth_vectorized(theta_for_smooth, smooth_eta, smooth_iterations);
    theta_target = mod(0.5 * angle(z), pi);
    if ~isempty(material_mask_core)
        theta_target(~material_mask_core) = NaN;
    end

    if nargin < 2 || isempty(theta_prev)
        theta_next = theta_target;
        return;
    end

    if ~isempty(material_mask_core)
        theta_prev(~material_mask_core) = NaN;
    end

    theta_diff = theta_target - theta_prev;
    theta_diff = atan2(sin(theta_diff), cos(theta_diff));
    theta_step = sign(theta_diff) .* min(abs(theta_diff), delta_theta_max);
    theta_next = mod(theta_prev + theta_step, pi);
    if ~isempty(material_mask_core)
        theta_next(~material_mask_core) = NaN;
    end
end
