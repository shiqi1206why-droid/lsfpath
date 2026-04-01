function [theta_next, theta_target, theta_transition_cache] = advance_theta_state(lsf, theta_prev, delta_theta_max, dx, dy, material_mask_core, smooth_eta, smooth_iterations)
    % 统一角度状态推进算子：
    % 1) theta_target = smooth(theta(lsf))
    % 2) theta_next   = limit(theta_prev -> theta_target, delta_theta_max)
    % 这样可保证当lsf不变时，theta状态不发生额外漂移

    transition = compute_theta_transition_forward(lsf, theta_prev, delta_theta_max, dx, dy, ...
        material_mask_core, smooth_eta, smooth_iterations);
    theta_next = transition.theta_next;
    theta_target = transition.theta_target;
    theta_transition_cache = transition.cache;
end
