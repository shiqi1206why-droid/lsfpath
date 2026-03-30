function candidate = evaluate_candidate_state(lsf_trial, theta_prev, delta_theta_max, dx, dy, ...
    nelx, nely, material_mask, E_L, E_T, nu_LT, G_LT, thickness, F_mag, smooth_eta, smooth_iterations)
% 评估从当前theta状态出发、在候选lsf上形成的一步候选状态

    if nargin < 15 || isempty(smooth_eta)
        smooth_eta = 0.10;
    end
    if nargin < 16 || isempty(smooth_iterations)
        smooth_iterations = 2;
    end

    [theta_trial, theta_target] = advance_theta_state(lsf_trial, theta_prev, delta_theta_max, ...
        dx, dy, material_mask, smooth_eta, smooth_iterations);
    state = evaluate_state_with_theta(lsf_trial, theta_trial, nelx, nely, material_mask, ...
        E_L, E_T, nu_LT, G_LT, thickness, F_mag, dx, dy);

    candidate = state;
    candidate.theta_target = theta_target;
end
