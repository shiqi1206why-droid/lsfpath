function trial_compliance = evaluate_candidate_compliance(lsf_trial, theta_ref, delta_theta_max, dx, dy, ...
    nelx, nely, material_mask, E_L, E_T, nu_LT, G_LT, thickness, F_mag, smooth_eta, smooth_iterations)
% 计算候选水平集更新在下一步对应的柔度（用于单步接受/拒绝）

    try
        state = evaluate_candidate_state(lsf_trial, theta_ref, delta_theta_max, dx, dy, ...
            nelx, nely, material_mask, E_L, E_T, nu_LT, G_LT, thickness, F_mag, ...
            smooth_eta, smooth_iterations);
        trial_compliance = state.compliance;
        if ~isfinite(trial_compliance)
            trial_compliance = inf;
        end
    catch
        trial_compliance = inf;
    end
end
