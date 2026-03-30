function state = evaluate_state_with_theta(lsf_state, theta_e, nelx, nely, material_mask, ...
    E_L, E_T, nu_LT, G_LT, thickness, F_mag, dx, dy)
% 评估给定(lsf, theta)状态对应的FE响应与目标量

    [U, K, F] = FE_analysis_cantilever(nelx, nely, theta_e, ...
        E_L, E_T, nu_LT, G_LT, thickness, F_mag, dx, dy, material_mask);

    compliance = full(U' * F);
    strain_energy = compute_strain_energy(nelx, nely, U, theta_e, ...
        E_L, E_T, nu_LT, G_LT, thickness, dx, dy, material_mask);
    FCS = compute_fiber_continuity(theta_e, material_mask);

    state = struct();
    state.lsf = lsf_state;
    state.theta = theta_e;
    state.U = U;
    state.K = K;
    state.F = F;
    state.compliance = compliance;
    state.strain_energy = strain_energy;
    state.FCS = FCS;
end
