function [strain, stress, C, S, dC_dtheta] = compute_element_strain_stress(theta, Ue, E_L, E_T, nu_LT, G_LT, t, dx, dy)
    % 计算单元的应变、应力、刚度及其导数

    [C, dC_dtheta, S] = orthotropic_constitutive_matrix(theta, E_L, E_T, nu_LT, G_LT);

    dN_dxi = 0.25 * [-1, 1, 1, -1];
    dN_deta = 0.25 * [-1, -1, 1, 1];
    J = [dx/2, 0; 0, dy/2];
    dN_dxy = J \ [dN_dxi; dN_deta];

    B = zeros(3, 8);
    for k = 1:4
        B(1, 2*k-1) = dN_dxy(1, k);
        B(2, 2*k) = dN_dxy(2, k);
        B(3, 2*k-1) = dN_dxy(2, k);
        B(3, 2*k) = dN_dxy(1, k);
    end

    strain = B * Ue;
    stress = C * strain;
end

