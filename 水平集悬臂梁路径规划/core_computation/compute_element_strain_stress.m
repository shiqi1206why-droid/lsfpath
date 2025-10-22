function [strain, stress, C, S, dC_dtheta] = compute_element_strain_stress(theta, Ue, E_L, E_T, nu_LT, G_LT, t, dx, dy)
    % 计算单元的应变、应力、刚度及其导数

    nu_TL = nu_LT * E_T / E_L;
    Q11 = E_L / (1 - nu_LT * nu_TL);
    Q22 = E_T / (1 - nu_LT * nu_TL);
    Q12 = nu_LT * E_T / (1 - nu_LT * nu_TL);
    Q66 = G_LT;

    Q = [Q11, Q12, 0;
         Q12, Q22, 0;
         0,   0,   Q66];

    c = cos(theta);
    s = sin(theta);

    T = [c^2,  s^2,   2*s*c;
         s^2,  c^2,  -2*s*c;
         -s*c, s*c,  c^2-s^2];

    dT_dtheta = [-2*c*s,  2*c*s,   2*(c^2-s^2);
                  2*c*s, -2*c*s,  -2*(c^2-s^2);
                 -c^2+s^2, c^2-s^2, -4*c*s];

    C = T' * Q * T;
    dC_dtheta = dT_dtheta' * Q * T + T' * Q * dT_dtheta;
    
    % === 本构矩阵求逆稳健化：防止奇异矩阵 ===
    % 参考：代码修改交流.md (Codex方案 2025-10-17)
    rc = rcond(C);
    if ~isfinite(rc) || rc < 1e-12
        % 轻微正则化避免奇异性
        C = C + 1e-9 * eye(3);
    end
    S = C \ eye(3);  % 用线性求解代替inv，更稳定

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

