function sensitivity = compute_sensitivity_adjoint(nelx, nely, U, theta_e, E_L, E_T, nu_LT, G_LT, t, dx, dy)
    % 纤维角度对应变能的伴随灵敏度分析

    sensitivity = zeros(nely, nelx);
    alpha = 0.5;

    for ely = 1:nely
        for elx = 1:nelx
            n1 = (nely+1)*(elx-1) + ely;
            n2 = (nely+1)*elx + ely;
            n3 = n2 + 1;
            n4 = n1 + 1;
            nodes = [n1, n2, n3, n4];
            edof = [];
            for n = nodes
                edof = [edof, 2*n-1, 2*n]; %#ok<AGROW>
            end
            Ue = U(edof);

            theta = theta_e(ely, elx);
            [strain, stress, C, S, dC_dtheta] = compute_element_strain_stress(theta, Ue, E_L, E_T, nu_LT, G_LT, t, dx, dy);
            dS_dtheta = -S * dC_dtheta * S;

            term1 = (-1 + 2*alpha - alpha^2) * (strain' * dC_dtheta * strain);
            term2 = alpha^2 * (stress' * dS_dtheta * stress);
            term3 = (2*alpha^2 - 2*alpha) * (strain' * dC_dtheta * stress);

            sensitivity(ely, elx) = term1 + term2 + term3;
        end
    end

    % 保留真实灵敏度量级，供时间步长公式使用（修改时间：2025-10-30）
    % 原归一化逻辑已注释，以符合 Δt = Δθ_max / max|∂E/∂φ| 公式要求
    % max_sensitivity = max(abs(sensitivity(:)));
    % if max_sensitivity > 1e-10
    %     sensitivity = sensitivity / max_sensitivity;
    % end
end

