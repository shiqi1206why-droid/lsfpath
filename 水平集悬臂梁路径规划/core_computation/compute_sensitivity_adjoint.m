function sensitivity = compute_sensitivity_adjoint(nelx, nely, U, theta_e, E_L, E_T, nu_LT, G_LT, t, dx, dy, normalize_flag)
    % 纤维角度对应变能的伴随灵敏度分析

    if nargin < 12 || isempty(normalize_flag)
        normalize_flag = true;
    end

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

    if normalize_flag
        max_sensitivity = max(abs(sensitivity(:)));
        if max_sensitivity > 1e-10
            sensitivity = sensitivity / max_sensitivity;
        end
    end
end

