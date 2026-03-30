function sensitivity = compute_sensitivity_adjoint(nelx, nely, U, theta_e, E_L, E_T, nu_LT, G_LT, t, dx, dy, arg12, arg13)
    % 纤维角度对应柔度 C = U'F 的伴随灵敏度分析
    % 对于线弹性平衡系统且载荷与边界条件固定，有
    %   dC/dtheta_e = -Ue' * (dKe/dtheta_e) * Ue
    % 这里使用与 FE 装配完全一致的单元刚度离散，避免目标与灵敏度使用不同离散。

    material_mask = true(nely, nelx);
    normalize_flag = true;

    if nargin >= 12 && ~isempty(arg12)
        if isscalar(arg12)
            normalize_flag = logical(arg12);
        elseif isequal(size(arg12), [nely, nelx])
            material_mask = logical(arg12);
        else
            error('arg12应为标量normalize_flag或[%d,%d]的material_mask。', nely, nelx);
        end
    end

    if nargin >= 13 && ~isempty(arg13)
        if ~isscalar(arg13)
            error('arg13应为标量normalize_flag。');
        end
        normalize_flag = logical(arg13);
    end

    sensitivity = zeros(nely, nelx);

    for ely = 1:nely
        for elx = 1:nelx
            if ~material_mask(ely, elx)
                sensitivity(ely, elx) = 0;
                continue;
            end
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
            [~, dKe_dtheta] = element_stiffness(theta, E_L, E_T, nu_LT, G_LT, t, dx, dy);
            sensitivity(ely, elx) = -full(Ue' * dKe_dtheta * Ue);
        end
    end

    if normalize_flag
        max_sensitivity = max(abs(sensitivity(:)));
        if max_sensitivity > 1e-10
            sensitivity = sensitivity / max_sensitivity;
        end
    end
end

