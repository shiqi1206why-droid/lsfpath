function node_sensitivity = aggregate_node_sensitivity(element_sensitivity, ~, lsf, nelx, nely, dx, dy, primary_update_mask)
    % 通过链式法则将单元灵敏度汇总到水平集节点

    global DIAG;
    node_sensitivity = zeros(size(lsf));

    if nargin < 8 || isempty(primary_update_mask)
        primary_update_mask = true(size(lsf));
    end
    if isequal(size(primary_update_mask), [nely, nelx])
        mask_full = false(size(lsf));
        mask_full(2:end-1, 2:end-1) = logical(primary_update_mask);
        primary_update_mask = mask_full;
    end
    if ~isequal(size(primary_update_mask), size(lsf))
        error('primary_update_mask尺寸错误：期望[%d,%d]或[%d,%d]，实际[%d,%d]。', ...
            size(lsf,1), size(lsf,2), nely, nelx, size(primary_update_mask,1), size(primary_update_mask,2));
    end
    primary_update_mask = logical(primary_update_mask);

    dN_dxi = [-0.25,  0.25,  0.25, -0.25];
    dN_deta = [-0.25, -0.25,  0.25,  0.25];
    dN_dx = dN_dxi * (2/dx);
    dN_dy = dN_deta * (2/dy);
    grad_threshold = 0.05;
    eps_denom = 1e-12;            % 正则化分母的下限

    for ely = 1:nely
        for elx = 1:nelx
            base_i = ely + 1;
            base_j = elx + 1;

            node_coords = [base_i,   base_j;
                           base_i+1, base_j;
                           base_i+1, base_j+1;
                           base_i,   base_j+1];
            node_active = false(1, 4);
            for k = 1:4
                node_active(k) = primary_update_mask(node_coords(k,1), node_coords(k,2));
            end
            if ~any(node_active)
                continue;
            end

            phi_nodes = zeros(1,4);
            for k = 1:4
                phi_nodes(k) = lsf(node_coords(k,1), node_coords(k,2));
            end

            dphi_dx = sum(dN_dx .* phi_nodes);
            dphi_dy = sum(dN_dy .* phi_nodes);
            sum_dx_phi = dphi_dx;
            sum_dy_phi = dphi_dy;
            grad_sq = dphi_dx^2 + dphi_dy^2;

            if ~isempty(DIAG)
                DIAG.theta_grad_samples = DIAG.theta_grad_samples + 1;
            end

            if grad_sq < grad_threshold^2
                if ~isempty(DIAG)
                    DIAG.theta_grad_small = DIAG.theta_grad_small + 1;
                end
                continue;
            end

            coeff = element_sensitivity(ely, elx);
            if coeff == 0
                continue;
            end

            for k = 1:4
                if ~node_active(k)
                    continue;
                end
                phi_i = phi_nodes(k);
                Pi = sum_dx_phi - dN_dx(k) * phi_i;
                Qi = sum_dy_phi - dN_dy(k) * phi_i;
                A = dN_dx(k) * phi_i + Pi;
                B = dN_dy(k) * phi_i + Qi;
                denom = max(A * A + B * B, eps_denom);
                numerator = dN_dy(k) * A - dN_dx(k) * B;
                if ~isfinite(numerator)
                    if ~isempty(DIAG)
                        DIAG.den_small = DIAG.den_small + 1;
                    end
                    continue;
                end
                contrib = coeff * (numerator / denom);
                if ~isfinite(contrib)
                    warning('aggregate_node_sensitivity: 非有限贡献 (ely=%d, elx=%d, node=%d)', ely, elx, k);
                    continue;
                end
                node_sensitivity(node_coords(k,1), node_coords(k,2)) = node_sensitivity(node_coords(k,1), node_coords(k,2)) + contrib;
                if ~isempty(DIAG) && contrib ~= 0
                    DIAG.contrib_nonzero = DIAG.contrib_nonzero + 1;
                end
            end

            if ~isempty(DIAG)
                DIAG.contrib_total = DIAG.contrib_total + 4;
            end
        end
    end

    % 强制掩膜外节点为零，避免任何数值泄漏
    node_sensitivity(~primary_update_mask) = 0;
end

