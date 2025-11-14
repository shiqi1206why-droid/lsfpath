function node_sensitivity = aggregate_node_sensitivity(element_sensitivity, theta_field, lsf, nelx, nely, dx, dy)
    % 通过链式法则将单元灵敏度汇总到水平集节点

    global DIAG;
    node_sensitivity = zeros(size(lsf));

    dN_dxi = [-0.25,  0.25,  0.25, -0.25];
    dN_deta = [-0.25, -0.25,  0.25,  0.25];
    dN_dx = dN_dxi * (2/dx);
    dN_dy = dN_deta * (2/dy);
    grad_threshold = 0.05;
    epsilon = 0.1 * min(dx, dy);

    for ely = 1:nely
        for elx = 1:nelx
            base_i = ely + 1;
            base_j = elx + 1;

            node_coords = [base_i,   base_j;
                           base_i+1, base_j;
                           base_i+1, base_j+1;
                           base_i,   base_j+1];

            phi_nodes = zeros(1,4);
            for k = 1:4
                phi_nodes(k) = lsf(node_coords(k,1), node_coords(k,2));
            end

            dphi_dx = sum(dN_dx .* phi_nodes);
            dphi_dy = sum(dN_dy .* phi_nodes);
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
                phi_i = phi_nodes(k);
                Pi = sum(dN_dx .* phi_nodes) - dN_dx(k) * phi_i;
                Qi = sum(dN_dy .* phi_nodes) - dN_dy(k) * phi_i;
                A = dN_dx(k) * phi_i + Pi;
                B = dN_dy(k) * phi_i + Qi;
                denom = A^2 * ((B/A)^2 + 1);
                if abs(denom) < epsilon
                    if ~isempty(DIAG)
                        DIAG.den_small = DIAG.den_small + 1;
                    end
                    continue;
                end
                numerator = dN_dy(k) * (dN_dx(k) * phi_i + Pi) - dN_dx(k) * (dN_dy(k) * phi_i + Qi);
                if ~isfinite(numerator) || ~isfinite(denom) || abs(denom) < epsilon
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
            end

            if ~isempty(DIAG)
                DIAG.contrib_total = DIAG.contrib_total + 4;
            end
        end
    end
end

