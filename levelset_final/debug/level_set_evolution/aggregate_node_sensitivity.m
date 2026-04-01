function node_sensitivity = aggregate_node_sensitivity(sensitivity_input, ~, lsf, nelx, nely, dx, dy, primary_update_mask)
%AGGREGATE_NODE_SENSITIVITY Assemble node sensitivity for legacy/exact chains.

    primary_update_mask = normalize_primary_update_mask(primary_update_mask, lsf, nelx, nely);

    if isstruct(sensitivity_input) && isfield(sensitivity_input, 'mode')
        switch lower(string(sensitivity_input.mode))
            case "exact"
                node_sensitivity = assemble_exact_pullback(sensitivity_input, lsf, primary_update_mask, nelx, nely);
            otherwise
                error('aggregate_node_sensitivity:UnsupportedMode', ...
                    '不支持的pullback模式: %s', sensitivity_input.mode);
        end
    else
        node_sensitivity = aggregate_node_sensitivity_legacy( ...
            sensitivity_input, lsf, nelx, nely, dx, dy, primary_update_mask);
    end
end

function primary_update_mask = normalize_primary_update_mask(primary_update_mask, lsf, nelx, nely)
    if nargin < 1 || isempty(primary_update_mask)
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
end

function node_sensitivity = aggregate_node_sensitivity_legacy(element_sensitivity, lsf, nelx, nely, dx, dy, primary_update_mask)
    % 旧版局部角度链近似，保留用于legacy/shadow对照。

    global DIAG;
    node_sensitivity = zeros(size(lsf));

    dN_dxi = [-0.25,  0.25,  0.25, -0.25];
    dN_deta = [-0.25, -0.25,  0.25,  0.25];
    dN_dx = dN_dxi * (2/dx);
    dN_dy = dN_deta * (2/dy);
    grad_threshold = 0.05;
    eps_denom = 1e-12;

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
                Pi = dphi_dx - dN_dx(k) * phi_i;
                Qi = dphi_dy - dN_dy(k) * phi_i;
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
                node_sensitivity(node_coords(k,1), node_coords(k,2)) = ...
                    node_sensitivity(node_coords(k,1), node_coords(k,2)) + contrib;
                if ~isempty(DIAG) && contrib ~= 0
                    DIAG.contrib_nonzero = DIAG.contrib_nonzero + 1;
                end
            end

            if ~isempty(DIAG)
                DIAG.contrib_total = DIAG.contrib_total + 4;
            end
        end
    end

    node_sensitivity(~primary_update_mask) = 0;
end

function node_sensitivity = assemble_exact_pullback(pullback, lsf, primary_update_mask, nelx, nely)
    node_sensitivity = zeros(size(lsf));
    material_mask_core = logical(pullback.material_mask_core);

    for ely = 1:nely
        for elx = 1:nelx
            if ~material_mask_core(ely, elx)
                continue;
            end

            base_i = ely + 1;
            base_j = elx + 1;

            node_sensitivity(base_i, base_j) = node_sensitivity(base_i, base_j) + ...
                pullback.top_left_contrib(ely, elx);
            node_sensitivity(base_i + 1, base_j) = node_sensitivity(base_i + 1, base_j) + ...
                pullback.bottom_left_contrib(ely, elx);
            node_sensitivity(base_i + 1, base_j + 1) = node_sensitivity(base_i + 1, base_j + 1) + ...
                pullback.bottom_right_contrib(ely, elx);
            node_sensitivity(base_i, base_j + 1) = node_sensitivity(base_i, base_j + 1) + ...
                pullback.top_right_contrib(ely, elx);
        end
    end

    node_sensitivity(~isfinite(node_sensitivity)) = 0;
    node_sensitivity(~primary_update_mask) = 0;
end
