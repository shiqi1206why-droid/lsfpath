function [audit_mask, diagnostics] = build_theta_transition_audit_support(transition_cache, opts)
%BUILD_THETA_TRANSITION_AUDIT_SUPPORT Build differentiable node support for exact-chain audits.

    if nargin < 2 || isempty(opts)
        opts = struct();
    end

    grad_floor = get_option(opts, 'grad_floor', 0.05);
    branch_cut_tol = get_option(opts, 'branch_cut_tol', 1e-6);
    support_mode = lower(string(get_option(opts, 'support_mode', 'full')));
    primary_update_mask = get_option(opts, 'primary_update_mask', []);

    phi_base = transition_cache.phi_base;
    material_mask_full = get_option(opts, 'material_mask_full', []);
    if isempty(material_mask_full)
        material_mask_full = expand_material_mask_to_full_grid( ...
            transition_cache.material_mask_core, size(phi_base));
    end

    grad_sq = transition_cache.angle_cache.dphi_dx .^ 2 + transition_cache.angle_cache.dphi_dy .^ 2;
    grad_guard = transition_cache.degenerate_grad_mask | (grad_sq < grad_floor * grad_floor);

    theta_raw = transition_cache.theta_raw;
    branch_dist = min(theta_raw, pi - theta_raw);
    branch_guard = branch_dist <= branch_cut_tol;
    branch_guard(~transition_cache.material_mask_core) = false;

    cell_guard = (grad_guard | branch_guard) & transition_cache.material_mask_core;
    node_guard = guard_cells_to_nodes(cell_guard, size(phi_base));

    audit_mask = material_mask_full & ~node_guard;
    if ~isempty(primary_update_mask) && support_mode == "opt"
        audit_mask = audit_mask & logical(primary_update_mask);
    end
    audit_mask(~isfinite(phi_base)) = false;

    diagnostics = struct();
    diagnostics.grad_guard_ratio = safe_mask_fraction(grad_guard, transition_cache.material_mask_core);
    diagnostics.branch_guard_ratio = safe_mask_fraction(branch_guard, transition_cache.material_mask_core);
    diagnostics.cell_guard_ratio = safe_mask_fraction(cell_guard, transition_cache.material_mask_core);
    diagnostics.node_guard_ratio = safe_mask_fraction(node_guard, material_mask_full);
    diagnostics.audit_support_fraction = safe_mask_fraction(audit_mask, material_mask_full);
    diagnostics.support_mode = char(support_mode);
    diagnostics.grad_floor = grad_floor;
    diagnostics.branch_cut_tol = branch_cut_tol;
end

function node_guard = guard_cells_to_nodes(cell_guard, target_size)
    node_guard = false(target_size);
    [guard_i, guard_j] = find(cell_guard);
    for idx = 1:numel(guard_i)
        base_i = guard_i(idx) + 1;
        base_j = guard_j(idx) + 1;
        node_guard(base_i, base_j) = true;
        node_guard(base_i + 1, base_j) = true;
        node_guard(base_i + 1, base_j + 1) = true;
        node_guard(base_i, base_j + 1) = true;
    end
end

function material_mask_full = expand_material_mask_to_full_grid(material_mask_core, target_size)
    material_mask_full = false(target_size);
    material_mask_full(2:end-1, 2:end-1) = logical(material_mask_core);
    material_mask_full(1, :) = material_mask_full(2, :);
    material_mask_full(end, :) = material_mask_full(end-1, :);
    material_mask_full(:, 1) = material_mask_full(:, 2);
    material_mask_full(:, end) = material_mask_full(:, end-1);
end

function value = get_option(opts, name, default_value)
    value = default_value;
    if isfield(opts, name) && ~isempty(opts.(name))
        value = opts.(name);
    end
end

function ratio = safe_mask_fraction(mask, base_mask)
    mask = logical(mask);
    base_mask = logical(base_mask);
    denom = nnz(base_mask);
    if denom == 0
        ratio = 0;
    else
        ratio = nnz(mask & base_mask) / denom;
    end
end
