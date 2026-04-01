function gradient_out = compute_gradient_chain_sensitivity(in)
%COMPUTE_GRADIENT_CHAIN_SENSITIVITY Build legacy/exact/chosen node gradients.

    legacy_element_sensitivity = compute_sensitivity_adjoint( ...
        in.nelx, in.nely, in.current_state.U, in.current_state.theta, ...
        in.E_L, in.E_T, in.nu_LT, in.G_LT, in.thickness, ...
        in.dx, in.dy, in.material_mask_core, in.normalize_sensitivity);
    legacy_node_sensitivity = aggregate_node_sensitivity( ...
        legacy_element_sensitivity, in.current_state.theta, in.lsf, ...
        in.nelx, in.nely, in.dx, in.dy, in.primary_update_mask);

    exact_theta_sensitivity = [];
    exact_node_sensitivity = [];
    exact_node_sensitivity_full = [];
    exact_node_sensitivity_opt = [];
    exact_pullback = [];
    exact_pullback_diag = struct();
    gradient_diag = struct();

    chain_mode = lower(string(in.gradient_opts.chain_mode));
    if any(strcmp(chain_mode, ["shadow", "exact"]))
        exact_theta_sensitivity = compute_sensitivity_adjoint( ...
            in.nelx, in.nely, in.theta_only_state.U, in.theta_only_state.theta, ...
            in.E_L, in.E_T, in.nu_LT, in.G_LT, in.thickness, ...
            in.dx, in.dy, in.material_mask_core, in.normalize_sensitivity);
        diff_opts = struct( ...
            'limiter_mode', in.gradient_opts.limiter_mode, ...
            'soft_limiter_beta', in.gradient_opts.soft_limiter_beta, ...
            'grad_floor', in.gradient_opts.theta_raw_grad_floor);
        [exact_pullback, exact_pullback_diag] = differentiate_theta_transition_exact( ...
            exact_theta_sensitivity, in.theta_only_state.theta_transition_cache, diff_opts);
        exact_node_sensitivity_full = aggregate_node_sensitivity( ...
            exact_pullback, [], in.lsf, in.nelx, in.nely, in.dx, in.dy, true(size(in.lsf)));
        exact_node_sensitivity_opt = aggregate_node_sensitivity( ...
            exact_pullback, [], in.lsf, in.nelx, in.nely, in.dx, in.dy, in.primary_update_mask);
        exact_node_sensitivity = exact_node_sensitivity_opt;
        gradient_diag = compute_gradient_chain_diagnostics( ...
            legacy_node_sensitivity, exact_node_sensitivity_opt, exact_node_sensitivity_full, ...
            in.primary_update_mask, ...
            in.material_mask_full, in.theta_only_state.theta_transition_cache, ...
            in.gradient_opts.shadow_topk);
        gradient_diag.pullback = exact_pullback_diag;
        gradient_diag.audit_mode = char(in.gradient_opts.audit_support_mode);
        gradient_diag.theta_raw_guard_ratio = get_struct_field_or_default( ...
            exact_pullback_diag, 'theta_raw_guard_ratio', NaN);
        gradient_diag.theta_raw_guard_nonzero_overlap = compute_theta_raw_guard_overlap( ...
            exact_theta_sensitivity, in.theta_only_state.theta_transition_cache, ...
            in.gradient_opts.theta_raw_grad_floor);
    end

    switch chain_mode
        case "exact"
            chosen_node_sensitivity = exact_node_sensitivity;
            chosen_source = 'exact';
        otherwise
            chosen_node_sensitivity = legacy_node_sensitivity;
            chosen_source = 'legacy';
    end

    gradient_out = struct();
    gradient_out.chosen_node_sensitivity = chosen_node_sensitivity;
    gradient_out.legacy_node_sensitivity = legacy_node_sensitivity;
    gradient_out.exact_node_sensitivity = exact_node_sensitivity;
    gradient_out.exact_node_sensitivity_full = exact_node_sensitivity_full;
    gradient_out.exact_node_sensitivity_opt = exact_node_sensitivity_opt;
    gradient_out.legacy_element_sensitivity = legacy_element_sensitivity;
    gradient_out.exact_theta_sensitivity = exact_theta_sensitivity;
    gradient_out.exact_pullback = exact_pullback;
    gradient_out.diagnostics = gradient_diag;
    gradient_out.chosen_source = chosen_source;
end

function overlap = compute_theta_raw_guard_overlap(theta_sensitivity, transition_cache, grad_floor)
    if isempty(theta_sensitivity)
        overlap = NaN;
        return;
    end

    grad_sq = transition_cache.angle_cache.dphi_dx .^ 2 + transition_cache.angle_cache.dphi_dy .^ 2;
    guard_mask = transition_cache.degenerate_grad_mask | (grad_sq < grad_floor * grad_floor);
    guard_mask = guard_mask & transition_cache.material_mask_core;
    support_mask = (abs(theta_sensitivity) > 1e-12) & transition_cache.material_mask_core;

    union_mask = guard_mask | support_mask;
    if any(union_mask(:))
        overlap = nnz(guard_mask & support_mask) / nnz(union_mask);
    else
        overlap = 1;
    end
end

function value = get_struct_field_or_default(data, field_name, default_value)
    if isstruct(data) && isfield(data, field_name)
        value = data.(field_name);
    else
        value = default_value;
    end
end
