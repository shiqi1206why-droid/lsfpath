function iter_out = fiber_run_optimization_iterations(runtime_ctx, problem_ctx)
%FIBER_RUN_OPTIMIZATION_ITERATIONS Run the main optimization loop.

    params = runtime_ctx.params;

    nelx = runtime_ctx.nelx;
    nely = runtime_ctx.nely;
    dx = runtime_ctx.dx;
    dy = runtime_ctx.dy;

    E_L = runtime_ctx.E_L;
    E_T = runtime_ctx.E_T;
    nu_LT = runtime_ctx.nu_LT;
    G_LT = runtime_ctx.G_LT;
    thickness = runtime_ctx.thickness;

    max_iter = runtime_ctx.max_iter;
    tol = runtime_ctx.tol;
    delta_theta_max = runtime_ctx.delta_theta_max;
    enable_step_acceptance = runtime_ctx.enable_step_acceptance;
    acceptance_tol = runtime_ctx.acceptance_tol;
    enable_current_state_guard = runtime_ctx.enable_current_state_guard;
    current_state_tol = runtime_ctx.current_state_tol;
    current_guard_start_iter = runtime_ctx.current_guard_start_iter;
    enable_reinit_current_guard = runtime_ctx.enable_reinit_current_guard;
    reinit_current_tol = runtime_ctx.reinit_current_tol;
    reinit_guard_start_iter = runtime_ctx.reinit_guard_start_iter;
    backtrack_factor = runtime_ctx.backtrack_factor;
    max_backtrack = runtime_ctx.max_backtrack;
    min_backtrack_dt = runtime_ctx.min_backtrack_dt;
    enable_best_state_guard = runtime_ctx.enable_best_state_guard;
    best_state_rel_tol = runtime_ctx.best_state_rel_tol;
    best_state_patience = runtime_ctx.best_state_patience;
    theta_only_fuse_limit = runtime_ctx.theta_only_fuse_limit;
    velocity_opts = runtime_ctx.velocity_opts;
    gradient_opts = runtime_ctx.gradient_opts;
    manufacturing_opts = runtime_ctx.manufacturing_opts;

    F_mag = runtime_ctx.F_mag;
    hj_update_opts = runtime_ctx.hj_update_opts;
    reinit_opts = runtime_ctx.reinit_opts;
    path_quality_opts = runtime_ctx.path_quality_opts;

    material_mask_core = problem_ctx.material_mask_core;
    material_mask_full = problem_ctx.material_mask_full;
    h_grid = problem_ctx.h_grid;
    h = problem_ctx.h;
    boundary_guard_band = problem_ctx.boundary_guard_band;
    lsf_target_global = problem_ctx.lsf_target_global;

    compliance_history = problem_ctx.compliance_history;
    FCS_history = problem_ctx.FCS_history;
    theta_only_compliance_history = problem_ctx.theta_only_compliance_history;
    hj_raw_compliance_history = problem_ctx.hj_raw_compliance_history;
    hj_trial_compliance_history = problem_ctx.hj_trial_compliance_history;
    hj_local_reinit_compliance_history = problem_ctx.hj_local_reinit_compliance_history;
    reinit_trial_compliance_history = problem_ctx.reinit_trial_compliance_history;
    accepted_source_history = problem_ctx.accepted_source_history;
    accepted_source_detail_history = problem_ctx.accepted_source_detail_history;
    accepted_steps = problem_ctx.accepted_steps;
    rejected_steps = problem_ctx.rejected_steps;
    theta_only_accept_count = problem_ctx.theta_only_accept_count;
    theta_only_reject_count = problem_ctx.theta_only_reject_count;
    theta_only_consecutive_count = problem_ctx.theta_only_consecutive_count;
    reject_due_next_guard = problem_ctx.reject_due_next_guard;
    reject_due_current_guard = problem_ctx.reject_due_current_guard;
    reinit_trigger_count = problem_ctx.reinit_trigger_count;
    reinit_skip_due_reject = problem_ctx.reinit_skip_due_reject;
    reinit_skip_due_objective = problem_ctx.reinit_skip_due_objective;
    reinit_skip_due_next_guard = problem_ctx.reinit_skip_due_next_guard;
    reinit_skip_due_current_guard = problem_ctx.reinit_skip_due_current_guard;
    iter_since_last_reinit = problem_ctx.iter_since_last_reinit;
    best_state = problem_ctx.best_state;
    no_improve_counter = problem_ctx.no_improve_counter;
    early_stop_triggered = problem_ctx.early_stop_triggered;
    early_stop_reason = problem_ctx.early_stop_reason;
    history_count = problem_ctx.history_count;
    current_state_recorded = problem_ctx.current_state_recorded;
    prev_theta_for_diag = problem_ctx.prev_theta_for_diag;
    candidate_select_tol = problem_ctx.candidate_select_tol;
    loop_iter_count = problem_ctx.loop_iter_count;
    quality_eval_stride = problem_ctx.quality_eval_stride;
    raw_turn_history = problem_ctx.raw_turn_history;
    raw_kappa_history = problem_ctx.raw_kappa_history;
    raw_spacing_error_history = problem_ctx.raw_spacing_error_history;
    raw_grad_dev_history = problem_ctx.raw_grad_dev_history;
    raw_near_zero_outlier_history = problem_ctx.raw_near_zero_outlier_history;
    hj_fallback_history = problem_ctx.hj_fallback_history;
    hj_second_order_history = problem_ctx.hj_second_order_history;
    hj_frozen_incomplete_history = problem_ctx.hj_frozen_incomplete_history;
    hj_first_order_complete_history = problem_ctx.hj_first_order_complete_history;
    reinit_method_history = problem_ctx.reinit_method_history;
    reinit_fallback_history = problem_ctx.reinit_fallback_history;
    reinit_reason_history = problem_ctx.reinit_reason_history;
    hj_update_diagnostics_history = problem_ctx.hj_update_diagnostics_history;
    reinit_diagnostics_history = problem_ctx.reinit_diagnostics_history;
    boundary_guard_ratio_history = problem_ctx.boundary_guard_ratio_history;
    frozen_boundary_point_history = problem_ctx.frozen_boundary_point_history;
    accepted_hj_reinit_count = problem_ctx.accepted_hj_reinit_count;
    local_reinit_shell_size_history = problem_ctx.local_reinit_shell_size_history;
    post_reinit_grad_dev_mean_history = problem_ctx.post_reinit_grad_dev_mean_history;
    post_reinit_grad_outlier_ratio_history = problem_ctx.post_reinit_grad_outlier_ratio_history;
    refresh_shell_size_history = problem_ctx.refresh_shell_size_history;
    refresh_count = problem_ctx.refresh_count;
    last_velocity_field = problem_ctx.last_velocity_field;
    last_propagation_mask = problem_ctx.last_propagation_mask;
    last_hj_info = problem_ctx.last_hj_info;
    last_reinit_info = problem_ctx.last_reinit_info;
    gradient_chain_cosine_history = problem_ctx.gradient_chain_cosine_history;
    gradient_chain_norm_ratio_history = problem_ctx.gradient_chain_norm_ratio_history;
    gradient_chain_topk_sign_history = problem_ctx.gradient_chain_topk_sign_history;
    gradient_chain_saturation_history = problem_ctx.gradient_chain_saturation_history;
    gradient_chain_degenerate_history = problem_ctx.gradient_chain_degenerate_history;
    gradient_chain_active_band_coverage_history = problem_ctx.gradient_chain_active_band_coverage_history;
    gradient_chain_zero_limiter_history = problem_ctx.gradient_chain_zero_limiter_history;
    gradient_chain_exact_nonzero_history = problem_ctx.gradient_chain_exact_nonzero_history;
    gradient_chain_support_overlap_history = problem_ctx.gradient_chain_support_overlap_history;
    gradient_chain_theta_raw_guard_history = problem_ctx.gradient_chain_theta_raw_guard_history;
    gradient_chain_theta_raw_guard_overlap_history = problem_ctx.gradient_chain_theta_raw_guard_overlap_history;
    gradient_chain_full_vs_opt_overlap_history = problem_ctx.gradient_chain_full_vs_opt_overlap_history;
    gradient_chain_selected_source_history = problem_ctx.gradient_chain_selected_source_history;
    gradient_chain_audit_mode = problem_ctx.gradient_chain_audit_mode;
    last_gradient_chain_diag = problem_ctx.last_gradient_chain_diag;
    manufacturing_grad_norm_history = problem_ctx.manufacturing_grad_norm_history;
    manufacturing_curvature_norm_history = problem_ctx.manufacturing_curvature_norm_history;
    manufacturing_gap_overlap_norm_history = problem_ctx.manufacturing_gap_overlap_norm_history;
    theta_only_vs_current_history = problem_ctx.theta_only_vs_current_history;
    hj_raw_vs_theta_only_history = problem_ctx.hj_raw_vs_theta_only_history;
    reinit_vs_theta_only_history = problem_ctx.reinit_vs_theta_only_history;
    last_manufacturing_diag = problem_ctx.last_manufacturing_diag;

    current_state = problem_ctx.current_state;
for iter = 1:max_iter
    lsf = current_state.lsf;
    theta_e = current_state.theta;
    theta_target = current_state.theta_target;
    U = current_state.U;
    K = current_state.K;
    F = current_state.F;
    compliance = current_state.compliance;
    strain_energy = current_state.strain_energy;
    FCS = current_state.FCS;

    history_state = struct();
    history_state.current_state_recorded = current_state_recorded;
    history_state.history_count = history_count;
    history_state.quality_eval_stride = quality_eval_stride;
    history_state.compliance_history = compliance_history;
    history_state.FCS_history = FCS_history;
    history_state.raw_turn_history = raw_turn_history;
    history_state.raw_kappa_history = raw_kappa_history;
    history_state.raw_spacing_error_history = raw_spacing_error_history;
    history_state.raw_grad_dev_history = raw_grad_dev_history;
    history_state.raw_near_zero_outlier_history = raw_near_zero_outlier_history;
    history_state = fiber_record_history_entry(history_state, current_state, dx, dy, material_mask_core, path_quality_opts);
    current_state_recorded = history_state.current_state_recorded;
    history_count = history_state.history_count;
    compliance_history = history_state.compliance_history;
    FCS_history = history_state.FCS_history;
    raw_turn_history = history_state.raw_turn_history;
    raw_kappa_history = history_state.raw_kappa_history;
    raw_spacing_error_history = history_state.raw_spacing_error_history;
    raw_grad_dev_history = history_state.raw_grad_dev_history;
    raw_near_zero_outlier_history = history_state.raw_near_zero_outlier_history;

    % === 优化1.3：预计算常用掩码（每次迭代开始） ===
    abs_lsf = abs(lsf);
    bands = struct();
    bands.narrow_05h = abs_lsf <= 0.5 * h_grid;
    bands.narrow_10h = abs_lsf <= 1.0 * h_grid;
    bands.narrow_15h = abs_lsf <= 1.5 * h_grid;

    if iter == 1
        center_i = round((nely+2)/2);
        center_j = round((nelx+2)/2);
        dphi_dx_center = (lsf(center_i, center_j+1) - lsf(center_i, center_j-1)) / (2*dx);
        dphi_dy_center = (lsf(center_i+1, center_j) - lsf(center_i-1, center_j)) / (2*dy);
        grad_mag = hypot(dphi_dx_center, dphi_dy_center);
        fprintf('  中心位置梯度：dx=%e，dy=%e，|grad|=%e\n', dphi_dx_center, dphi_dy_center, grad_mag);
    end

    if iter == 1 || mod(iter, 10) == 0
        compliance_Ck = full(U' * K * U);
        energy_error = abs(compliance - compliance_Ck);
        relative_error = energy_error / max(1e-12, abs(compliance));
        if relative_error > 1e-6
            warning('迭代%d: U''*F 与 U''*K*U 不一致 (误差=%.3e, 相对误差=%.2e)', ...
                iter, energy_error, relative_error);
        end
        fprintf('  [能量自检] U''*F = %.6e, U''*K*U = %.6e, 差异 = %.2e\n', ...
            compliance, compliance_Ck, energy_error);
    end

    if iter == 1
        fprintf('\n=== 优化诊断 ===\n');
        fprintf('初始柔度: %e\n', compliance);
        fprintf('施加载荷大小: %e\n', F_mag);
    end

    if ~isempty(prev_theta_for_diag)
        theta_delta = atan2(sin(theta_e(:) - prev_theta_for_diag(:)), cos(theta_e(:) - prev_theta_for_diag(:)));
        fprintf('  角度变化：最大=%.2f 度，平均=%.2f 度\n', ...
            max(abs(theta_delta)) * 180/pi, mean(abs(theta_delta)) * 180/pi);
    end

    if ~isfinite(best_state.compliance) || compliance < best_state.compliance * (1 - best_state_rel_tol)
        best_state.compliance = compliance;
        best_state.iter = history_count;
        best_state.lsf = lsf;
        best_state.theta = theta_e;
        best_state.strain_energy = strain_energy;
        best_state.FCS = FCS;
        best_state.theta_target = theta_target;
        no_improve_counter = 0;
    else
        no_improve_counter = no_improve_counter + 1;
    end

    if enable_best_state_guard && iter >= 20 && no_improve_counter >= best_state_patience
        early_stop_triggered = true;
        early_stop_reason = sprintf('连续%d步无改进(best_iter=%d)', no_improve_counter, best_state.iter);
        log_message('INFO', params, '触发最优守护早停: %s', early_stop_reason);
        save_checkpoint(iter, lsf, theta_e, compliance_history(1:history_count), FCS_history(1:history_count), params);
        break;
    end

    theta_only_state = evaluate_candidate_state(lsf, theta_e, delta_theta_max, dx, dy, ...
        nelx, nely, material_mask_core, E_L, E_T, nu_LT, G_LT, thickness, F_mag, ...
        params.smooth.eta, params.smooth.iterations);
    theta_only_compliance = theta_only_state.compliance;
    theta_only_compliance_history(iter) = theta_only_compliance;
    theta_only_vs_current_history(iter) = compute_candidate_delta(theta_only_compliance, compliance);
    theta_only_delta = atan2(sin(theta_only_state.theta(:) - theta_e(:)), cos(theta_only_state.theta(:) - theta_e(:)));
    theta_only_state_changed = max(abs(theta_only_delta)) > 1e-12;
    theta_only_tol = compute_theta_only_tol(iter, max_iter, current_state_tol);
    theta_only_fuse_active = theta_only_consecutive_count >= theta_only_fuse_limit && no_improve_counter > 5;
    if theta_only_fuse_active
        if theta_only_consecutive_count == theta_only_fuse_limit
            log_message('INFO', params, ...
                'theta_only 熔断触发：连续%d步 theta_only 接管，后续仅允许不恶化。', ...
                theta_only_consecutive_count);
        end
        theta_only_tol = 0;
    end
    theta_only_ok_current = isfinite(theta_only_compliance) && ...
        theta_only_compliance <= compliance * (1 + theta_only_tol);
    theta_only_accepted = theta_only_state_changed && theta_only_ok_current;
    if theta_only_state_changed
        if theta_only_accepted
            theta_only_accept_count = theta_only_accept_count + 1;
        else
            theta_only_reject_count = theta_only_reject_count + 1;
        end
    end

    % 3.5 灵敏度与速度场
    primary_update_mask = bands.narrow_15h & material_mask_full & ~boundary_guard_band;
    stencil_mask = dilate_binary_mask(primary_update_mask, params.levelset.stencil_buffer_cells) & material_mask_full;
    boundary_frozen_mask = bands.narrow_15h & material_mask_full & boundary_guard_band;
    boundary_guard_ratio_history(iter) = safe_fraction(nnz(boundary_frozen_mask), nnz(bands.narrow_15h & material_mask_full));
    frozen_boundary_point_history(iter) = nnz(boundary_frozen_mask);

    gradient_in = struct();
    gradient_in.current_state = current_state;
    gradient_in.theta_only_state = theta_only_state;
    gradient_in.lsf = lsf;
    gradient_in.nelx = nelx;
    gradient_in.nely = nely;
    gradient_in.dx = dx;
    gradient_in.dy = dy;
    gradient_in.material_mask_core = material_mask_core;
    gradient_in.material_mask_full = material_mask_full;
    gradient_in.primary_update_mask = primary_update_mask;
    gradient_in.gradient_opts = gradient_opts;
    gradient_in.normalize_sensitivity = params.opt.normalize_sensitivity;
    gradient_in.E_L = E_L;
    gradient_in.E_T = E_T;
    gradient_in.nu_LT = nu_LT;
    gradient_in.G_LT = G_LT;
    gradient_in.thickness = thickness;
    gradient_out = compute_gradient_chain_sensitivity(gradient_in);
    node_sensitivity = gradient_out.chosen_node_sensitivity;
    legacy_node_sensitivity = gradient_out.legacy_node_sensitivity;
    exact_node_sensitivity = gradient_out.exact_node_sensitivity;
    exact_node_sensitivity_full = gradient_out.exact_node_sensitivity_full;
    gradient_chain_diag = gradient_out.diagnostics;
    last_gradient_chain_diag = gradient_chain_diag;

    if iter == 1 || mod(iter, 10) == 0
        fprintf('迭代 %d - 节点灵敏度统计：最大=%e，最小=%e，平均=%e，source=%s\n', iter, ...
            max(node_sensitivity(:)), min(node_sensitivity(:)), mean(abs(node_sensitivity(:))), ...
            gradient_out.chosen_source);
    end

    gradient_chain_selected_source_history(iter) = string(gradient_out.chosen_source);
    if ~isempty(fieldnames(gradient_chain_diag))
        gradient_chain_cosine_history(iter) = get_struct_field_or_default(gradient_chain_diag, 'cosine_similarity', NaN);
        gradient_chain_norm_ratio_history(iter) = get_struct_field_or_default(gradient_chain_diag, 'norm_ratio', NaN);
        gradient_chain_topk_sign_history(iter) = get_struct_field_or_default(gradient_chain_diag, 'topk_sign_agreement', NaN);
        gradient_chain_saturation_history(iter) = get_struct_field_or_default(gradient_chain_diag, 'saturation_ratio', NaN);
        gradient_chain_degenerate_history(iter) = get_struct_field_or_default(gradient_chain_diag, 'degenerate_ratio', NaN);
        gradient_chain_active_band_coverage_history(iter) = get_struct_field_or_default(gradient_chain_diag, 'active_band_coverage', NaN);
        gradient_chain_zero_limiter_history(iter) = get_struct_field_or_default(gradient_chain_diag, 'zero_gradient_due_to_limiter_ratio', NaN);
        gradient_chain_exact_nonzero_history(iter) = get_struct_field_or_default(gradient_chain_diag, 'exact_gradient_active_nonzero_ratio', NaN);
        gradient_chain_support_overlap_history(iter) = get_struct_field_or_default(gradient_chain_diag, 'exact_vs_legacy_band_support_overlap', NaN);
        gradient_chain_theta_raw_guard_history(iter) = get_struct_field_or_default(gradient_chain_diag, 'theta_raw_guard_ratio', NaN);
        gradient_chain_theta_raw_guard_overlap_history(iter) = get_struct_field_or_default(gradient_chain_diag, 'theta_raw_guard_nonzero_overlap', NaN);
        gradient_chain_full_vs_opt_overlap_history(iter) = get_struct_field_or_default(gradient_chain_diag, 'full_vs_opt_support_overlap', NaN);
    end

    if iter == 1 || mod(iter, 10) == 0
        sens_diag = abs(node_sensitivity(:));
        sens_diag = sens_diag(isfinite(sens_diag));
        if ~isempty(sens_diag)
            p95 = prctile(sens_diag, 95);
            p99 = prctile(sens_diag, 99);
            maxv = max(sens_diag);
            ratio_p95 = mean(sens_diag > p95) * 100;
            ratio_p99 = mean(sens_diag > p99) * 100;
            fprintf('  [灵敏度分位] p95=%.2e, p99=%.2e, max=%.2e, >p95=%.2f%%, >p99=%.2f%%\n', ...
                p95, p99, maxv, ratio_p95, ratio_p99);
        end
    end

    [manufacturing_gradient, manufacturing_diag] = compute_manufacturing_penalty_gradient( ...
        lsf, dx, dy, material_mask_full, lsf_target_global, primary_update_mask, manufacturing_opts);
    node_sensitivity = node_sensitivity + manufacturing_gradient;
    last_manufacturing_diag = manufacturing_diag;
    manufacturing_grad_norm_history(iter) = manufacturing_diag.grad_norm_grad_norm;
    manufacturing_curvature_norm_history(iter) = manufacturing_diag.curvature_grad_norm;
    manufacturing_gap_overlap_norm_history(iter) = manufacturing_diag.gap_overlap_grad_norm;

    sens_scale_value = NaN;
    sens_abs_band = abs(node_sensitivity(primary_update_mask));
    if ~isempty(sens_abs_band) && velocity_opts.scale_quantile > 0
        sens_scale_value = prctile(sens_abs_band, velocity_opts.scale_quantile);
        if sens_scale_value > 0
            node_sensitivity = node_sensitivity / sens_scale_value;
        end
    end
    if isfinite(velocity_opts.clip_abs)
        node_sensitivity = max(min(node_sensitivity, velocity_opts.clip_abs), -velocity_opts.clip_abs);
    end

    gamma_curv = 0.5 * h_grid;
    if params.opt.enable_curvature
        curvature_gamma = gamma_curv;
    else
        curvature_gamma = 0;
    end
    if strcmp(gradient_out.chosen_source, 'exact')
        [velocity, velocity_stats] = build_velocity_exact(node_sensitivity, lsf, dx, dy, 1.5*h_grid, ...
            velocity_opts.enable_bias_removal, curvature_gamma, iter, primary_update_mask, velocity_opts);
    else
        [velocity, velocity_stats] = build_velocity_field(node_sensitivity, lsf, dx, dy, 1.5*h_grid, ...
            velocity_opts.enable_bias_removal, curvature_gamma, iter, primary_update_mask, velocity_opts);
    end
    last_velocity_field = velocity;
    last_propagation_mask = primary_update_mask;

    log_diag_in = struct();
    log_diag_in.iter = iter;
    log_diag_in.lsf = lsf;
    log_diag_in.dx = dx;
    log_diag_in.dy = dy;
    log_diag_in.h = h;
    log_diag_in.lsf_target_global = lsf_target_global;
    log_diag_in.node_sensitivity = node_sensitivity;
    log_diag_in.sens_scale_value = sens_scale_value;
    log_diag_in.velocity_stats = velocity_stats;
    log_diag_in.bands = bands;
    log_diag_in.material_mask_full = material_mask_full;
    log_diag_in.material_mask_core = material_mask_core;
    log_diag_in.path_quality_opts = path_quality_opts;
    fiber_log_velocity_path_diagnostics(log_diag_in);

    narrow_band_mask = primary_update_mask;
    band_velocity = velocity;
    band_velocity(~narrow_band_mask) = 0;
    dt_cfl = compute_adaptive_timestep(band_velocity, dx, dy);
    dt_adaptive = dt_cfl;
    if velocity_stats.max_band > 1e-12
        dt_angle = delta_theta_max / velocity_stats.max_band;
        dt_adaptive = min(dt_adaptive, dt_angle);
    else
        dt_angle = inf;
    end

    log_summary_in = struct();
    log_summary_in.iter = iter;
    log_summary_in.node_sensitivity = node_sensitivity;
    log_summary_in.velocity = velocity;
    log_summary_in.velocity_stats = velocity_stats;
    log_summary_in.dt_adaptive = dt_adaptive;
    log_summary_in.dt_cfl = dt_cfl;
    log_summary_in.dt_angle = dt_angle;
    log_summary_in.gradient_chain_diag = gradient_chain_diag;
    log_summary_in.gradient_chain_source = gradient_out.chosen_source;
    log_summary_in.legacy_node_sensitivity = legacy_node_sensitivity;
    log_summary_in.exact_node_sensitivity = exact_node_sensitivity;
    log_summary_in.exact_node_sensitivity_full = exact_node_sensitivity_full;
    log_summary_in.manufacturing_diag = manufacturing_diag;
    fiber_log_sensitivity_velocity_summary(log_summary_in);

    % 3.6 HJ候选
    lsf_before = lsf;
    dt_trial = dt_adaptive;
    hj_raw_compliance = NaN;
    hj_trial_compliance = NaN;
    hj_local_reinit_compliance = NaN;
    reinit_trial_compliance = NaN;
    hj_raw_state = [];
    hj_trial_state = [];
    hj_local_reinit_state = [];
    reinit_trial_state = [];
    step_accepted = true;
    backtrack_used = 0;
    fail_next_guard = false;
    fail_current_guard = false;
    current_guard_active = enable_current_state_guard && (iter >= current_guard_start_iter);
    ref_next_compliance = compliance;
    if theta_only_accepted
        ref_next_compliance = theta_only_compliance;
    end

    step_in = struct();
    step_in.enable_step_acceptance = enable_step_acceptance;
    step_in.max_backtrack = max_backtrack;
    step_in.dt_adaptive = dt_adaptive;
    step_in.backtrack_factor = backtrack_factor;
    step_in.min_backtrack_dt = min_backtrack_dt;
    step_in.lsf_before = lsf_before;
    step_in.velocity = velocity;
    step_in.dx = dx;
    step_in.dy = dy;
    step_in.primary_update_mask = primary_update_mask;
    step_in.stencil_mask = stencil_mask;
    step_in.hj_update_opts = hj_update_opts;
    step_in.theta_e = theta_e;
    step_in.delta_theta_max = delta_theta_max;
    step_in.nelx = nelx;
    step_in.nely = nely;
    step_in.material_mask_core = material_mask_core;
    step_in.E_L = E_L;
    step_in.E_T = E_T;
    step_in.nu_LT = nu_LT;
    step_in.G_LT = G_LT;
    step_in.thickness = thickness;
    step_in.F_mag = F_mag;
    step_in.smooth_eta = params.smooth.eta;
    step_in.smooth_iterations = params.smooth.iterations;
    step_in.ref_next_compliance = ref_next_compliance;
    step_in.acceptance_tol = acceptance_tol;
    step_in.current_guard_active = current_guard_active;
    step_in.compliance = compliance;
    step_in.current_state_tol = current_state_tol;
    step_in.iter = iter;
    step_in.accepted_steps = accepted_steps;
    step_in.rejected_steps = rejected_steps;
    step_in.reject_due_next_guard = reject_due_next_guard;
    step_in.reject_due_current_guard = reject_due_current_guard;
    step_out = fiber_run_hj_backtracking_step(step_in);
    hj_raw_compliance = step_out.hj_trial_compliance;
    hj_raw_state = step_out.hj_trial_state;
    hj_trial_compliance = hj_raw_compliance;
    hj_trial_state = hj_raw_state;
    step_accepted = step_out.step_accepted;
    backtrack_used = step_out.backtrack_used;
    fail_next_guard = step_out.fail_next_guard;
    fail_current_guard = step_out.fail_current_guard;
    accepted_steps = step_out.accepted_steps;
    rejected_steps = step_out.rejected_steps;
    reject_due_next_guard = step_out.reject_due_next_guard;
    reject_due_current_guard = step_out.reject_due_current_guard;
    if ~isempty(fieldnames(step_out.hj_info_trial))
        hj_info_trial = step_out.hj_info_trial;
        hj_fallback_history(iter) = hj_info_trial.fallback_count;
        hj_second_order_history(iter) = hj_info_trial.second_order_count;
        hj_frozen_incomplete_history(iter) = get_struct_field_or_default(hj_info_trial, 'frozen_incomplete_godunov_count', NaN);
        hj_first_order_complete_history(iter) = get_struct_field_or_default(hj_info_trial, 'used_first_order_complete_count', NaN);
        hj_update_diagnostics_history{iter} = hj_info_trial;
        last_hj_info = hj_info_trial;
    end
    hj_raw_compliance_history(iter) = hj_raw_compliance;
    hj_trial_compliance_history(iter) = hj_trial_compliance;
    hj_local_reinit_compliance_history(iter) = hj_local_reinit_compliance;
    hj_raw_vs_theta_only_history(iter) = compute_candidate_delta(hj_raw_compliance, theta_only_compliance);

    % 接受HJ候选后，先执行局部强制重初始化，再按原有条件决定是否做更大范围重初始化。
    reinit_candidate_selected = false;
    iter_since_last_reinit_candidate = iter_since_last_reinit;
    if step_accepted && ~isempty(hj_trial_state)
        iter_since_last_reinit_candidate = iter_since_last_reinit + 1;
        if strcmpi(params.levelset.reinit_domain, 'masked')
            reinit_mask = material_mask_core;
        else
            reinit_mask = [];
        end

        zero_mask_dynamic = compute_zero_mask_from_lsf(hj_raw_state.lsf, h_grid);
        zero_mask_dynamic = zero_mask_dynamic & material_mask_full;
        zero_band_local = abs(hj_raw_state.lsf) <= params.levelset.local_reinit_zero_band_factor * h_grid;
        local_reinit_shell = dilate_binary_mask(zero_band_local, params.levelset.local_reinit_buffer_cells) & material_mask_full;
        local_reinit_shell_size_history(iter) = nnz(local_reinit_shell);

        if any(local_reinit_shell(:))
            accepted_hj_reinit_count = accepted_hj_reinit_count + 1;
            local_reinit_opts = reinit_opts;
            local_reinit_opts.local_shell_mask = local_reinit_shell;
            local_reinit_opts.preserve_outside_shell = true;
            [lsf_local_reinit, local_reinit_diag] = fmm_reinitialize(hj_raw_state.lsf, dx, dy, ...
                zero_mask_dynamic, reinit_mask, local_reinit_opts);
            local_reinit_state = evaluate_candidate_state(lsf_local_reinit, theta_e, delta_theta_max, dx, dy, ...
                nelx, nely, material_mask_core, E_L, E_T, nu_LT, G_LT, thickness, F_mag, ...
                params.smooth.eta, params.smooth.iterations);
            local_reinit_compliance = local_reinit_state.compliance;
            post_local_metrics = compute_raw_path_quality_metrics(lsf_local_reinit, dx, dy, material_mask_core, path_quality_opts);
            post_reinit_grad_dev_mean_history(iter) = post_local_metrics.grad_dev_mean;
            post_reinit_grad_outlier_ratio_history(iter) = post_local_metrics.near_zero_grad_outlier_ratio;

            reinit_ok_next = isfinite(local_reinit_compliance) && isfinite(hj_raw_compliance) && ...
                local_reinit_compliance <= hj_raw_compliance * (1 + acceptance_tol);
            reinit_guard_active = enable_reinit_current_guard && (iter >= reinit_guard_start_iter);
            if reinit_guard_active
                reinit_ok_current = isfinite(local_reinit_compliance) && ...
                    local_reinit_compliance <= compliance * (1 + reinit_current_tol);
            else
                reinit_ok_current = true;
            end

            reinit_method_history(iter) = string(local_reinit_diag.method_used);
            reinit_fallback_history(iter) = logical(local_reinit_diag.fallback_used);
            reinit_reason_history(iter) = "accepted_hj_local_shell";
            reinit_diagnostics_history{iter} = local_reinit_diag;
            last_reinit_info = local_reinit_diag;
            if ~(reinit_ok_next && reinit_ok_current)
                log_message('INFO', params, ...
                    ['局部重初始化已并入HJ状态: C_local=%.4e, C_hj=%.4e, C_current=%.4e, ' ...
                    'next_ok=%d, current_ok=%d, guard_active=%d'], ...
                    local_reinit_compliance, hj_raw_compliance, compliance, ...
                    reinit_ok_next, reinit_ok_current, reinit_guard_active);
            end
            hj_local_reinit_state = local_reinit_state;
            hj_local_reinit_compliance = local_reinit_compliance;
            hj_trial_state = hj_local_reinit_state;
            hj_trial_compliance = hj_local_reinit_compliance;
            hj_trial_compliance_history(iter) = hj_trial_compliance;
            hj_local_reinit_compliance_history(iter) = hj_local_reinit_compliance;
            iter_since_last_reinit_candidate = 0;
        end

        reinit_reference_state = hj_trial_state;
        reinit_reference_compliance = hj_trial_compliance;
        if reinit_candidate_selected && ~isempty(reinit_trial_state)
            reinit_reference_state = reinit_trial_state;
            reinit_reference_compliance = reinit_trial_compliance;
        end

        [do_reinit, reinit_reason] = should_reinitialize(reinit_reference_state.lsf, lsf_before, ...
            iter_since_last_reinit_candidate, iter, params);
        if do_reinit
            log_message('INFO', params, '触发重初始化: %s', reinit_reason);
            zero_mask_dynamic = compute_zero_mask_from_lsf(reinit_reference_state.lsf, h_grid);
            zero_mask_dynamic = zero_mask_dynamic & material_mask_full;
            global_reinit_opts = reinit_opts;
            global_reinit_opts.local_shell_mask = [];
            global_reinit_opts.preserve_outside_shell = false;
            [lsf_reinit_candidate, reinit_diag] = fmm_reinitialize(reinit_reference_state.lsf, dx, dy, ...
                zero_mask_dynamic, reinit_mask, global_reinit_opts);
            reinit_method_history(iter) = string(reinit_diag.method_used);
            reinit_fallback_history(iter) = logical(reinit_diag.fallback_used);
            reinit_reason_history(iter) = string(reinit_reason);
            reinit_diagnostics_history{iter} = reinit_diag;
            last_reinit_info = reinit_diag;
            global_reinit_state = evaluate_candidate_state(lsf_reinit_candidate, theta_e, delta_theta_max, dx, dy, ...
                nelx, nely, material_mask_core, E_L, E_T, nu_LT, G_LT, thickness, F_mag, ...
                params.smooth.eta, params.smooth.iterations);
            global_reinit_compliance = global_reinit_state.compliance;
            post_global_metrics = compute_raw_path_quality_metrics(lsf_reinit_candidate, dx, dy, material_mask_core, path_quality_opts);
            post_reinit_grad_dev_mean_history(iter) = post_global_metrics.grad_dev_mean;
            post_reinit_grad_outlier_ratio_history(iter) = post_global_metrics.near_zero_grad_outlier_ratio;

            reinit_ok_next = isfinite(global_reinit_compliance) && isfinite(reinit_reference_compliance) && ...
                global_reinit_compliance <= reinit_reference_compliance * (1 + acceptance_tol);
            reinit_guard_active = enable_reinit_current_guard && (iter >= reinit_guard_start_iter);
            if reinit_guard_active
                reinit_ok_current = isfinite(global_reinit_compliance) && ...
                    global_reinit_compliance <= compliance * (1 + reinit_current_tol);
            else
                reinit_ok_current = true;
            end

            if ~(reinit_ok_next && reinit_ok_current)
                reinit_skip_due_objective = reinit_skip_due_objective + 1;
                if ~reinit_ok_next
                    reinit_skip_due_next_guard = reinit_skip_due_next_guard + 1;
                end
                if ~reinit_ok_current
                    reinit_skip_due_current_guard = reinit_skip_due_current_guard + 1;
                end
                log_message('INFO', params, ...
                    ['重初始化回退: C_reinit=%.4e, C_ref=%.4e, C_current=%.4e, ' ...
                    'next_ok=%d, current_ok=%d, reinit_guard_active=%d'], ...
                    global_reinit_compliance, reinit_reference_compliance, compliance, ...
                    reinit_ok_next, reinit_ok_current, reinit_guard_active);
            else
                reinit_candidate_selected = true;
                reinit_trial_state = global_reinit_state;
                reinit_trial_compliance = global_reinit_compliance;
            end
        end
    else
        reinit_skip_due_reject = reinit_skip_due_reject + 1;
        if iter == 1 || mod(iter, 10) == 0
            fprintf('  [重初始化] 本步拒绝更新，跳过固定频率重初始化\n');
        end
    end
    reinit_trial_compliance_history(iter) = reinit_trial_compliance;
    reinit_vs_theta_only_history(iter) = compute_candidate_delta(reinit_trial_compliance, theta_only_compliance);

    sel_in = struct();
    sel_in.current_state = current_state;
    sel_in.theta_only_accepted = theta_only_accepted;
    sel_in.theta_only_compliance = theta_only_compliance;
    sel_in.theta_only_state = theta_only_state;
    sel_in.step_accepted = step_accepted;
    sel_in.hj_raw_state = hj_raw_state;
    sel_in.hj_local_reinit_state = hj_local_reinit_state;
    sel_in.hj_trial_state = hj_trial_state;
    sel_in.reinit_candidate_selected = reinit_candidate_selected;
    sel_in.reinit_trial_state = reinit_trial_state;
    sel_in.candidate_select_tol = candidate_select_tol;
    sel_in.acceptance_tol = acceptance_tol;
    [next_state, accepted_source, accepted_source_detail] = fiber_select_candidate_state(sel_in);

    if step_accepted && any(strcmp(accepted_source, {'hj', 'reinit'})) && ...
            mod(accepted_steps, params.levelset.refresh_interval) == 0
        refresh_shell = dilate_binary_mask(primary_update_mask, params.levelset.refresh_buffer_cells) & ...
            material_mask_full & ~primary_update_mask;
        refresh_shell_size_history(iter) = nnz(refresh_shell);
        if any(refresh_shell(:))
            if strcmpi(params.levelset.reinit_domain, 'masked')
                refresh_mask = material_mask_core;
            else
                refresh_mask = [];
            end
            refresh_zero_mask = compute_zero_mask_from_lsf(next_state.lsf, h_grid);
            refresh_zero_mask = refresh_zero_mask & material_mask_full;
            refresh_opts = reinit_opts;
            refresh_opts.local_shell_mask = refresh_shell;
            refresh_opts.preserve_outside_shell = true;
            [lsf_refresh, refresh_diag] = fmm_reinitialize(next_state.lsf, dx, dy, refresh_zero_mask, refresh_mask, refresh_opts);
            refresh_state = evaluate_candidate_state(lsf_refresh, theta_e, delta_theta_max, dx, dy, ...
                nelx, nely, material_mask_core, E_L, E_T, nu_LT, G_LT, thickness, F_mag, ...
                params.smooth.eta, params.smooth.iterations);
            refresh_ok_next = isfinite(refresh_state.compliance) && ...
                refresh_state.compliance <= next_state.compliance * (1 + acceptance_tol);
            if enable_reinit_current_guard && (iter >= reinit_guard_start_iter)
                refresh_ok_current = isfinite(refresh_state.compliance) && ...
                    refresh_state.compliance <= compliance * (1 + reinit_current_tol);
            else
                refresh_ok_current = true;
            end
            if refresh_ok_next && refresh_ok_current
                next_state = refresh_state;
                accepted_source = 'reinit';
                refresh_count = refresh_count + 1;
                reinit_method_history(iter) = string(refresh_diag.method_used);
                reinit_fallback_history(iter) = logical(refresh_diag.fallback_used);
                reinit_reason_history(iter) = "accepted_hj_refresh_shell";
                reinit_diagnostics_history{iter} = refresh_diag;
                last_reinit_info = refresh_diag;
                refresh_metrics = compute_raw_path_quality_metrics(lsf_refresh, dx, dy, material_mask_core, path_quality_opts);
                post_reinit_grad_dev_mean_history(iter) = refresh_metrics.grad_dev_mean;
                post_reinit_grad_outlier_ratio_history(iter) = refresh_metrics.near_zero_grad_outlier_ratio;
            end
        end
    end
    accepted_source_history{iter} = accepted_source;
    accepted_source_detail_history(iter) = string(accepted_source_detail);
    if strcmp(accepted_source, 'theta_only')
        theta_only_consecutive_count = theta_only_consecutive_count + 1;
    else
        theta_only_consecutive_count = 0;
    end
    loop_iter_count = iter;

    ENABLE_HARD_PROJECTION = logical(params.projection.enable);
    if ENABLE_HARD_PROJECTION && ~strcmp(accepted_source, 'hold')
        if iter <= 100
            omega_proj = 0.7;
            proj_band = bands.narrow_15h;
        else
            omega_proj = 0.5;
            proj_band = bands.narrow_10h;
        end
        proj_band = proj_band & material_mask_full;
        deviation_proj = next_state.lsf(proj_band) - lsf_target_global(proj_band);
        next_state.lsf(proj_band) = next_state.lsf(proj_band) - omega_proj * deviation_proj;
        next_state = evaluate_candidate_state(next_state.lsf, next_state.theta, delta_theta_max, dx, dy, ...
            nelx, nely, material_mask_core, E_L, E_T, nu_LT, G_LT, thickness, F_mag, ...
            params.smooth.eta, params.smooth.iterations);
        if iter == 1 || mod(iter, 10) == 0
            max_correction = max(abs(omega_proj * deviation_proj));
            phase_str = '是';
            if iter > 100
                phase_str = '否';
            end
            fprintf('  [硬约束投影] omega=%.2f, 最大修正=%.4f (前期=%s)\n', ...
                omega_proj, max_correction, phase_str);
        end
    end

    if strcmp(accepted_source, 'reinit')
        reinit_trigger_count = reinit_trigger_count + 1;
        iter_since_last_reinit = 0;
    elseif strcmp(accepted_source, 'hj')
        iter_since_last_reinit = iter_since_last_reinit_candidate;
    end

    lsf_change = max(abs(next_state.lsf(:) - lsf_before(:)));
    log_in = struct();
    log_in.iter = iter;
    log_in.lsf_change = lsf_change;
    log_in.compliance = compliance;
    log_in.theta_only_compliance = theta_only_compliance;
    log_in.hj_trial_compliance = hj_trial_compliance;
    log_in.reinit_trial_compliance = reinit_trial_compliance;
    log_in.accepted_source = accepted_source;
    log_in.bands_narrow_10h = bands.narrow_10h;
    log_in.material_mask_full = material_mask_full;
    log_in.next_state_lsf = next_state.lsf;
    log_in.lsf_target_global = lsf_target_global;
    log_in.next_state_FCS = next_state.FCS;
    log_in.h_grid = h_grid;
    log_in.next_state_theta = next_state.theta;
    log_in.next_state_compliance = next_state.compliance;
    fiber_log_iteration_snapshot(log_in);

    recent = compliance_history(max(1, history_count-18):history_count);
    recent = [recent; next_state.compliance];
    recent = recent(isfinite(recent) & recent > 0);
    if iter >= 20 && numel(recent) >= 10
        rel_change = std(recent) / max(mean(recent), eps);
        if rel_change < tol
            current_state = next_state;
            current_state_recorded = false;
            prev_theta_for_diag = theta_e;
            save_checkpoint(iter, current_state.lsf, current_state.theta, compliance_history(1:history_count), FCS_history(1:history_count), params);
            fprintf('优化在第 %d 次迭代时收敛\n', iter);
            break;
        end
    end

    if iter == 1
        fprintf('\n=== 系统诊断 ===\n');
        fprintf('材料各向异性 E_L/E_T = %.1f\n', E_L / E_T);
        fprintf('网格: %dx%d，单元尺寸 %.3fx%.3f\n', nelx, nely, dx, dy);
        zero_angle_elements = sum(abs(theta_e(:)) < 0.01);
        fprintf('初始接近零角度单元数: %d\n', zero_angle_elements);
        test_angle = 0;
        c = cos(test_angle);
        s = sin(test_angle);
        fprintf('当 θ=0 时: cos=%.3f，sin=%.3f，-2cs=%.3f\n', c, s, -2*c*s);
    end

    prev_theta_for_diag = theta_e;
    current_state = next_state;
    current_state_recorded = false;
    save_checkpoint(iter, current_state.lsf, current_state.theta, compliance_history(1:history_count), FCS_history(1:history_count), params);
end


    fin_in = struct();
    fin_in.params = params;
    fin_in.dx = dx;
    fin_in.dy = dy;
    fin_in.h_grid = h_grid;
    fin_in.material_mask_core = material_mask_core;
    fin_in.material_mask_full = material_mask_full;
    fin_in.boundary_guard_band = boundary_guard_band;
    fin_in.path_quality_opts = path_quality_opts;
    fin_in.current_state = current_state;
    fin_in.best_state = best_state;
    fin_in.best_state_rel_tol = best_state_rel_tol;
    fin_in.enable_best_state_guard = enable_best_state_guard;
    fin_in.compliance_history = compliance_history;
    fin_in.FCS_history = FCS_history;
    fin_in.raw_turn_history = raw_turn_history;
    fin_in.raw_kappa_history = raw_kappa_history;
    fin_in.raw_spacing_error_history = raw_spacing_error_history;
    fin_in.raw_grad_dev_history = raw_grad_dev_history;
    fin_in.raw_near_zero_outlier_history = raw_near_zero_outlier_history;
    fin_in.theta_only_compliance_history = theta_only_compliance_history;
    fin_in.hj_raw_compliance_history = hj_raw_compliance_history;
    fin_in.hj_trial_compliance_history = hj_trial_compliance_history;
    fin_in.hj_local_reinit_compliance_history = hj_local_reinit_compliance_history;
    fin_in.reinit_trial_compliance_history = reinit_trial_compliance_history;
    fin_in.accepted_source_history = accepted_source_history;
    fin_in.accepted_source_detail_history = accepted_source_detail_history;
    fin_in.hj_fallback_history = hj_fallback_history;
    fin_in.hj_second_order_history = hj_second_order_history;
    fin_in.hj_frozen_incomplete_history = hj_frozen_incomplete_history;
    fin_in.hj_first_order_complete_history = hj_first_order_complete_history;
    fin_in.reinit_method_history = reinit_method_history;
    fin_in.reinit_fallback_history = reinit_fallback_history;
    fin_in.reinit_reason_history = reinit_reason_history;
    fin_in.hj_update_diagnostics_history = hj_update_diagnostics_history;
    fin_in.reinit_diagnostics_history = reinit_diagnostics_history;
    fin_in.boundary_guard_ratio_history = boundary_guard_ratio_history;
    fin_in.frozen_boundary_point_history = frozen_boundary_point_history;
    fin_in.local_reinit_shell_size_history = local_reinit_shell_size_history;
    fin_in.post_reinit_grad_dev_mean_history = post_reinit_grad_dev_mean_history;
    fin_in.post_reinit_grad_outlier_ratio_history = post_reinit_grad_outlier_ratio_history;
    fin_in.refresh_shell_size_history = refresh_shell_size_history;
    fin_in.gradient_chain_cosine_history = gradient_chain_cosine_history;
    fin_in.gradient_chain_norm_ratio_history = gradient_chain_norm_ratio_history;
    fin_in.gradient_chain_topk_sign_history = gradient_chain_topk_sign_history;
    fin_in.gradient_chain_saturation_history = gradient_chain_saturation_history;
    fin_in.gradient_chain_degenerate_history = gradient_chain_degenerate_history;
    fin_in.gradient_chain_active_band_coverage_history = gradient_chain_active_band_coverage_history;
    fin_in.gradient_chain_zero_limiter_history = gradient_chain_zero_limiter_history;
    fin_in.gradient_chain_exact_nonzero_history = gradient_chain_exact_nonzero_history;
    fin_in.gradient_chain_support_overlap_history = gradient_chain_support_overlap_history;
    fin_in.gradient_chain_theta_raw_guard_history = gradient_chain_theta_raw_guard_history;
    fin_in.gradient_chain_theta_raw_guard_overlap_history = gradient_chain_theta_raw_guard_overlap_history;
    fin_in.gradient_chain_full_vs_opt_overlap_history = gradient_chain_full_vs_opt_overlap_history;
    fin_in.gradient_chain_selected_source_history = gradient_chain_selected_source_history;
    fin_in.gradient_chain_audit_mode = gradient_chain_audit_mode;
    fin_in.manufacturing_grad_norm_history = manufacturing_grad_norm_history;
    fin_in.manufacturing_curvature_norm_history = manufacturing_curvature_norm_history;
    fin_in.manufacturing_gap_overlap_norm_history = manufacturing_gap_overlap_norm_history;
    fin_in.theta_only_vs_current_history = theta_only_vs_current_history;
    fin_in.hj_raw_vs_theta_only_history = hj_raw_vs_theta_only_history;
    fin_in.reinit_vs_theta_only_history = reinit_vs_theta_only_history;
    fin_in.history_count = history_count;
    fin_in.current_state_recorded = current_state_recorded;
    fin_in.loop_iter_count = loop_iter_count;
    fin_in.accepted_hj_reinit_count = accepted_hj_reinit_count;
    fin_in.refresh_count = refresh_count;
    fin_in.last_hj_info = last_hj_info;
    fin_in.last_reinit_info = last_reinit_info;
    fin_in.last_velocity_field = last_velocity_field;
    fin_in.last_propagation_mask = last_propagation_mask;
    fin_in.last_gradient_chain_diag = last_gradient_chain_diag;
    fin_in.last_manufacturing_diag = last_manufacturing_diag;

    final_out = fiber_finalize_iteration_outputs(fin_in);

    iter_out = final_out;
    iter_out.accepted_steps = accepted_steps;
    iter_out.rejected_steps = rejected_steps;
    iter_out.theta_only_accept_count = theta_only_accept_count;
    iter_out.theta_only_reject_count = theta_only_reject_count;
    iter_out.theta_only_consecutive_count_final = theta_only_consecutive_count;
    iter_out.reject_due_next_guard = reject_due_next_guard;
    iter_out.reject_due_current_guard = reject_due_current_guard;
    iter_out.reinit_trigger_count = reinit_trigger_count;
    iter_out.reinit_skip_due_reject = reinit_skip_due_reject;
    iter_out.reinit_skip_due_objective = reinit_skip_due_objective;
    iter_out.reinit_skip_due_next_guard = reinit_skip_due_next_guard;
    iter_out.reinit_skip_due_current_guard = reinit_skip_due_current_guard;
    iter_out.accepted_hj_reinit_count = accepted_hj_reinit_count;
    iter_out.refresh_count = refresh_count;
    iter_out.loop_iter_count = loop_iter_count;
    iter_out.early_stop_triggered = early_stop_triggered;
    iter_out.early_stop_reason = early_stop_reason;
    iter_out.no_improve_counter = no_improve_counter;
end

function tol = compute_theta_only_tol(iter, max_iter, base_tol)
    if base_tol <= 0 || max_iter <= 0
        tol = 0;
        return;
    end

    frac = iter / max_iter;
    if frac <= 0.15
        tol = base_tol;
    elseif frac <= 0.40
        tol = base_tol * (0.40 - frac) / 0.25;
    else
        tol = 0;
    end
end

function delta = compute_candidate_delta(candidate_value, reference_value)
    if ~isfinite(candidate_value) || ~isfinite(reference_value)
        delta = NaN;
    else
        delta = candidate_value - reference_value;
    end
end
