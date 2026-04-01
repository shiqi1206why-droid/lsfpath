function out = fiber_finalize_iteration_outputs(fin_in)
%FIBER_FINALIZE_ITERATION_OUTPUTS Consolidate histories, rollback, and diagnostics after loop.

    params = fin_in.params;
    dx = fin_in.dx;
    dy = fin_in.dy;
    h_grid = fin_in.h_grid;
    material_mask_core = fin_in.material_mask_core;
    material_mask_full = fin_in.material_mask_full;
    boundary_guard_band = fin_in.boundary_guard_band;
    path_quality_opts = fin_in.path_quality_opts;

    current_state = fin_in.current_state;
    best_state = fin_in.best_state;
    best_state_rel_tol = fin_in.best_state_rel_tol;

    compliance_history = fin_in.compliance_history;
    FCS_history = fin_in.FCS_history;
    raw_turn_history = fin_in.raw_turn_history;
    raw_kappa_history = fin_in.raw_kappa_history;
    raw_spacing_error_history = fin_in.raw_spacing_error_history;
    raw_grad_dev_history = fin_in.raw_grad_dev_history;
    raw_near_zero_outlier_history = fin_in.raw_near_zero_outlier_history;
    theta_only_compliance_history = fin_in.theta_only_compliance_history;
    hj_raw_compliance_history = fin_in.hj_raw_compliance_history;
    hj_trial_compliance_history = fin_in.hj_trial_compliance_history;
    hj_local_reinit_compliance_history = fin_in.hj_local_reinit_compliance_history;
    reinit_trial_compliance_history = fin_in.reinit_trial_compliance_history;
    accepted_source_history = fin_in.accepted_source_history;
    accepted_source_detail_history = fin_in.accepted_source_detail_history;
    hj_fallback_history = fin_in.hj_fallback_history;
    hj_second_order_history = fin_in.hj_second_order_history;
    hj_frozen_incomplete_history = fin_in.hj_frozen_incomplete_history;
    hj_first_order_complete_history = fin_in.hj_first_order_complete_history;
    reinit_method_history = fin_in.reinit_method_history;
    reinit_fallback_history = fin_in.reinit_fallback_history;
    reinit_reason_history = fin_in.reinit_reason_history;
    hj_update_diagnostics_history = fin_in.hj_update_diagnostics_history;
    reinit_diagnostics_history = fin_in.reinit_diagnostics_history;
    boundary_guard_ratio_history = fin_in.boundary_guard_ratio_history;
    frozen_boundary_point_history = fin_in.frozen_boundary_point_history;
    local_reinit_shell_size_history = fin_in.local_reinit_shell_size_history;
    post_reinit_grad_dev_mean_history = fin_in.post_reinit_grad_dev_mean_history;
    post_reinit_grad_outlier_ratio_history = fin_in.post_reinit_grad_outlier_ratio_history;
    refresh_shell_size_history = fin_in.refresh_shell_size_history;
    gradient_chain_cosine_history = fin_in.gradient_chain_cosine_history;
    gradient_chain_norm_ratio_history = fin_in.gradient_chain_norm_ratio_history;
    gradient_chain_topk_sign_history = fin_in.gradient_chain_topk_sign_history;
    gradient_chain_saturation_history = fin_in.gradient_chain_saturation_history;
    gradient_chain_degenerate_history = fin_in.gradient_chain_degenerate_history;
    gradient_chain_active_band_coverage_history = fin_in.gradient_chain_active_band_coverage_history;
    gradient_chain_zero_limiter_history = fin_in.gradient_chain_zero_limiter_history;
    gradient_chain_exact_nonzero_history = fin_in.gradient_chain_exact_nonzero_history;
    gradient_chain_support_overlap_history = fin_in.gradient_chain_support_overlap_history;
    gradient_chain_theta_raw_guard_history = fin_in.gradient_chain_theta_raw_guard_history;
    gradient_chain_theta_raw_guard_overlap_history = fin_in.gradient_chain_theta_raw_guard_overlap_history;
    gradient_chain_full_vs_opt_overlap_history = fin_in.gradient_chain_full_vs_opt_overlap_history;
    gradient_chain_selected_source_history = fin_in.gradient_chain_selected_source_history;
    gradient_chain_audit_mode = fin_in.gradient_chain_audit_mode;
    manufacturing_grad_norm_history = fin_in.manufacturing_grad_norm_history;
    manufacturing_curvature_norm_history = fin_in.manufacturing_curvature_norm_history;
    manufacturing_gap_overlap_norm_history = fin_in.manufacturing_gap_overlap_norm_history;
    theta_only_vs_current_history = fin_in.theta_only_vs_current_history;
    hj_raw_vs_theta_only_history = fin_in.hj_raw_vs_theta_only_history;
    reinit_vs_theta_only_history = fin_in.reinit_vs_theta_only_history;

    history_count = fin_in.history_count;
    current_state_recorded = fin_in.current_state_recorded;

    if ~current_state_recorded
        history_count = history_count + 1;
        compliance_history(history_count) = current_state.compliance;
        FCS_history(history_count) = current_state.FCS;
        raw_metrics_iter = compute_raw_path_quality_metrics(current_state.lsf, dx, dy, material_mask_core, path_quality_opts);
        raw_turn_history(history_count) = raw_metrics_iter.mean_abs_turn_deg;
        raw_kappa_history(history_count) = raw_metrics_iter.max_abs_kappa;
        raw_spacing_error_history(history_count) = raw_metrics_iter.parallel_spacing_error_percent;
        raw_grad_dev_history(history_count) = raw_metrics_iter.grad_dev_mean;
        raw_near_zero_outlier_history(history_count) = raw_metrics_iter.near_zero_grad_outlier_ratio;
        if current_state.compliance < best_state.compliance * (1 - best_state_rel_tol)
            best_state.compliance = current_state.compliance;
            best_state.iter = history_count;
            best_state.lsf = current_state.lsf;
            best_state.theta = current_state.theta;
            best_state.strain_energy = current_state.strain_energy;
            best_state.FCS = current_state.FCS;
            best_state.theta_target = current_state.theta_target;
        end
    end

    final_iter = history_count;
    history_input = struct();
    history_input.compliance_history = compliance_history;
    history_input.FCS_history = FCS_history;
    history_input.raw_turn_history = raw_turn_history;
    history_input.raw_kappa_history = raw_kappa_history;
    history_input.raw_spacing_error_history = raw_spacing_error_history;
    history_input.raw_grad_dev_history = raw_grad_dev_history;
    history_input.raw_near_zero_outlier_history = raw_near_zero_outlier_history;
    history_input.theta_only_compliance_history = theta_only_compliance_history;
    history_input.hj_raw_compliance_history = hj_raw_compliance_history;
    history_input.hj_trial_compliance_history = hj_trial_compliance_history;
    history_input.hj_local_reinit_compliance_history = hj_local_reinit_compliance_history;
    history_input.reinit_trial_compliance_history = reinit_trial_compliance_history;
    history_input.accepted_source_history = accepted_source_history;
    history_input.accepted_source_detail_history = accepted_source_detail_history;
    history_input.hj_fallback_history = hj_fallback_history;
    history_input.hj_second_order_history = hj_second_order_history;
    history_input.hj_frozen_incomplete_history = hj_frozen_incomplete_history;
    history_input.hj_first_order_complete_history = hj_first_order_complete_history;
    history_input.reinit_method_history = reinit_method_history;
    history_input.reinit_fallback_history = reinit_fallback_history;
    history_input.reinit_reason_history = reinit_reason_history;
    history_input.hj_update_diagnostics_history = hj_update_diagnostics_history;
    history_input.reinit_diagnostics_history = reinit_diagnostics_history;
    history_input.boundary_guard_ratio_history = boundary_guard_ratio_history;
    history_input.frozen_boundary_point_history = frozen_boundary_point_history;
    history_input.local_reinit_shell_size_history = local_reinit_shell_size_history;
    history_input.post_reinit_grad_dev_mean_history = post_reinit_grad_dev_mean_history;
    history_input.post_reinit_grad_outlier_ratio_history = post_reinit_grad_outlier_ratio_history;
    history_input.refresh_shell_size_history = refresh_shell_size_history;
    history_input.gradient_chain_cosine_history = gradient_chain_cosine_history;
    history_input.gradient_chain_norm_ratio_history = gradient_chain_norm_ratio_history;
    history_input.gradient_chain_topk_sign_history = gradient_chain_topk_sign_history;
    history_input.gradient_chain_saturation_history = gradient_chain_saturation_history;
    history_input.gradient_chain_degenerate_history = gradient_chain_degenerate_history;
    history_input.gradient_chain_active_band_coverage_history = gradient_chain_active_band_coverage_history;
    history_input.gradient_chain_zero_limiter_history = gradient_chain_zero_limiter_history;
    history_input.gradient_chain_exact_nonzero_history = gradient_chain_exact_nonzero_history;
    history_input.gradient_chain_support_overlap_history = gradient_chain_support_overlap_history;
    history_input.gradient_chain_theta_raw_guard_history = gradient_chain_theta_raw_guard_history;
    history_input.gradient_chain_theta_raw_guard_overlap_history = gradient_chain_theta_raw_guard_overlap_history;
    history_input.gradient_chain_full_vs_opt_overlap_history = gradient_chain_full_vs_opt_overlap_history;
    history_input.gradient_chain_selected_source_history = gradient_chain_selected_source_history;
    history_input.manufacturing_grad_norm_history = manufacturing_grad_norm_history;
    history_input.manufacturing_curvature_norm_history = manufacturing_curvature_norm_history;
    history_input.manufacturing_gap_overlap_norm_history = manufacturing_gap_overlap_norm_history;
    history_input.theta_only_vs_current_history = theta_only_vs_current_history;
    history_input.hj_raw_vs_theta_only_history = hj_raw_vs_theta_only_history;
    history_input.reinit_vs_theta_only_history = reinit_vs_theta_only_history;
    history_data = finalize_history_data(final_iter, fin_in.loop_iter_count, history_input);

    compliance_history = history_data.compliance_history;
    FCS_history = history_data.FCS_history;
    raw_turn_history = history_data.raw_turn_history;
    raw_kappa_history = history_data.raw_kappa_history;
    raw_spacing_error_history = history_data.raw_spacing_error_history;
    raw_grad_dev_history = history_data.raw_grad_dev_history;
    raw_near_zero_outlier_history = history_data.raw_near_zero_outlier_history;
    theta_only_compliance_history = history_data.theta_only_compliance_history;
    hj_raw_compliance_history = history_data.hj_raw_compliance_history;
    hj_trial_compliance_history = history_data.hj_trial_compliance_history;
    hj_local_reinit_compliance_history = history_data.hj_local_reinit_compliance_history;
    reinit_trial_compliance_history = history_data.reinit_trial_compliance_history;
    accepted_source_history = history_data.accepted_source_history;
    accepted_source_detail_history = history_data.accepted_source_detail_history;
    hj_fallback_history = history_data.hj_fallback_history;
    hj_second_order_history = history_data.hj_second_order_history;
    hj_frozen_incomplete_history = history_data.hj_frozen_incomplete_history;
    hj_first_order_complete_history = history_data.hj_first_order_complete_history;
    reinit_method_history = history_data.reinit_method_history;
    reinit_fallback_history = history_data.reinit_fallback_history;
    reinit_reason_history = history_data.reinit_reason_history;
    hj_update_diagnostics_history = history_data.hj_update_diagnostics_history;
    reinit_diagnostics_history = history_data.reinit_diagnostics_history;
    boundary_guard_ratio_history = history_data.boundary_guard_ratio_history;
    frozen_boundary_point_history = history_data.frozen_boundary_point_history;
    local_reinit_shell_size_history = history_data.local_reinit_shell_size_history;
    post_reinit_grad_dev_mean_history = history_data.post_reinit_grad_dev_mean_history;
    post_reinit_grad_outlier_ratio_history = history_data.post_reinit_grad_outlier_ratio_history;
    refresh_shell_size_history = history_data.refresh_shell_size_history;
    gradient_chain_cosine_history = history_data.gradient_chain_cosine_history;
    gradient_chain_norm_ratio_history = history_data.gradient_chain_norm_ratio_history;
    gradient_chain_topk_sign_history = history_data.gradient_chain_topk_sign_history;
    gradient_chain_saturation_history = history_data.gradient_chain_saturation_history;
    gradient_chain_degenerate_history = history_data.gradient_chain_degenerate_history;
    gradient_chain_active_band_coverage_history = history_data.gradient_chain_active_band_coverage_history;
    gradient_chain_zero_limiter_history = history_data.gradient_chain_zero_limiter_history;
    gradient_chain_exact_nonzero_history = history_data.gradient_chain_exact_nonzero_history;
    gradient_chain_support_overlap_history = history_data.gradient_chain_support_overlap_history;
    gradient_chain_theta_raw_guard_history = history_data.gradient_chain_theta_raw_guard_history;
    gradient_chain_theta_raw_guard_overlap_history = history_data.gradient_chain_theta_raw_guard_overlap_history;
    gradient_chain_full_vs_opt_overlap_history = history_data.gradient_chain_full_vs_opt_overlap_history;
    gradient_chain_selected_source_history = history_data.gradient_chain_selected_source_history;
    manufacturing_grad_norm_history = history_data.manufacturing_grad_norm_history;
    manufacturing_curvature_norm_history = history_data.manufacturing_curvature_norm_history;
    manufacturing_gap_overlap_norm_history = history_data.manufacturing_gap_overlap_norm_history;
    theta_only_vs_current_history = history_data.theta_only_vs_current_history;
    hj_raw_vs_theta_only_history = history_data.hj_raw_vs_theta_only_history;
    reinit_vs_theta_only_history = history_data.reinit_vs_theta_only_history;
    executed_iter = history_data.executed_iter;

    raw_final_compliance = current_state.compliance;
    raw_final_FCS = current_state.FCS;
    raw_final_improvement_ratio = NaN;
    final_to_best_gap_percent = NaN;
    if final_iter >= 1 && compliance_history(1) > 0
        raw_final_improvement_ratio = (compliance_history(1) - raw_final_compliance) / compliance_history(1) * 100;
    end
    if isfinite(best_state.compliance) && best_state.compliance > 0
        final_to_best_gap_percent = (raw_final_compliance - best_state.compliance) / best_state.compliance * 100;
    end

    rollback_to_best = false;
    if fin_in.enable_best_state_guard && best_state.iter > 0
        raw_final_worse_than_best = ~isfinite(raw_final_compliance) || ...
            raw_final_compliance > best_state.compliance * (1 + best_state_rel_tol);
        if raw_final_worse_than_best
            rollback_to_best = true;
            log_message('INFO', params, ...
                '输出状态回滚到历史最优: best_iter=%d/%d, C_best=%.4e, C_last=%.4e', ...
                best_state.iter, executed_iter, best_state.compliance, compliance_history(end));
            lsf = best_state.lsf;
            theta_e = best_state.theta;
            strain_energy = best_state.strain_energy;
            compliance = best_state.compliance;
            FCS = best_state.FCS;
            theta_target = best_state.theta_target;
        else
            lsf = current_state.lsf;
            theta_e = current_state.theta;
            strain_energy = current_state.strain_energy;
            compliance = raw_final_compliance;
            FCS = raw_final_FCS;
            theta_target = current_state.theta_target;
        end
    else
        if isempty(FCS_history)
            FCS = NaN;
        else
            FCS = FCS_history(end);
        end
        lsf = current_state.lsf;
        theta_e = current_state.theta;
        strain_energy = current_state.strain_energy;
        compliance = raw_final_compliance;
        theta_target = current_state.theta_target;
    end

    final_raw_metrics = compute_raw_path_quality_metrics(lsf, dx, dy, material_mask_core, path_quality_opts);
    path_quality_history = build_path_quality_history( ...
        raw_turn_history, raw_kappa_history, raw_spacing_error_history, ...
        raw_grad_dev_history, raw_near_zero_outlier_history);

    interface_input = struct();
    interface_input.params = params;
    interface_input.path_quality_history = path_quality_history;
    interface_input.hj_fallback_history = hj_fallback_history;
    interface_input.hj_second_order_history = hj_second_order_history;
    interface_input.hj_frozen_incomplete_history = hj_frozen_incomplete_history;
    interface_input.hj_first_order_complete_history = hj_first_order_complete_history;
    interface_input.reinit_method_history = reinit_method_history;
    interface_input.reinit_fallback_history = reinit_fallback_history;
    interface_input.reinit_reason_history = reinit_reason_history;
    interface_input.hj_update_diagnostics_history = hj_update_diagnostics_history;
    interface_input.reinit_diagnostics_history = reinit_diagnostics_history;
    interface_input.boundary_guard_ratio_history = boundary_guard_ratio_history;
    interface_input.frozen_boundary_point_history = frozen_boundary_point_history;
    interface_input.accepted_hj_reinit_count = fin_in.accepted_hj_reinit_count;
    interface_input.local_reinit_shell_size_history = local_reinit_shell_size_history;
    interface_input.post_reinit_grad_dev_mean_history = post_reinit_grad_dev_mean_history;
    interface_input.post_reinit_grad_outlier_ratio_history = post_reinit_grad_outlier_ratio_history;
    interface_input.refresh_shell_size_history = refresh_shell_size_history;
    interface_input.refresh_count = fin_in.refresh_count;
    interface_input.last_hj_info = fin_in.last_hj_info;
    interface_input.last_reinit_info = fin_in.last_reinit_info;
    interface_input.last_velocity_field = fin_in.last_velocity_field;
    interface_input.last_propagation_mask = fin_in.last_propagation_mask;
    interface_input.last_gradient_chain_diag = fin_in.last_gradient_chain_diag;
    interface_input.last_manufacturing_diag = fin_in.last_manufacturing_diag;
    interface_input.gradient_chain_cosine_history = gradient_chain_cosine_history;
    interface_input.gradient_chain_norm_ratio_history = gradient_chain_norm_ratio_history;
    interface_input.gradient_chain_topk_sign_history = gradient_chain_topk_sign_history;
    interface_input.gradient_chain_saturation_history = gradient_chain_saturation_history;
    interface_input.gradient_chain_degenerate_history = gradient_chain_degenerate_history;
    interface_input.gradient_chain_active_band_coverage_history = gradient_chain_active_band_coverage_history;
    interface_input.gradient_chain_zero_limiter_history = gradient_chain_zero_limiter_history;
    interface_input.gradient_chain_exact_nonzero_history = gradient_chain_exact_nonzero_history;
    interface_input.gradient_chain_support_overlap_history = gradient_chain_support_overlap_history;
    interface_input.gradient_chain_theta_raw_guard_history = gradient_chain_theta_raw_guard_history;
    interface_input.gradient_chain_theta_raw_guard_overlap_history = gradient_chain_theta_raw_guard_overlap_history;
    interface_input.gradient_chain_full_vs_opt_overlap_history = gradient_chain_full_vs_opt_overlap_history;
    interface_input.gradient_chain_selected_source_history = gradient_chain_selected_source_history;
    interface_input.gradient_chain_audit_mode = gradient_chain_audit_mode;
    interface_input.manufacturing_grad_norm_history = manufacturing_grad_norm_history;
    interface_input.manufacturing_curvature_norm_history = manufacturing_curvature_norm_history;
    interface_input.manufacturing_gap_overlap_norm_history = manufacturing_gap_overlap_norm_history;
    interface_input.theta_only_vs_current_history = theta_only_vs_current_history;
    interface_input.hj_raw_vs_theta_only_history = hj_raw_vs_theta_only_history;
    interface_input.reinit_vs_theta_only_history = reinit_vs_theta_only_history;
    interface_input.boundary_guard_band = boundary_guard_band;
    interface_input.material_mask_full = material_mask_full;
    interface_input.lsf = lsf;
    interface_diagnostics = build_interface_diagnostics(interface_input);

    out = struct();
    out.current_state = current_state;
    out.best_state = best_state;
    out.final_iter = final_iter;
    out.executed_iter = executed_iter;

    out.compliance_history = compliance_history;
    out.FCS_history = FCS_history;
    out.raw_turn_history = raw_turn_history;
    out.raw_kappa_history = raw_kappa_history;
    out.raw_spacing_error_history = raw_spacing_error_history;
    out.raw_grad_dev_history = raw_grad_dev_history;
    out.raw_near_zero_outlier_history = raw_near_zero_outlier_history;
    out.theta_only_compliance_history = theta_only_compliance_history;
    out.hj_raw_compliance_history = hj_raw_compliance_history;
    out.hj_trial_compliance_history = hj_trial_compliance_history;
    out.hj_local_reinit_compliance_history = hj_local_reinit_compliance_history;
    out.reinit_trial_compliance_history = reinit_trial_compliance_history;
    out.accepted_source_history = accepted_source_history;
    out.accepted_source_detail_history = accepted_source_detail_history;
    out.hj_fallback_history = hj_fallback_history;
    out.hj_second_order_history = hj_second_order_history;
    out.hj_frozen_incomplete_history = hj_frozen_incomplete_history;
    out.hj_first_order_complete_history = hj_first_order_complete_history;
    out.reinit_method_history = reinit_method_history;
    out.reinit_fallback_history = reinit_fallback_history;
    out.reinit_reason_history = reinit_reason_history;
    out.hj_update_diagnostics_history = hj_update_diagnostics_history;
    out.reinit_diagnostics_history = reinit_diagnostics_history;
    out.boundary_guard_ratio_history = boundary_guard_ratio_history;
    out.frozen_boundary_point_history = frozen_boundary_point_history;
    out.local_reinit_shell_size_history = local_reinit_shell_size_history;
    out.post_reinit_grad_dev_mean_history = post_reinit_grad_dev_mean_history;
    out.post_reinit_grad_outlier_ratio_history = post_reinit_grad_outlier_ratio_history;
    out.refresh_shell_size_history = refresh_shell_size_history;
    out.gradient_chain_cosine_history = gradient_chain_cosine_history;
    out.gradient_chain_norm_ratio_history = gradient_chain_norm_ratio_history;
    out.gradient_chain_topk_sign_history = gradient_chain_topk_sign_history;
    out.gradient_chain_saturation_history = gradient_chain_saturation_history;
    out.gradient_chain_degenerate_history = gradient_chain_degenerate_history;
    out.gradient_chain_active_band_coverage_history = gradient_chain_active_band_coverage_history;
    out.gradient_chain_zero_limiter_history = gradient_chain_zero_limiter_history;
    out.gradient_chain_exact_nonzero_history = gradient_chain_exact_nonzero_history;
    out.gradient_chain_support_overlap_history = gradient_chain_support_overlap_history;
    out.gradient_chain_theta_raw_guard_history = gradient_chain_theta_raw_guard_history;
    out.gradient_chain_theta_raw_guard_overlap_history = gradient_chain_theta_raw_guard_overlap_history;
    out.gradient_chain_full_vs_opt_overlap_history = gradient_chain_full_vs_opt_overlap_history;
    out.gradient_chain_selected_source_history = gradient_chain_selected_source_history;
    out.gradient_chain_audit_mode = gradient_chain_audit_mode;
    out.manufacturing_grad_norm_history = manufacturing_grad_norm_history;
    out.manufacturing_curvature_norm_history = manufacturing_curvature_norm_history;
    out.manufacturing_gap_overlap_norm_history = manufacturing_gap_overlap_norm_history;
    out.theta_only_vs_current_history = theta_only_vs_current_history;
    out.hj_raw_vs_theta_only_history = hj_raw_vs_theta_only_history;
    out.reinit_vs_theta_only_history = reinit_vs_theta_only_history;

    out.lsf = lsf;
    out.theta_e = theta_e;
    out.strain_energy = strain_energy;
    out.compliance = compliance;
    out.FCS = FCS;
    out.theta_target = theta_target;

    out.final_raw_metrics = final_raw_metrics;
    out.path_quality_history = path_quality_history;
    out.interface_diagnostics = interface_diagnostics;

    out.rollback_to_best = rollback_to_best;
    out.raw_final_compliance = raw_final_compliance;
    out.raw_final_FCS = raw_final_FCS;
    out.raw_final_improvement_ratio = raw_final_improvement_ratio;
    out.final_to_best_gap_percent = final_to_best_gap_percent;
end
