function results = fiber_finalize_pipeline_results(runtime_ctx, problem_ctx, iter_out)
%FIBER_FINALIZE_PIPELINE_RESULTS Visualize, log, and build final results struct.

    params = runtime_ctx.params;
    paths = runtime_ctx.paths;

    lsf = iter_out.lsf;
    theta_e = iter_out.theta_e;
    strain_energy = iter_out.strain_energy;
    compliance_history = iter_out.compliance_history;
    FCS_history = iter_out.FCS_history;
    final_raw_metrics = iter_out.final_raw_metrics;
    path_quality_history = iter_out.path_quality_history;
    interface_diagnostics = iter_out.interface_diagnostics;

    plot_if_enabled(params, @() visualize_results_article(lsf, theta_e, strain_energy, ...
        compliance_history, FCS_history, runtime_ctx.nelx, runtime_ctx.nely, ...
        runtime_ctx.Lx, runtime_ctx.Ly, runtime_ctx.dx, runtime_ctx.dy, ...
        problem_ctx.material_mask_core, problem_ctx.material_mask_full, struct( ...
            'raw_history', path_quality_history, ...
            'path_quality_raw', final_raw_metrics, ...
            'interface_diagnostics', interface_diagnostics, ...
            'runtime_paths', paths, ...
            'material_constraint_mode', params.constraint.material_constraint_mode, ...
            'boundary_contact_policy', params.constraint.boundary_contact_policy, ...
            'init_boundary_geometry', problem_ctx.init_info.boundary_geometry, ...
            'init_phi_boundary_full', problem_ctx.init_info.phi_boundary_full)));

    log_message('INFO', params, '\n优化完成!');
    log_message('INFO', params, '最终柔度: %.4e', iter_out.compliance);
    log_message('INFO', params, '纤维连续性评分: %.2f%%', iter_out.FCS*100);
    log_message('INFO', params, 'raw路径质量: mean|turn|=%.3f deg, max|kappa|=%.3e, spacing err=%.2f%%, mean||grad|-1|=%.3e, near-zero outlier=%.2f%%', ...
        final_raw_metrics.mean_abs_turn_deg, final_raw_metrics.max_abs_kappa, ...
        final_raw_metrics.parallel_spacing_error_percent, final_raw_metrics.grad_dev_mean, ...
        100 * final_raw_metrics.near_zero_grad_outlier_ratio);
    if iter_out.early_stop_triggered
        log_message('INFO', params, '早停原因: %s', iter_out.early_stop_reason);
    end
    log_message('INFO', params, '步长接受/拒绝统计: 接受=%d, 拒绝=%d', iter_out.accepted_steps, iter_out.rejected_steps);
    if isfield(interface_diagnostics, 'gradient_chain') && isfield(interface_diagnostics, 'last_gradient_chain_diag')
        log_message('INFO', params, '梯度链模式: %s', params.gradient.chain_mode);
    end
    log_message('INFO', params, 'theta-only统计: 接受=%d, 拒绝=%d', iter_out.theta_only_accept_count, iter_out.theta_only_reject_count);
    log_message('INFO', params, '拒绝原因统计: next_guard=%d, current_guard=%d', ...
        iter_out.reject_due_next_guard, iter_out.reject_due_current_guard);
    log_message('INFO', params, ...
        ['重初始化统计: 触发=%d, 因拒绝跳过=%d, 因目标守护回退=%d, ' ...
        'next_guard回退=%d, current_guard回退=%d'], ...
        iter_out.reinit_trigger_count, iter_out.reinit_skip_due_reject, iter_out.reinit_skip_due_objective, ...
        iter_out.reinit_skip_due_next_guard, iter_out.reinit_skip_due_current_guard);
    if isfinite(iter_out.final_to_best_gap_percent)
        log_message('INFO', params, '末轮原始状态相对历史最优差距: %.4f%%', iter_out.final_to_best_gap_percent);
    end

    base_compliance = compliance_history(1);
    current_compliance = iter_out.compliance;
    if base_compliance > 0
        improve_ratio = (base_compliance - current_compliance) / base_compliance * 100;
        log_message('INFO', params, '柔度降低比例: %.2f%% (初始=%.4e, 最终=%.4e)', ...
            improve_ratio, base_compliance, current_compliance);
    else
        log_message('WARN', params, '警告：初始柔度<=0，无法计算降低比例');
    end

    if base_compliance > 0
        improvement_ratio = improve_ratio;
    else
        improvement_ratio = NaN;
    end

    results_input = struct();
    results_input.lsf = lsf;
    results_input.theta_e = theta_e;
    results_input.compliance_history = compliance_history;
    results_input.FCS_history = FCS_history;
    results_input.final_compliance = iter_out.compliance;
    results_input.final_FCS = iter_out.FCS;
    results_input.final_iter = iter_out.final_iter;
    results_input.executed_iter = iter_out.executed_iter;
    results_input.strain_energy = strain_energy;
    results_input.params = params;
    results_input.accepted_steps = iter_out.accepted_steps;
    results_input.rejected_steps = iter_out.rejected_steps;
    results_input.reject_due_next_guard = iter_out.reject_due_next_guard;
    results_input.reject_due_current_guard = iter_out.reject_due_current_guard;
    results_input.reinit_trigger_count = iter_out.reinit_trigger_count;
    results_input.reinit_skip_due_reject = iter_out.reinit_skip_due_reject;
    results_input.reinit_skip_due_objective = iter_out.reinit_skip_due_objective;
    results_input.reinit_skip_due_next_guard = iter_out.reinit_skip_due_next_guard;
    results_input.reinit_skip_due_current_guard = iter_out.reinit_skip_due_current_guard;
    results_input.accepted_hj_reinit_count = iter_out.accepted_hj_reinit_count;
    results_input.refresh_count = iter_out.refresh_count;
    results_input.theta_target = iter_out.theta_target;
    results_input.theta_only_compliance_history = iter_out.theta_only_compliance_history;
    results_input.hj_raw_compliance_history = iter_out.hj_raw_compliance_history;
    results_input.hj_trial_compliance_history = iter_out.hj_trial_compliance_history;
    results_input.hj_local_reinit_compliance_history = iter_out.hj_local_reinit_compliance_history;
    results_input.reinit_trial_compliance_history = iter_out.reinit_trial_compliance_history;
    results_input.accepted_source_history = iter_out.accepted_source_history;
    results_input.accepted_source_detail_history = iter_out.accepted_source_detail_history;
    results_input.theta_only_accept_count = iter_out.theta_only_accept_count;
    results_input.theta_only_reject_count = iter_out.theta_only_reject_count;
    results_input.loop_iter_count = iter_out.loop_iter_count;
    results_input.best_iter = iter_out.best_state.iter;
    results_input.best_compliance = iter_out.best_state.compliance;
    results_input.best_FCS = iter_out.best_state.FCS;
    results_input.rollback_to_best = iter_out.rollback_to_best;
    results_input.early_stop_triggered = iter_out.early_stop_triggered;
    results_input.early_stop_reason = iter_out.early_stop_reason;
    results_input.no_improve_counter_final = iter_out.no_improve_counter;
    results_input.raw_final_compliance = iter_out.raw_final_compliance;
    results_input.raw_final_FCS = iter_out.raw_final_FCS;
    results_input.raw_final_improvement_ratio = iter_out.raw_final_improvement_ratio;
    results_input.final_to_best_gap_percent = iter_out.final_to_best_gap_percent;
    results_input.theta_only_vs_current_history = iter_out.theta_only_vs_current_history;
    results_input.hj_raw_vs_theta_only_history = iter_out.hj_raw_vs_theta_only_history;
    results_input.reinit_vs_theta_only_history = iter_out.reinit_vs_theta_only_history;
    results_input.material_mask_core = problem_ctx.material_mask_core;
    results_input.material_mask_full = problem_ctx.material_mask_full;
    results_input.path_quality_raw = final_raw_metrics;
    results_input.path_quality_history = path_quality_history;
    results_input.interface_diagnostics = interface_diagnostics;
    results_input.gradient_chain_history = interface_diagnostics.gradient_chain;
    results_input.constraint = params.constraint;
    results_input.init_boundary_geometry = problem_ctx.init_info.boundary_geometry;
    results_input.init_info = problem_ctx.init_info;
    results_input.improvement_ratio = improvement_ratio;
    results = build_results_struct(results_input);

    log_message('INFO', params, '结果已保存至输出结构体');
end
