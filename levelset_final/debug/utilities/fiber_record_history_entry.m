function history_state = fiber_record_history_entry(history_state, current_state, dx, dy, material_mask_core, path_quality_opts)
%FIBER_RECORD_HISTORY_ENTRY Record one accepted/current state into history arrays.

    if history_state.current_state_recorded
        return;
    end

    history_state.history_count = history_state.history_count + 1;
    idx = history_state.history_count;
    history_state.compliance_history(idx) = current_state.compliance;
    history_state.FCS_history(idx) = current_state.FCS;

    if idx == 1 || mod(idx - 1, history_state.quality_eval_stride) == 0
        raw_metrics_iter = compute_raw_path_quality_metrics(current_state.lsf, dx, dy, material_mask_core, path_quality_opts);
        history_state.raw_turn_history(idx) = raw_metrics_iter.mean_abs_turn_deg;
        history_state.raw_kappa_history(idx) = raw_metrics_iter.max_abs_kappa;
        history_state.raw_spacing_error_history(idx) = raw_metrics_iter.parallel_spacing_error_percent;
        history_state.raw_grad_dev_history(idx) = raw_metrics_iter.grad_dev_mean;
        history_state.raw_near_zero_outlier_history(idx) = raw_metrics_iter.near_zero_grad_outlier_ratio;
    elseif idx > 1
        history_state.raw_turn_history(idx) = history_state.raw_turn_history(idx - 1);
        history_state.raw_kappa_history(idx) = history_state.raw_kappa_history(idx - 1);
        history_state.raw_spacing_error_history(idx) = history_state.raw_spacing_error_history(idx - 1);
        history_state.raw_grad_dev_history(idx) = history_state.raw_grad_dev_history(idx - 1);
        history_state.raw_near_zero_outlier_history(idx) = history_state.raw_near_zero_outlier_history(idx - 1);
    end

    history_state.current_state_recorded = true;
end
