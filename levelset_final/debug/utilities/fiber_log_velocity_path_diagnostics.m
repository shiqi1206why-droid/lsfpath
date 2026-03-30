function fiber_log_velocity_path_diagnostics(log_in)
%FIBER_LOG_VELOCITY_PATH_DIAGNOSTICS Emit gradient/path diagnostics for one iteration.

    if ~(log_in.iter == 1 || mod(log_in.iter, 10) == 0)
        return;
    end

    [grad_y, grad_x] = gradient(log_in.lsf, log_in.dy, log_in.dx);
    grad_mag = hypot(grad_x, grad_y);
    deviation = abs(log_in.lsf - log_in.lsf_target_global);
    raw_metrics_diag = compute_raw_path_quality_metrics(log_in.lsf, log_in.dx, log_in.dy, ...
        log_in.material_mask_core, log_in.path_quality_opts);

    fprintf('  [诊断] 梯度模：最小=%.3e，最大=%.3e，平均=%.3e\n', ...
        min(grad_mag(:)), max(grad_mag(:)), mean(grad_mag(:)));
    fprintf('  [诊断] 偏离量：最小=%.3e，最大=%.3e，平均=%.3e\n', ...
        min(deviation(:)), max(deviation(:)), mean(deviation(:)));
    fprintf('  [诊断] node_sensitivity：最小=%.2e，最大=%.2e，平均=%.2e，scale=%.2e\n', ...
        min(log_in.node_sensitivity(:)), max(log_in.node_sensitivity(:)), ...
        mean(abs(log_in.node_sensitivity(:))), log_in.sens_scale_value);
    fprintf('  [诊断] 速度场max_band=%.2e\n', log_in.velocity_stats.max_band);

    zero_band = log_in.bands.narrow_05h & log_in.material_mask_full;
    if any(zero_band(:))
        deviation_zero = deviation(zero_band);
        consistency_ratio = sum(deviation_zero < log_in.h) / numel(deviation_zero) * 100;
        fprintf('  [诊断] 路径一致性：零线附近%.1f%%节点偏离<1网格（平均偏离=%.3e）\n', ...
            consistency_ratio, mean(deviation_zero));
    end

    fprintf('  [raw路径] mean|turn|=%.3f deg, max|kappa|=%.3e, spacing err=%.2f%%, mean||grad|-1|=%.3e, near-zero outlier=%.2f%%\n', ...
        raw_metrics_diag.mean_abs_turn_deg, raw_metrics_diag.max_abs_kappa, ...
        raw_metrics_diag.parallel_spacing_error_percent, raw_metrics_diag.grad_dev_mean, ...
        100 * raw_metrics_diag.near_zero_grad_outlier_ratio);
end
