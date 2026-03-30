function fiber_log_iteration_snapshot(log_in)
%FIBER_LOG_ITERATION_SNAPSHOT Emit periodic per-iteration diagnostics.

    if mod(log_in.iter, 10) == 0
        fprintf('  水平集变化: %.3e\n', log_in.lsf_change);
        fprintf('  [候选诊断] C_current=%.4e, C_theta_only=%.4e, C_hj=%.4e, C_reinit=%.4e, accept=%s\n', ...
            log_in.compliance, log_in.theta_only_compliance, log_in.hj_trial_compliance, ...
            log_in.reinit_trial_compliance, log_in.accepted_source);
    end

    if mod(log_in.iter, 10) == 0
        band = log_in.bands_narrow_10h & log_in.material_mask_full;
        if any(band(:))
            dev95 = prctile(abs(log_in.next_state_lsf(band) - log_in.lsf_target_global(band)), 95);
            mean_off = mean(abs(log_in.next_state_lsf(band) - log_in.lsf_target_global(band)));
            fprintf('  ✅ [验收] dev95=%.3e  mean_off=%.3e  FCS=%.1f%%\n', ...
                dev95, mean_off, log_in.next_state_FCS*100);
            if dev95 < 1.0*log_in.h_grid && log_in.next_state_FCS >= 0.80
                fprintf('  🎉 已达标！（dev95<h, FCS≥80%%）\n');
            end
        end
    end

    if mod(log_in.iter, 5) == 0
        theta_valid = log_in.next_state_theta(isfinite(log_in.next_state_theta));
        if isempty(theta_valid)
            theta_valid = 0;
        end
        fprintf('\n迭代 %d: 柔度 = %.4e，FCS = %.2f%%\n', log_in.iter, log_in.next_state_compliance, log_in.next_state_FCS*100);
        fprintf('  角度统计：最小=%.1f 度，最大=%.1f 度，平均=%.1f 度\n', ...
            min(theta_valid(:))*180/pi, max(theta_valid(:))*180/pi, mean(theta_valid(:))*180/pi);
        fprintf('  接近 0 度单元: %d，接近 90 度单元: %d\n', ...
            sum(abs(log_in.next_state_theta(:)) < 0.1), sum(abs(log_in.next_state_theta(:) - pi/2) < 0.1));
    end
end
