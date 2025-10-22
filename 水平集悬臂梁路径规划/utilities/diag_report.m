function diag_report(iter, node_sensitivity, velocity_stats)
    % 输出诊断信息的辅助函数

    global DIAG;
    try
        nz = nnz(node_sensitivity);
        total = numel(node_sensitivity);
        fprintf('  诊断：节点灵敏度非零比例 = %.1f%%\n', 100*nz/max(1,total));

        if DIAG.theta_grad_samples > 0
            fprintf('  诊断：θ 梯度接近零比例 = %.1f%% (%d/%d)\n', ...
                100*DIAG.theta_grad_small/DIAG.theta_grad_samples, ...
                DIAG.theta_grad_small, DIAG.theta_grad_samples);
        end

        if DIAG.den_small > 0
            fprintf('  诊断：链式法分母退化次数 = %d\n', DIAG.den_small);
        end

        if DIAG.contrib_total > 0
            fprintf('  诊断：链式法非零贡献比例 = %.1f%% (%d/%d)\n', ...
                100*DIAG.contrib_nonzero/DIAG.contrib_total, ...
                DIAG.contrib_nonzero, DIAG.contrib_total);
        end

        if isstruct(velocity_stats)
            fprintf('  诊断：窄带覆盖=%.2f%%，窄带平均=%.3e，加权平均=%.3e\n', ...
                100*velocity_stats.band_fraction, velocity_stats.mean_band, velocity_stats.weighted_mean);
        end
    catch ME
        fprintf('  诊断输出错误: %s\n', ME.message);
    end
end

