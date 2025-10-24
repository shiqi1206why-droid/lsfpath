function [should_reinit, reason] = should_reinitialize(lsf_for_change, lsf_for_gradient, lsf_before, ...
    iter_since_last_reinit, iter, params)
    % 自适应重初始化判断
    % 
    % 输入：
    %   lsf_for_change   - 用于变化量判据的水平集（通常是HJ演化后、投影前）
    %   lsf_for_gradient - 用于梯度偏差判据的水平集（通常是投影后，最终使用的状态）
    %   lsf_before       - 更新前水平集
    %   iter_since_last_reinit - 距上次重初始化的迭代数
    %   iter - 当前迭代次数
    %   params - 参数结构
    % 
    % 输出：
    %   should_reinit - 是否需要重初始化（逻辑值）
    %   reason - 触发原因（字符串，用于日志）
    %
    % 判据（按优先级）：
    %   1. 水平集变化过大（基于lsf_for_change，排除投影影响）
    %   2. 梯度模偏离1过多（基于lsf_for_gradient，使用最终状态）
    %   3. 超过最大允许间隔
    %   4. 前期/后期固定频率（保持原有逻辑）
    
    should_reinit = false;
    reason = '';
    
    h = params.grid.h;
    dx = params.grid.dx;
    dy = params.grid.dy;
    
    % 判据1：水平集变化过大（监控HJ演化，不包括投影修正）
    lsf_change = max(abs(lsf_for_change(:) - lsf_before(:)));
    threshold_change = params.levelset.reinit_threshold * h;
    if lsf_change > threshold_change
        should_reinit = true;
        reason = sprintf('水平集变化过大 (%.3e > %.3e)', lsf_change, threshold_change);
        return;
    end
    
    % 判据2：梯度模偏离1过多（符号距离性质退化，仅在窄带统计）
    % 设计意图：只关心主路径附近的梯度质量，远离主路径的区域梯度质量本就较差
    % 使用投影后的lsf计算梯度，因为这是最终使用的状态
    [grad_y, grad_x] = gradient(lsf_for_gradient, dy, dx);
    grad_mag = hypot(grad_x, grad_y);
    
    % 窄带统计（默认1.5h，可配置）
    reinit_grad_bandwidth = 1.5;
    if isfield(params.levelset, 'reinit_grad_bandwidth')
        reinit_grad_bandwidth = params.levelset.reinit_grad_bandwidth;
    end
    % 基于投影后的lsf确定窄带
    band = abs(lsf_for_gradient) <= reinit_grad_bandwidth * h;
    
    if any(band(:))
        % 使用P95统计替代均值，更抗噪，忽略极值点
        grad_deviation = prctile(abs(grad_mag(band) - 1), 95);
        band_info = sprintf('窄带%.1fh,P95', reinit_grad_bandwidth);
    else
        % 兜底：窄带为空时使用全域
        grad_deviation = prctile(abs(grad_mag(:) - 1), 95);
        band_info = '全域,P95';
    end
    
    grad_tol = params.levelset.gradient_deviation_tol;
    if grad_deviation > grad_tol
        should_reinit = true;
        reason = sprintf('梯度偏差过大 (%.3f > %.3f, %s)', grad_deviation, grad_tol, band_info);
        return;
    end
    
    % 判据3：超过最大允许间隔
    if iter_since_last_reinit >= params.levelset.reinit_max_interval
        should_reinit = true;
        reason = sprintf('达到最大间隔 (%d)', params.levelset.reinit_max_interval);
        return;
    end
    
    % 判据4：前期/后期固定频率（保持原有逻辑）
    is_early_phase = (iter <= params.levelset.transition_iter);
    if is_early_phase && iter_since_last_reinit >= params.levelset.reinit_freq_early
        should_reinit = true;
        reason = sprintf('前期固定频率 (%d次)', params.levelset.reinit_freq_early);
        return;
    elseif ~is_early_phase && iter_since_last_reinit >= params.levelset.reinit_freq_late
        should_reinit = true;
        reason = sprintf('后期固定频率 (%d次)', params.levelset.reinit_freq_late);
        return;
    end
end
