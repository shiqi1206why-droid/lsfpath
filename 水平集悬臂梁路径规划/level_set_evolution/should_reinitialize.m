function [should_reinit, reason] = should_reinitialize(lsf, lsf_before, ...
    iter_since_last_reinit, iter, params)
    % 自适应重初始化判断
    % 
    % 输入：
    %   lsf - 当前水平集 [nely+2, nelx+2]
    %   lsf_before - 更新前水平集
    %   iter_since_last_reinit - 距上次重初始化的迭代数
    %   iter - 当前迭代次数
    %   params - 参数结构
    % 
    % 输出：
    %   should_reinit - 是否需要重初始化（逻辑值）
    %   reason - 触发原因（字符串，用于日志）
    %
    % 判据（按优先级）：
    %   1. 水平集变化过大
    %   2. 梯度模偏离1过多（符号距离性质退化）
    %   3. 超过最大允许间隔
    %   4. 前期/后期固定频率（保持原有逻辑）
    
    should_reinit = false;
    reason = '';
    
    h = params.grid.h;
    dx = params.grid.dx;
    dy = params.grid.dy;
    
    % 判据1：水平集变化过大
    lsf_change = max(abs(lsf(:) - lsf_before(:)));
    threshold_change = params.levelset.reinit_threshold * h;
    if lsf_change > threshold_change
        should_reinit = true;
        reason = sprintf('水平集变化过大 (%.3e > %.3e)', lsf_change, threshold_change);
        return;
    end
    
    % 判据2：梯度模偏离1过多（符号距离性质退化）
    [grad_y, grad_x] = gradient(lsf, dy, dx);
    grad_mag = hypot(grad_x, grad_y);
    grad_deviation = mean(abs(grad_mag(:) - 1));
    grad_tol = params.levelset.gradient_deviation_tol;
    if grad_deviation > grad_tol
        should_reinit = true;
        reason = sprintf('梯度偏差过大 (%.3f > %.3f)', grad_deviation, grad_tol);
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

