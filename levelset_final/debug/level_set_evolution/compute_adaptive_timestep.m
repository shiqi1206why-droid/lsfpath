function dt_adaptive = compute_adaptive_timestep(velocity_field, dx, dy, cfl_factor)
    % 根据速度稳健幅值按 CFL 原则计算自适应时间步长
    % velocity_field: 水平集窄带上的速度
    % dx, dy: 网格尺寸
    % cfl_factor: CFL 系数（可选，默认 0.08）

    % 如果未显式传入系数，则使用经验值 0.08（比 0.1 更保守）
    if nargin < 4 || isempty(cfl_factor)
        cfl_factor = 0.08;
    end

    % 提取速度绝对值并剔除 NaN/Inf，保证统计稳定
    v_abs = abs(velocity_field(:));
    v_abs = v_abs(isfinite(v_abs));
    if isempty(v_abs)
        % 没有有效速度时给一个温和的步长（避免除零）
        dt_adaptive = 0.1;
        return;
    end

    % 按绝对值排序，取 95% 分位数作为“稳健幅值”，对抗极端 outlier
    v_sorted = sort(v_abs);
    k = max(1, round(0.95 * numel(v_sorted)));
    v_robust = v_sorted(k);

    if v_robust < 1e-12
        % 速度几乎全 0 时直接返回默认步长
        dt_adaptive = 0.1;
        return;
    end

    % CFL 时间步：Δt = c * h / v，其中 h 为最小网格间距
    grid_spacing = min(dx, dy);
    dt_adaptive = cfl_factor * grid_spacing / v_robust;

    % 限幅：下限 0.01h -> 至少移动 1% 网格；上限 0.5 -> 避免过大
    dt_min = 0.01 * grid_spacing;
    dt_adaptive = max(min(dt_adaptive, 0.5), dt_min);
end
