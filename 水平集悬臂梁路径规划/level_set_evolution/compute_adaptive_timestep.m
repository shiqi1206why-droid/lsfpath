function dt_adaptive = compute_adaptive_timestep(velocity_field, dx, dy, cfl_factor)
    % 根据速度稳健幅值按 CFL 原则计算自适应时间步长

    if nargin < 4 || isempty(cfl_factor)
        cfl_factor = 0.08;  % 强力稳住4：从0.1降至0.08
    end

    v_abs = abs(velocity_field(:));
    v_abs = v_abs(isfinite(v_abs));
    if isempty(v_abs)
        dt_adaptive = 0.1;
        return;
    end

    v_sorted = sort(v_abs);
    k = max(1, round(0.99 * numel(v_sorted)));
    v_robust = v_sorted(k);

    if v_robust < 1e-12
        dt_adaptive = 0.1;
        return;
    end

    grid_spacing = min(dx, dy);
    dt_adaptive = cfl_factor * grid_spacing / v_robust;
    dt_adaptive = max(min(dt_adaptive, 0.5), 1e-6);
end

