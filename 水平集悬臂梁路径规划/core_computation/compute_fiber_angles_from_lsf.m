function theta_e = compute_fiber_angles_from_lsf(lsf, dx, dy)
    % 根据水平集梯度计算局部纤维方向

    global DIAG;
    [nely_plus, nelx_plus] = size(lsf);
    nely = nely_plus - 2;
    nelx = nelx_plus - 2;
    theta_e = zeros(nely, nelx);

    for i = 2:nely_plus-1
        for j = 2:nelx_plus-1
            dphi_dx = (lsf(i, j+1) - lsf(i, j-1)) / (2*dx);
            dphi_dy = (lsf(i+1, j) - lsf(i-1, j)) / (2*dy);

            if ~isempty(DIAG)
                DIAG.theta_grad_samples = DIAG.theta_grad_samples + 1;
                if (dphi_dx*dphi_dx + dphi_dy*dphi_dy) < 1e-12
                    DIAG.theta_grad_small = DIAG.theta_grad_small + 1;
                end
            end

            % === 角度计算容错：防止NaN传播 ===
            % 参考：代码修改交流.md (Codex方案 2025-10-17)
            % 检测非有限值或过小梯度，赋默认角度
            if ~isfinite(dphi_dx) || ~isfinite(dphi_dy) || ...
               (dphi_dx*dphi_dx + dphi_dy*dphi_dy) < 1e-12
                theta_e(i-1, j-1) = 0;  % 默认角度（可选：继承上一次迭代）
            else
                theta_e(i-1, j-1) = pi/2 + atan2(dphi_dy, dphi_dx);
            end

            theta_e(i-1, j-1) = mod(theta_e(i-1, j-1), pi);
        end
    end
end

