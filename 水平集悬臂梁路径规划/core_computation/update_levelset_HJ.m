function lsf_new = update_levelset_HJ(lsf, velocity, dt, dx, dy)
    % 采用迎风格式求解 Hamilton-Jacobi 方程更新水平集

    [ny, nx] = size(lsf);
    lsf_new = lsf;
    lsf_old = lsf_new;

    for i = 2:ny-1
        for j = 2:nx-1
            dpx = (lsf_old(i, j+1) - lsf_old(i, j)) / dx;
            dmx = (lsf_old(i, j) - lsf_old(i, j-1)) / dx;
            dpy = (lsf_old(i+1, j) - lsf_old(i, j)) / dy;
            dmy = (lsf_old(i, j) - lsf_old(i-1, j)) / dy;

            if velocity(i, j) > 0
                % 正速度：使用标准Godunov迎风格式（取最大值而非累加）
                % X方向：dmx为plus，dpx为minus；Y方向：dmy为plus，dpy为minus
                a_plus  = max(dmx, 0);
                a_minus = min(dpx, 0);
                b_plus  = max(dmy, 0);
                b_minus = min(dpy, 0);
                grad = sqrt( max(a_plus^2, (-a_minus)^2) + ...
                             max(b_plus^2, (-b_minus)^2) );
                lsf_new(i, j) = lsf_old(i, j) - dt * velocity(i, j) * grad;
            else
                % 负速度：使用标准Godunov迎风格式（取最大值而非累加）
                % X方向：dpx为plus，dmx为minus；Y方向：dpy为plus，dmy为minus
                a_plus  = max(dpx, 0);
                a_minus = min(dmx, 0);
                b_plus  = max(dpy, 0);
                b_minus = min(dmy, 0);
                grad = sqrt( max(a_plus^2, (-a_minus)^2) + ...
                             max(b_plus^2, (-b_minus)^2) );
                lsf_new(i, j) = lsf_old(i, j) - dt * velocity(i, j) * grad;
            end
        end
    end
end

