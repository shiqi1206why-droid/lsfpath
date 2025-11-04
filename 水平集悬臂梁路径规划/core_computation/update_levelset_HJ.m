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
                grad_plus = sqrt(max(dmx, 0)^2 + min(dpx, 0)^2 + ...
                    max(dmy, 0)^2 + min(dpy, 0)^2);
                lsf_new(i, j) = lsf_old(i, j) - dt * velocity(i, j) * grad_plus;
            else
                grad_minus = sqrt(min(dmx, 0)^2 + max(dpx, 0)^2 + ...
                    min(dmy, 0)^2 + max(dpy, 0)^2);
                lsf_new(i, j) = lsf_old(i, j) - dt * velocity(i, j) * grad_minus;
            end
        end
    end
end

