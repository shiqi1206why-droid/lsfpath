function actual_spacing = verify_path_spacing(lsf, level1, level2, dx, dy)
    % 沿法向采样，估计两条水平集等值线的间距

    [nely_plus, nelx_plus] = size(lsf);
    spacings = [];
    sample_points = 10;

    for sample = 1:sample_points
        i = randi([2, nely_plus-1]);
        j = randi([2, nelx_plus-1]);

        if abs(lsf(i, j) - level1) < 0.1
            [grad_y, grad_x] = gradient(lsf, dy, dx);
            grad_norm = sqrt(grad_x(i, j)^2 + grad_y(i, j)^2);

            if grad_norm > 1e-6
                gx = grad_x(i, j) / grad_norm;
                gy = grad_y(i, j) / grad_norm;
                max_search = 20;
                for step = 1:max_search
                    ni = round(i + step * gy);
                    nj = round(j + step * gx);

                    if ni >= 1 && ni <= nely_plus && nj >= 1 && nj <= nelx_plus
                        if abs(lsf(ni, nj) - level2) < 0.1
                            physical_dist = step * sqrt((dx * gx)^2 + (dy * gy)^2);
                            spacings = [spacings, physical_dist]; %#ok<*AGROW>
                            break;
                        end
                    else
                        break;
                    end
                end
            end
        end
    end

    if isempty(spacings)
        actual_spacing = abs(level2 - level1);
    else
        actual_spacing = mean(spacings);
    end
end

