function consistency = check_sign_consistency(lsf)
    % 检查沿随机射线时水平集符号是否最多变化一次

    [rows, cols] = size(lsf);
    num_rays = 20;
    consistent_rays = 0;

    for ray = 1:num_rays
        start_i = randi([2, rows-1]);
        start_j = randi([2, cols-1]);
        angle = rand() * 2 * pi;
        dir_i = sin(angle);
        dir_j = cos(angle);

        max_steps = min(rows, cols) / 2;
        sign_changes = 0;
        prev_sign = sign(lsf(start_i, start_j));

        for step = 1:max_steps
            curr_i = round(start_i + step * dir_i);
            curr_j = round(start_j + step * dir_j);

            if curr_i >= 1 && curr_i <= rows && curr_j >= 1 && curr_j <= cols
                curr_sign = sign(lsf(curr_i, curr_j));
                if curr_sign ~= prev_sign && curr_sign ~= 0 && prev_sign ~= 0
                    sign_changes = sign_changes + 1;
                end
                prev_sign = curr_sign;
            else
                break;
            end
        end

        if sign_changes <= 1
            consistent_rays = consistent_rays + 1;
        end
    end

    consistency = consistent_rays / num_rays;
end

