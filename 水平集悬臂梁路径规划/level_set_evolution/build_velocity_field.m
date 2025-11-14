function [velocity_field, stats] = build_velocity_field(node_sensitivity, lsf, dx, dy, bandwidth, remove_bias, gamma_curv, iter)
    % 构造用于水平集演化的窄带法向速度场

    if nargin < 5 || isempty(bandwidth)
        bandwidth = 2 * min(dx, dy);
    end
    if nargin < 6 || isempty(remove_bias)
        remove_bias = true;
    end
    if nargin < 7 || isempty(gamma_curv)
        gamma_curv = 0;  % 默认不启用曲率正则化
    end
    if nargin < 8 || isempty(iter)
        iter = 1;  % 默认迭代1
    end

    [ny_lsf, nx_lsf] = size(lsf);
    [ny_sens, nx_sens] = size(node_sensitivity);

    if ny_lsf ~= ny_sens || nx_lsf ~= nx_sens
        sensitivity_expanded = zeros(ny_lsf, nx_lsf);
        sensitivity_expanded(2:ny_sens+1, 2:nx_sens+1) = node_sensitivity;
        sensitivity_expanded(1, :) = sensitivity_expanded(2, :);
        sensitivity_expanded(end, :) = sensitivity_expanded(end-1, :);
        sensitivity_expanded(:, 1) = sensitivity_expanded(:, 2);
        sensitivity_expanded(:, end) = sensitivity_expanded(:, end-1);
        node_sensitivity = sensitivity_expanded;
    end


    % 步骤：1) v_shape去均值 → 2) 合成 → 3) 轻微抑制净平移
    
    stats = struct('bandwidth', bandwidth, 'band_fraction', 0, 'mean_band', 0, ...
                   'weighted_mean', 0, 'max_band', 0, 'pos_ratio', 0, ...
                   'neg_ratio', 0, 'length', 0, 'centroid', [NaN, NaN]);

    % 步骤1：v_shape去均值（仅对形状项）
    v_shape = -node_sensitivity;
    
    if remove_bias
        band_mask_temp = abs(lsf) <= bandwidth;
        if any(band_mask_temp(:))
            [grad_y, grad_x] = gradient(lsf, dy, dx);
            grad_magnitude = hypot(grad_x, grad_y);
            grad_magnitude = max(grad_magnitude, 1e-12);
            weight_field = grad_magnitude;
            weight_field(~band_mask_temp) = 0;
            weight_values_temp = weight_field(band_mask_temp);
            v_shape_band = v_shape(band_mask_temp);
            
            if any(weight_values_temp)
                weighted_mean_shape = sum(v_shape_band .* weight_values_temp) / sum(weight_values_temp);
                v_shape(band_mask_temp) = v_shape_band - weighted_mean_shape;
            else
                weighted_mean_shape = 0;
            end
        else
            weighted_mean_shape = 0;
        end
    else
        weighted_mean_shape = 0;
    end
    
    % 初始化velocity_field为去均值后的v_shape
    velocity_field = v_shape;

    % === 重新计算掩膜和梯度（用于后续操作） ===
    [grad_y, grad_x] = gradient(lsf, dy, dx);
    grad_magnitude = hypot(grad_x, grad_y);
    grad_magnitude = max(grad_magnitude, 1e-12);

    % 仅在带宽内更新，其余位置速度置零
    band_mask = abs(lsf) <= bandwidth;
    if ~any(band_mask(:))
        velocity_field(:) = 0;
        return;
    end
    velocity_field(~band_mask) = 0;

    % === 幅值缩放（Bug 4修复：仅缩放v_shape，不缩放约束项） ===
    % 在加入v_fid之前进行缩放，确保约束力不被削弱
    % 注意：此时velocity_field = v_shape（去均值后）
    weight_field_scale = grad_magnitude;
    weight_field_scale(~band_mask) = 0;
    weight_values_scale = weight_field_scale(band_mask);
    band_values = velocity_field(band_mask);
    
    v_abs = abs(band_values);
    if any(~isfinite(v_abs))
        warning('velocity_field: v_abs 包含非有限值，已强制置零');
        velocity_field(~isfinite(velocity_field)) = 0;
    end

    % 等距约束v_fid已禁用，依赖FMM保持等距信息

    % === 曲率正则化（纤维路径-修改思路.txt 三-5)、三-6)） ===
    % 平滑路径，抑制锯齿状震荡和局部尖角
    if gamma_curv > 0
        % 计算归一化梯度
        grad_mag_curv = hypot(grad_x, grad_y);
        grad_mag_curv = max(grad_mag_curv, 1e-12);  % 防止除零
        
        nx = grad_x ./ grad_mag_curv;
        ny = grad_y ./ grad_mag_curv;
        
        % 计算曲率：κ = div(∇φ/|∇φ|)
        [nx_x, ~] = gradient(nx, dx, dy);
        [~, ny_y] = gradient(ny, dx, dy);
        kappa = nx_x + ny_y;
        
        % 曲率正则化项
        v_curv = -gamma_curv * kappa;
        v_curv(~band_mask) = 0;
        
        % 添加到速度场
        velocity_field = velocity_field + v_curv;
        
        % 统计曲率信息
        kappa_band = kappa(band_mask);
        if ~isempty(kappa_band)
            fprintf('  [曲率正则] gamma=%.4f, 曲率范围=[%.3f, %.3f], 平均=%.3f\n', ...
                gamma_curv, min(kappa_band), max(kappa_band), mean(kappa_band));
        end
    end

    % === 修改20-A3：暂停净平移抑制（让v_fid充分发挥） ===
    % 参考：解决方案-分点总结.txt A3
    % beta=0.2可能削弱了v_fid锚定效果，暂停去偏（稳定后可恢复0.1）
    beta = 0.0;  % 暂停去偏
    if any(band_mask(:))
        weight_field_final = grad_magnitude;
        weight_field_final(~band_mask) = 0;
        weight_values_final = weight_field_final(band_mask);
        velocity_band_final = velocity_field(band_mask);
        
        if sum(weight_values_final) > 0
            mu_total = sum(velocity_band_final .* weight_values_final) / sum(weight_values_final);
            velocity_field(band_mask) = velocity_band_final - beta * mu_total;
        end
    end
    
    % === 更新band_values（用于统计） ===
    band_values = velocity_field(band_mask);

    stats.band_fraction = nnz(band_mask) / numel(band_mask);
    stats.mean_band = mean(band_values);
    stats.max_band = max(abs(band_values));
    stats.pos_ratio = mean(band_values > 0);
    stats.neg_ratio = mean(band_values < 0);

    % 使用最终的weight_field计算统计量
    weight_field_stats = grad_magnitude;
    weight_field_stats(~band_mask) = 0;
    weight_values_stats = weight_field_stats(band_mask);
    weights_sum = sum(weight_values_stats);
    
    if weights_sum > 0
        [ny, nx] = size(lsf);
        x_coords = linspace(0, dx*(nx-1), nx);
        y_coords = linspace(0, dy*(ny-1), ny);
        [X_grid, Y_grid] = meshgrid(x_coords, y_coords);
        stats.length = weights_sum * min(dx, dy);
        stats.centroid = [sum(X_grid(band_mask).*weight_values_stats) / weights_sum, ...
                          sum(Y_grid(band_mask).*weight_values_stats) / weights_sum];
        stats.weighted_mean = sum(band_values .* weight_values_stats) / weights_sum;
    else
        stats.length = 0;
        stats.centroid = [NaN, NaN];
        stats.weighted_mean = 0;
    end
end

