function [velocity_field, stats] = build_velocity_field(node_sensitivity, lsf, dx, dy, bandwidth, remove_bias, lsf_target, lambda_fid, material_mask, gamma_curv, iter)
    % 构造用于水平集演化的窄带法向速度场
    % 参数：lsf_target为边界等距目标场，lambda_fid为等距约束权重
    %       material_mask用于边界邻域限制，gamma_curv为曲率正则化系数
    %       iter为当前迭代次数（用于动态调整v_target等参数）

    if nargin < 5 || isempty(bandwidth)
        bandwidth = 2 * min(dx, dy);
    end
    if nargin < 6 || isempty(remove_bias)
        remove_bias = true;
    end
    if nargin < 7 || isempty(lsf_target)
        lsf_target = lsf;  % 若未提供，使用当前lsf（退化为无约束）
    end
    if nargin < 8 || isempty(lambda_fid)
        lambda_fid = 0;  % 默认不启用等距约束
    end
    if nargin < 9
        material_mask = [];  % 默认不启用边界邻域限制
    end
    if nargin < 10 || isempty(gamma_curv)
        gamma_curv = 0;  % 默认不启用曲率正则化
    end
    if nargin < 11 || isempty(iter)
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

    % === 强力稳住7（细节A）：去均值顺序优化 ===
    % 参考：强力稳住-总结清单.txt 二-A)
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

    band_mask = abs(lsf) <= bandwidth;
    if ~any(band_mask(:))
        return;
    end

    % === 强力稳住3：运动域再收紧（环带+目标护栏）===
    % 参考：强力稳住-总结清单.txt 一-3)
    % 双重限制：材料环带 + 目标护栏
    if ~isempty(material_mask)
        % 提取材料边界外周
        boundary_perim = bwperim(material_mask);
        
        % 形态学膨胀（环带更紧）
        r = 1;  % 从2缩减到1，更紧
        boundary_region = imdilate(boundary_perim, strel('disk', r));
        
        % 扩展到含ghost cells的尺寸
        if size(boundary_region, 1) == size(lsf, 1) - 2
            boundary_region_expanded = false(size(lsf));
            boundary_region_expanded(2:end-1, 2:end-1) = boundary_region;
            boundary_region_expanded(1, :) = boundary_region_expanded(2, :);
            boundary_region_expanded(end, :) = boundary_region_expanded(end-1, :);
            boundary_region_expanded(:, 1) = boundary_region_expanded(:, 2);
            boundary_region_expanded(:, end) = boundary_region_expanded(:, end-1);
            boundary_region = boundary_region_expanded;
        end
        
        % 叠加材料环带
        band_mask = band_mask & boundary_region;
        
        % 增加目标护栏：限制在目标层附近3h范围
        h = min(dx, dy);
        target_ring = abs(lsf_target) <= 3*h;
        band_mask = band_mask & target_ring;
        
        % 统计运动域限制效果
        boundary_nodes = nnz(band_mask);
        fprintf('  [运动域收紧] 环带r=%d，目标护栏=±%.2fh，窄带节点数=%d\n', r, 3.0, boundary_nodes);
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
    if ~isempty(v_abs)
        v_sorted = sort(v_abs(:));
        k = max(1, round(0.99 * numel(v_sorted)));
        v_robust = v_sorted(k);
        if isfinite(v_robust) && v_robust > 0
            % === 强力稳住4：步幅更小（前期慢稳）===
            % 参考：强力稳住-总结清单.txt 一-4)
            % 前100步极慢，后期正常
            if iter <= 100
                v_target = 3;  % 前100步：极慢步幅
            else
                v_target = 5;  % 后期：正常步幅
            end
            scale = min(1, v_target / v_robust);
            if scale < 1
                velocity_field(band_mask) = band_values * scale;
                band_values = band_values * scale;
            end
        end
    end

    % === 强力稳住5：去掉v_fid限幅 ===
    % 参考：强力稳住-总结清单.txt 一-5)
    % V_fid = -λ_fid · (φ - φ_target)，将路径拉回到边界等距位置
    if lambda_fid > 0
        deviation = lsf - lsf_target;  % 偏离量
        % 大幅放宽限幅：让约束项充分发挥作用
        h = min(dx, dy);
        max_deviation = 100 * h;  % 从20*h大幅放宽至100*h，几乎不限幅
        deviation = sign(deviation) .* min(abs(deviation), max_deviation);
        v_fid = -lambda_fid * deviation;  % Bug 1修复：加负号
        v_fid(~band_mask) = 0;  % 仅在窄带内生效
        velocity_field = velocity_field + v_fid;  % Bug 1+4修复：在缩放后加入
    end

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

