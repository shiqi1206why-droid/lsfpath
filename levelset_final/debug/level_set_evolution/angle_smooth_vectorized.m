function [z_smooth, smooth_cache] = angle_smooth_vectorized(theta_e, eta, num_iters)
    % 二倍角向量化平滑（替代三重嵌套循环）
    % 
    % 输入：
    %   theta_e - 角度场 [nely, nelx]
    %   eta - 平滑系数 (0.05~0.15，推荐0.10)
    %   num_iters - 平滑迭代次数 (通常2次)
    % 
    % 输出：
    %   z_smooth - 平滑后的二倍角复向量 [nely, nelx]
    %
    % 性能优化：使用conv2矢量化替代三重嵌套循环，预期加速80%
    % 参考：fiber_levelset优化方案-融合版.md 阶段1.1
    
    % 转换为二倍角复向量（避免0/π环绕问题）
    z = cos(2*theta_e) + 1i*sin(2*theta_e);
    smooth_cache = struct();
    smooth_cache.theta_input = theta_e;
    smooth_cache.eta = eta;
    smooth_cache.num_iters = num_iters;
    smooth_cache.iterations = repmat(struct( ...
        'z_before_conv', [], ...
        'z_after_conv', [], ...
        'z_after_bc', [], ...
        'norm_before_clamp', [], ...
        'clamped_norm', [], ...
        'clamped_mask', [], ...
        'z_after_normalize', []), num_iters, 1);
    
    % 拉普拉斯核（五点模板）
    laplacian_kernel = [0, 1, 0; 
                        1, -4, 1; 
                        0, 1, 0];
    
    for k = 1:num_iters
        z_before_conv = z;
        % 矢量化卷积（替代嵌套循环）
        lap = conv2(z_before_conv, laplacian_kernel, 'same');
        z_after_conv = z_before_conv + eta * lap;
        
        % Neumann边界条件（保持原有逻辑）
        z_after_bc = apply_neumann_boundary(z_after_conv);
        
        % 归一化（防止幅度漂移）
        norm_before_clamp = abs(z_after_bc);
        clamped_norm = max(norm_before_clamp, 1e-12);
        clamped_mask = norm_before_clamp <= 1e-12;
        z = z_after_bc ./ clamped_norm;

        smooth_cache.iterations(k).z_before_conv = z_before_conv;
        smooth_cache.iterations(k).z_after_conv = z_after_conv;
        smooth_cache.iterations(k).z_after_bc = z_after_bc;
        smooth_cache.iterations(k).norm_before_clamp = norm_before_clamp;
        smooth_cache.iterations(k).clamped_norm = clamped_norm;
        smooth_cache.iterations(k).clamped_mask = clamped_mask;
        smooth_cache.iterations(k).z_after_normalize = z;
    end
    
    z_smooth = z;
    smooth_cache.z_output = z_smooth;
end
