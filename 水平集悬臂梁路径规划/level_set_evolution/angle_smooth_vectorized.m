function z_smooth = angle_smooth_vectorized(theta_e, eta, num_iters)
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
    
    % 拉普拉斯核（五点模板）
    laplacian_kernel = [0, 1, 0; 
                        1, -4, 1; 
                        0, 1, 0];
    
    for k = 1:num_iters
        % 矢量化卷积（替代嵌套循环）
        lap = conv2(z, laplacian_kernel, 'same');
        z = z + eta * lap;
        
        % Neumann边界条件（保持原有逻辑）
        z(1,:) = z(2,:);
        z(end,:) = z(end-1,:);
        z(:,1) = z(:,2);
        z(:,end) = z(:,end-1);
        
        % 归一化（防止幅度漂移）
        z = z ./ max(abs(z), 1e-12);
    end
    
    z_smooth = z;
end

