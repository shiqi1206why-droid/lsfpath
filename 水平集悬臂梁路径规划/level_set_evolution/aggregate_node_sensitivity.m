function node_sensitivity = aggregate_node_sensitivity(element_sensitivity, lsf_nodes, nelx, nely, dx, dy, band_mask_nodes)
    % 通过链式法则将单元灵敏度汇总到水平集节点
    % 
    % 修改时间：2025-10-31
    % 改造内容：严格遵循《链式求和严格实现方案.md》
    %          - 使用角点水平集值（节点网格）
    %          - 严格的零水平集穿越判定
    %          - 精确的形函数梯度公式
    %          - 无近似和宽松筛选
    %
    % 输入：
    %   element_sensitivity - (nely x nelx) 单元灵敏度
    %   lsf_nodes          - (nely+1 x nelx+1) 节点水平集值
    %   nelx, nely         - 网格尺寸
    %   dx, dy             - 网格步长
    %   band_mask_nodes    - (nely+1 x nelx+1) 窄带掩膜
    %
    % 输出：
    %   node_sensitivity   - (nely+1 x nelx+1) 节点灵敏度

    global DIAG;
    
    % 初始化节点灵敏度（与 lsf_nodes 同尺寸）
    node_sensitivity = zeros(size(lsf_nodes));
    
    % 统计计数器
    elements_processed = 0;
    elements_skipped_no_crossing = 0;
    elements_skipped_grad = 0;
    
    % 网格步长
    h = min(dx, dy);
    
    % 参数设置（严格模式）
    phi_zero_tol = 0.01 * h;      % 近零容差（捕捉等值线贴在角点的情况）
    grad_eps = 0.1 * h;          % 梯度保底值（防止分母为零）

    for ely = 1:nely
        for elx = 1:nelx
            % === 严格的交叉判定 ===
            % 四个角点的 φ 值（Q4单元标准索引）
            phi11 = lsf_nodes(ely,   elx  );  % 左下
            phi12 = lsf_nodes(ely,   elx+1);  % 右下
            phi21 = lsf_nodes(ely+1, elx  );  % 左上
            phi22 = lsf_nodes(ely+1, elx+1);  % 右上
            
            % 判据1：边穿越（任一边符号异号）
            edge_crossed = (phi11 * phi12 < 0) || (phi21 * phi22 < 0) || ...
                           (phi11 * phi21 < 0) || (phi12 * phi22 < 0);
            
            % 判据2：近零容差（等值线贴在角点）
            near_zero = min(abs([phi11, phi12, phi21, phi22])) < phi_zero_tol;
            
            % 判据3：至少一个角点在窄带内
            in_band = any(any(band_mask_nodes(ely:ely+1, elx:elx+1)));
            
            % 组合判据：必须在窄带内，且（穿越 或 近零）
            if ~(in_band && (edge_crossed || near_zero))
                elements_skipped_no_crossing = elements_skipped_no_crossing + 1;
                continue;  % 不满足穿越条件，跳过
            end
            
            % === 精确的梯度与导数 ===
            % 形函数导数（Q4单元，标准单元中心 ξ=0, η=0）
            dN_dxi  = [-1, 1, -1, 1] / 4;   % 对应节点顺序：11, 12, 21, 22
            dN_deta = [-1, -1, 1, 1] / 4;
            
            % Jacobi 逆矩阵（均匀网格）
            invJ = [2/dx, 0; 0, 2/dy];
            
            % 物理坐标导数
            dN_dx_row = invJ(1,1) * dN_dxi;   % = 2/dx * dN_dxi
            dN_dy_row = invJ(2,2) * dN_deta;  % = 2/dy * dN_deta
            
            % 四个角点的 φ 值（按形函数节点顺序）
            phi_nodes = [phi11, phi12, phi21, phi22];
            
            % 计算梯度
            phi_x = sum(dN_dx_row .* phi_nodes);
            phi_y = sum(dN_dy_row .* phi_nodes);
            
            if ~isempty(DIAG)
                DIAG.theta_grad_samples = DIAG.theta_grad_samples + 1;
            end
            
            % 梯度模平方
            grad_sq = phi_x^2 + phi_y^2;
            
            % 保底检查：防止分母为零（使用极小值，不引入经验阈值）
            if grad_sq < grad_eps^2
                elements_skipped_grad = elements_skipped_grad + 1;
                if ~isempty(DIAG)
                    DIAG.theta_grad_small = DIAG.theta_grad_small + 1;
                    DIAG.den_small = DIAG.den_small + 1;
                end
                continue;
            end
            
            % === 计算 ∂θ_e/∂φ_i（精确公式，无正则化） ===
            % θ_e = π/2 + atan2(∂φ/∂y, ∂φ/∂x)
            % ∂θ_e/∂φ_i = [∂N_i/∂y * (∂φ/∂x) - (∂φ/∂y) * ∂N_i/∂x] / [(∂φ/∂x)^2 + (∂φ/∂y)^2]
            
            dtheta_dphi = zeros(1, 4);
            for k = 1:4
                dtheta_dphi(k) = (dN_dy_row(k) * phi_x - phi_y * dN_dx_row(k)) / grad_sq;
            end
            
            % === 链式求和 ===
            coeff = element_sensitivity(ely, elx);
            contribs = coeff * dtheta_dphi;
            
            % 鲁棒性检查：只累加有效值
            if isfinite(coeff)
                % 累加到四个角点（直接使用节点索引）
                if isfinite(contribs(1))
                    node_sensitivity(ely,   elx  ) = node_sensitivity(ely,   elx  ) + contribs(1);
                end
                if isfinite(contribs(2))
                    node_sensitivity(ely,   elx+1) = node_sensitivity(ely,   elx+1) + contribs(2);
                end
                if isfinite(contribs(3))
                    node_sensitivity(ely+1, elx  ) = node_sensitivity(ely+1, elx  ) + contribs(3);
                end
                if isfinite(contribs(4))
                    node_sensitivity(ely+1, elx+1) = node_sensitivity(ely+1, elx+1) + contribs(4);
                end
            end
            
            elements_processed = elements_processed + 1;
            
            if ~isempty(DIAG)
                DIAG.contrib_total = DIAG.contrib_total + 4;
                DIAG.contrib_nonzero = DIAG.contrib_nonzero + nnz(contribs);
            end
        end
    end
    
    % 输出统计信息（仅首次调用）
    persistent first_call;
    if isempty(first_call)
        first_call = false;
        fprintf('  [节点灵敏度聚合-严格模式] 处理=%d, 跳过(未穿越)=%d, 跳过(梯度极小)=%d\n', ...
            elements_processed, elements_skipped_no_crossing, elements_skipped_grad);
        fprintf('  [严格参数] 近零容差=%.2e·h, 梯度保底=%.2e·h\n', ...
            phi_zero_tol/h, grad_eps/h);
    end
end
