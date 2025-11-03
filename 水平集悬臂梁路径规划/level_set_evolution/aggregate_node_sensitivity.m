function node_sensitivity = aggregate_node_sensitivity(element_sensitivity, theta_field, lsf, nelx, nely, dx, dy, band_mask)
    % 通过链式法则将单元灵敏度汇总到水平集节点
    % 
    % 修改时间：2025-10-31
    % 改造内容：按照《链式求和_论文规范.txt》要求，使用精确形函数梯度公式
    %          和零水平集邻域过滤（仅对主路径穿过的单元求和）

    % 参数默认值
    if nargin < 8 || isempty(band_mask)
        h = min(dx, dy);
        band_mask = abs(lsf) <= 1.0 * h;  % 默认1.0h窄带
    end

    global DIAG;
    node_sensitivity = zeros(size(lsf));
    
    % 统计计数器
    elements_processed = 0;
    elements_skipped_band = 0;
    elements_skipped_grad = 0;

    for ely = 1:nely
        for elx = 1:nelx
            ci = ely + 1;
            cj = elx + 1;

            % === 零水平集邻域过滤（论文规范第2条）===
            % 仅对主路径穿过的单元求和：至少一个节点在窄带内
            in_band = false;
            for ni = [ci-1, ci, ci+1]
                for nj = [cj-1, cj, cj+1]
                    if band_mask(ni, nj)
                        in_band = true;
                        break;
                    end
                end
                if in_band, break; end
            end
            
            if ~in_band
                elements_skipped_band = elements_skipped_band + 1;
                continue;  % 跳过远离主路径的单元
            end

            % === 形函数梯度（四节点双线性单元，单元中心点） ===
            % 论文规范第3条：使用形函数求和
            % 节点编号（逆时针）：
            %   3 ------- 4
            %   |         |
            %   |    o    |  (中心点 ξ=0, η=0)
            %   |         |
            %   1 ------- 2
            %
            % 形函数及导数（在单元中心 ξ=0, η=0）：
            % N1 = (1-ξ)(1-η)/4,  ∂N1/∂ξ = -(1-η)/4,  ∂N1/∂η = -(1-ξ)/4
            % N2 = (1+ξ)(1-η)/4,  ∂N2/∂ξ =  (1-η)/4,  ∂N2/∂η = -(1+ξ)/4
            % N3 = (1-ξ)(1+η)/4,  ∂N3/∂ξ = -(1+η)/4,  ∂N3/∂η =  (1-ξ)/4
            % N4 = (1+ξ)(1+η)/4,  ∂N4/∂ξ =  (1+η)/4,  ∂N4/∂η =  (1+ξ)/4
            %
            % 在中心点 (ξ=0, η=0)：
            dN_dxi  = [-1, 1, -1, 1] / 4;
            dN_deta = [-1, -1, 1, 1] / 4;
            
            % Jacobi 矩阵（均匀网格简化版）
            % J = [dx/2, 0; 0, dy/2];
            % det_J = dx * dy / 4;
            % inv_J = [2/dx, 0; 0, 2/dy];
            
            % 物理坐标导数
            % [∂N/∂x; ∂N/∂y] = inv(J) * [∂N/∂ξ; ∂N/∂η]
            dN_dx = (2/dx) * dN_dxi;   % = 2/dx * dN_dxi
            dN_dy = (2/dy) * dN_deta;  % = 2/dy * dN_deta
            
            % 四个节点的 φ 值（按编号顺序）
            phi_nodes = [lsf(ci, cj-1), lsf(ci, cj+1), ...   % 节点1：左下, 节点2：右下
                         lsf(ci-1, cj), lsf(ci+1, cj)];      % 节点3：左上, 节点4：右上
            
            % 计算梯度（论文公式：∂φ/∂x = Σ_k (∂N_k/∂x) φ_k）
            phi_x = sum(dN_dx .* phi_nodes);
            phi_y = sum(dN_dy .* phi_nodes);

            if ~isempty(DIAG)
                DIAG.theta_grad_samples = DIAG.theta_grad_samples + 1;
            end

            % === 梯度过滤与正则化（保留原有稳定性机制）===
            grad_sq = phi_x^2 + phi_y^2;
            grad_threshold = 0.05;
            
            if grad_sq < grad_threshold^2
                % |∇φ| < 0.05，影响极小，直接跳过
                elements_skipped_grad = elements_skipped_grad + 1;
                if ~isempty(DIAG)
                    DIAG.theta_grad_small = DIAG.theta_grad_small + 1;
                    DIAG.den_small = DIAG.den_small + 1;
                end
                continue;
            end
            
            % ε²正则化（防止分母过小导致数值爆炸）
            h = min(dx, dy);
            epsilon = 0.1 * h;  % 论文推荐：0.1倍网格步长
            grad_sq_reg = grad_sq + epsilon^2;  % 正则化分母

            % === 计算 ∂θ_e/∂φ_i（论文规范第3条，包含 P_i, Q_i 项）===
            % θ_e = π/2 + atan2(∂φ/∂y, ∂φ/∂x)
            % ∂θ_e/∂φ_i = [∂(∂φ/∂y)/∂φ_i * (∂φ/∂x) - (∂φ/∂y) * ∂(∂φ/∂x)/∂φ_i] / [(∂φ/∂x)^2 + (∂φ/∂y)^2]
            %
            % 其中：
            %   ∂(∂φ/∂x)/∂φ_i = ∂N_i/∂x
            %   ∂(∂φ/∂y)/∂φ_i = ∂N_i/∂y
            
            dtheta_dphi = zeros(1, 4);
            for k = 1:4
                % 对第 k 个节点求导
                dtheta_dphi(k) = (dN_dy(k) * phi_x - phi_y * dN_dx(k)) / grad_sq_reg;
            end
            
            % === 链式求和（论文规范第2条）===
            coeff = element_sensitivity(ely, elx);
            contribs = coeff * dtheta_dphi;
            
            % 分配到四个节点（按编号映射回网格索引）
            node_positions = [ci, cj-1;    % 节点1：左下
                              ci, cj+1;    % 节点2：右下
                              ci-1, cj;    % 节点3：左上
                              ci+1, cj];   % 节点4：右上
            
            for k = 1:4
                ni = node_positions(k, 1);
                nj = node_positions(k, 2);
                node_sensitivity(ni, nj) = node_sensitivity(ni, nj) + contribs(k);
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
        fprintf('  [节点灵敏度聚合] 处理单元=%d, 跳过(窄带外)=%d, 跳过(梯度小)=%d\n', ...
            elements_processed, elements_skipped_band, elements_skipped_grad);
    end
end
