function node_sensitivity = aggregate_node_sensitivity(element_sensitivity, theta_field, lsf, nelx, nely, dx, dy)
    % 通过链式法则将单元灵敏度汇总到水平集节点

    global DIAG;
    node_sensitivity = zeros(size(lsf));

    for ely = 1:nely
        for elx = 1:nelx
            ci = ely + 1;
            cj = elx + 1;

            phi_x = (lsf(ci, cj+1) - lsf(ci, cj-1)) / (2*dx);
            phi_y = (lsf(ci+1, cj) - lsf(ci-1, cj)) / (2*dy);

            if ~isempty(DIAG)
                DIAG.theta_grad_samples = DIAG.theta_grad_samples + 1;
            end

            % === 修改14.1：链式法则梯度正则化 ===
            % 参考：代码修改交流.md 2025-10-20，论文推荐方法
            
            % 第1步：梯度阈值过滤（跳过影响极小的单元）
            h = min(dx, dy);
            grad_sq = phi_x^2 + phi_y^2;
            grad_threshold = 0.05;  % 论文推荐值
            
            if grad_sq < grad_threshold^2
                % |∇φ| < 0.05，影响极小，直接跳过
                if ~isempty(DIAG)
                    DIAG.theta_grad_small = DIAG.theta_grad_small + 1;
                    DIAG.den_small = DIAG.den_small + 1;
                end
                continue;
            end
            
            % 第2步：ε²正则化（防止分母过小导致数值爆炸）
            epsilon = 0.1 * h;  % 论文推荐：0.1倍网格步长
            grad_sq_reg = grad_sq + epsilon^2;  % 正则化分母

            coeff = element_sensitivity(ely, elx);
            % 使用正则化后的分母，防止灵敏度爆炸至10⁵级
            dtheta_dphi_y_plus = (phi_x / grad_sq_reg) * (1/(2*dy));
            dtheta_dphi_y_minus = (phi_x / grad_sq_reg) * (-1/(2*dy));
            dtheta_dphi_x_plus = (-phi_y / grad_sq_reg) * (1/(2*dx));
            dtheta_dphi_x_minus = (-phi_y / grad_sq_reg) * (-1/(2*dx));

            contribs = coeff * [dtheta_dphi_y_plus, dtheta_dphi_y_minus, dtheta_dphi_x_plus, dtheta_dphi_x_minus];

            node_sensitivity(ci+1, cj) = node_sensitivity(ci+1, cj) + contribs(1);
            node_sensitivity(ci-1, cj) = node_sensitivity(ci-1, cj) + contribs(2);
            node_sensitivity(ci, cj+1) = node_sensitivity(ci, cj+1) + contribs(3);
            node_sensitivity(ci, cj-1) = node_sensitivity(ci, cj-1) + contribs(4);

            if ~isempty(DIAG)
                DIAG.contrib_total = DIAG.contrib_total + 4;
                DIAG.contrib_nonzero = DIAG.contrib_nonzero + nnz(contribs);
            end
        end
    end
end

