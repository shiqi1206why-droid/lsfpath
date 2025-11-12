function results = fiber_levelset(config_name)
    % 纤维路径优化主函数（阶段2+3优化版）
    % 
    % 输入：
    %   config_name - 配置名称（可选，默认'default'）
    %                支持: 'default', 'fast', 'precise', 'debug'
    % 
    % 输出：
    %   results - 优化结构体
    
    % 处理输入参数
    if nargin < 1
        config_name = 'default';
    end
    
    % 清理（保留输入参数）
    clc; close all;
    clearvars -except config_name;
    
    % 自动添加所有辅助函数子文件夹到路径
    script_dir = fileparts(mfilename('fullpath'));
    addpath(genpath(script_dir));
    cleanup = onCleanup(@() rmpath(genpath(script_dir)));
    
    %% 1. 加载和验证参数配置
    fprintf('=== 纤维路径优化（配置: %s） ===\n', config_name);
    params = get_fiber_optimization_params(config_name);
    validate_params(params);
    
    % 解包常用参数（简化代码）
    nelx = params.grid.nelx;
    nely = params.grid.nely;
    Lx = params.grid.Lx;
    Ly = params.grid.Ly;
    dx = params.grid.dx;
    dy = params.grid.dy;
    
    E_L = params.material.E_L;
    E_T = params.material.E_T;
    nu_LT = params.material.nu_LT;
    nu_TL = params.material.nu_TL;
    G_LT = params.material.G_LT;
    G_LW = params.material.G_LW;
    G_TW = params.material.G_TW;
    thickness = params.material.thickness;
    
    max_iter = params.opt.max_iter;
    tol = params.opt.tol;
    alpha = params.opt.alpha;
    dt = params.opt.dt;
    delta_theta_max = params.opt.delta_theta_max;
    fidelity_weight = params.opt.fidelity_weight;

% 初始化参数
    delta_phi = params.levelset.delta_phi_factor * params.grid.h;
    init_smooth_opts = struct('morph_radius', params.init.morph_radius);

% 载荷参数
    F_mag = params.load.F_mag;

global DIAG; %#ok<GVMIS>
DIAG = struct();
diag_reset();

% 载入拓扑优化结果
topo_file = 'topo_result.mat';
if ~exist(topo_file, 'file')
    error('未找到拓扑优化结果文件: %s', topo_file);
end

log_message('INFO', params, '正在加载拓扑优化结果...');
topo_data = load(topo_file);

if ~isfield(topo_data, 'struc')
    error('拓扑结果文件缺少 struc 字段。');
end

struc = topo_data.struc;
log_message('INFO', params, '拓扑网格尺寸: %dx%d', size(struc, 2), size(struc, 1));

if isfield(topo_data, 'nelx') && isfield(topo_data, 'nely')
    if topo_data.nelx ~= nelx || topo_data.nely ~= nely
        log_message('WARN', params, '网格尺寸不一致：拓扑(%dx%d) vs 当前(%dx%d)，正在重新采样...', ...
            topo_data.nelx, topo_data.nely, nelx, nely);
        struc = imresize(struc, [nely, nelx], 'nearest');
        log_message('INFO', params, '重采样后拓扑网格尺寸: %dx%d', size(struc, 2), size(struc, 1));
    end
end

% 初始化与诊断
log_message('INFO', params, '开始执行边界偏移初始化...');
figure('Name', '拓扑与初始化检查', 'Position', [100, 100, 1200, 400]);

subplot(1,3,1);
imagesc(struc);
colormap(gray);
axis equal; axis tight;
title('原始拓扑');
xlabel('x方向单元索引');
ylabel('y方向单元索引');
log_message('INFO', params, '正在清理材料掩膜...');
[material_mask, mask_info] = clean_material_mask(struc, params.init.min_component_size, init_smooth_opts.morph_radius);

subplot(1,3,2);
imagesc(material_mask);
colormap(gray);
axis equal; axis tight;
title('清理后的材料掩膜');
xlabel('x方向单元索引');
ylabel('y方向单元索引');
hold on;
boundary_mask = bwperim(material_mask);
[boundary_y, boundary_x] = find(boundary_mask);
plot(boundary_x, boundary_y, 'r.', 'MarkerSize', 2, 'DisplayName', '材料边界');
if ~isempty(boundary_x)
    legend('材料边界', 'Location', 'best');
end
log_message('DEBUG', params, '  连通区域数: %d，总像素数: %d', mask_info.num_components, mask_info.total_area);

log_message('INFO', params, '正在构建基于边界偏移的符号距离场...');
[lsf, parallel_paths, init_info] = construct_boundary_offset_levelset_with_parallel( ...
    material_mask, nelx, nely, dx, dy, delta_phi, init_smooth_opts);

initial_zero_mask = init_info.zero_mask;
lsf_initial = lsf;

% === 保存边界等距目标场（用于速度场约束） ===
% 参考：论文关键要点.md 第4节、代码修改交流.md 2025-10-17方案
log_message('INFO', params, '正在构建边界等距目标场...');
phi_boundary_global = zeros(nely+2, nelx+2);

% === 优化1.3：预计算网格特征尺寸（全局使用，消除重复计算） ===
h_grid = min(dx, dy);  % 只计算一次，后续全部使用h_grid
h = h_grid;  % 保持向后兼容
d_in = bwdist(~material_mask) * h;
d_out = bwdist(material_mask) * h;
phi_boundary_core = d_out - d_in;
phi_boundary_global(2:end-1, 2:end-1) = phi_boundary_core;
% 边界外推（Neumann边界条件）
phi_boundary_global(1, :) = phi_boundary_global(2, :);
phi_boundary_global(end, :) = phi_boundary_global(end-1, :);
phi_boundary_global(:, 1) = phi_boundary_global(:, 2);
phi_boundary_global(:, end) = phi_boundary_global(:, end-1);
% 目标场：φ_target = φ_b + Δφ_used（与初始化严格一致）
lsf_target_global = phi_boundary_global + init_info.delta_phi_used;
log_message('INFO', params, '  边界等距目标场已构建（用于优化约束，Δφ_used=%.4f）', init_info.delta_phi_used);

% === 调试信息：检查init_info.contour内容 ===
if isfield(init_info, 'contour') && ~isempty(init_info.contour.x)
    fprintf('  [调试] init_info.contour采样点数量: %d\n', length(init_info.contour.x));
    fprintf('  [调试] X范围: [%.3f, %.3f] m\n', min(init_info.contour.x), max(init_info.contour.x));
    fprintf('  [调试] Y范围: [%.3f, %.3f] m\n', min(init_info.contour.y), max(init_info.contour.y));
else
    fprintf('  [警告] init_info.contour为空或不存在！\n');
end

subplot(1,3,3);
% 使用索引坐标系统，与图一图二保持一致（Codex方案）
[X_idx, Y_idx] = meshgrid(0:nelx+1, 0:nely+1);
contour(X_idx, Y_idx, lsf, 20, 'LineWidth', 0.5, 'DisplayName', '等值线');  % 灰色背景等值线
hold on;

% 直接从当前lsf重新提取零等值线，确保完整性
C_zero = contourc(0:nelx+1, 0:nely+1, lsf, [0 0]);
[zero_x, zero_y] = contourc_to_points(C_zero);
if ~isempty(zero_x)
    fprintf('  [调试] 直接提取的零等值线采样点数量: %d\n', length(zero_x));
    % 先绘制连线
    plot(zero_x, zero_y, 'r-', 'LineWidth', 2, 'DisplayName', '主路径 \phi=0');
    % 再绘制采样点
    plot(zero_x, zero_y, 'mo', 'MarkerSize', 3, 'MarkerFaceColor', 'm', ...
        'LineStyle', 'none', 'DisplayName', '主路径采样点');
else
    fprintf('  [警告] 直接提取的零等值线也为空！\n');
end

axis equal; axis tight;
set(gca,'YDir','reverse');          % Y轴向下，与imagesc一致
title('初始化水平集等值线');
xlabel('x方向单元索引');
ylabel('y方向单元索引');
legend('Location','best');
colorbar;


fprintf('  目标Δφ = %.4f m，抽样均值 = %.4f m (标准差 = %.4f m，样本数 = %d)\n', ...
    delta_phi, init_info.mean_offset, init_info.std_offset, init_info.num_samples);
fprintf('  最大内部距离 = %.4f m\n', init_info.max_inner_distance);
if isfield(init_info, 'thin_ratio') && init_info.thin_ratio > 0
    fprintf('  薄壁警告：%.2f%% 的材料单元到边界距离小于 Δφ\n', 100 * init_info.thin_ratio);
end

% fprintf('边界偏移初始化完成，按任意键继续...\n');
% pause;
enhanced_visualization_check(lsf, material_mask, parallel_paths, struc, nelx, nely, Lx, Ly, dx, dy, delta_phi, init_info);
% === 优化1.2：历史数组预分配（消除动态扩容） ===
% 预分配固定大小数组，循环结束后裁剪到实际长度
compliance_history = zeros(max_iter, 1);
FCS_history = zeros(max_iter, 1);

theta_old = zeros(nely, nelx);
for iter = 1:max_iter
    % === 优化1.3：预计算常用掩码（每次迭代开始） ===
    abs_lsf = abs(lsf);  % 只计算一次
    
    % 统一管理所有频繁使用的掩码
    bands = struct();
    bands.narrow_05h = abs_lsf <= 0.5 * h_grid;  % 零线附近
    bands.narrow_10h = abs_lsf <= 1.0 * h_grid;  % 标准窄带
    bands.narrow_15h = abs_lsf <= 1.5 * h_grid;  % 宽窄带
    % 3.1 根据水平集计算纤维角度
    theta_new = compute_fiber_angles_from_lsf(lsf, dx, dy);

    if iter == 1
        center_i = round((nely+2)/2);
        center_j = round((nelx+2)/2);
        dphi_dx_center = (lsf(center_i, center_j+1) - lsf(center_i, center_j-1)) / (2*dx);
        dphi_dy_center = (lsf(center_i+1, center_j) - lsf(center_i-1, center_j)) / (2*dy);
        grad_mag = hypot(dphi_dx_center, dphi_dy_center);
        fprintf('  中心位置梯度：dx=%e，dy=%e，|grad|=%e\n', dphi_dx_center, dphi_dy_center, grad_mag);
    end

    if iter > 1
        max_change = delta_theta_max;
        theta_change = theta_new - theta_old;
        theta_change = atan2(sin(theta_change), cos(theta_change));
        theta_change = sign(theta_change) .* min(abs(theta_change), max_change);
        theta_e = theta_old + theta_change;
    else
        theta_e = theta_new;
        theta_old = zeros(nely, nelx);
    end
    
    % === 优化1.1：矢量化角度平滑（替代三重嵌套循环，加速80%） ===
    % 参考：解决方案-分点总结.txt B1 + 优化方案阶段1.1
    eta = 0.10;  % 平滑系数 0.05~0.15
    z = angle_smooth_vectorized(theta_e, eta, 2);
    
    % 回到角度
    theta_smooth = 0.5 * angle(z);
    theta_smooth = mod(theta_smooth, pi);  % 归到[0,π)
    % 使用平滑后的角度（只影响FE用的θ，不改φ）
    theta_e = theta_smooth;

    % 3.2 悬臂梁有限元分析（返回统一的载荷向量F）
    % === 阶段3：异常处理 ===
    try
    [U, K, F] = FE_analysis_cantilever(nelx, nely, theta_e, E_L, E_T, nu_LT, G_LT, thickness, F_mag, dx, dy);
    catch ME
        log_message('ERROR', params, '有限元分析失败（迭代 %d）: %s', iter, ME.message);
        
        % 保存错误状态（如启用检查点）
        if params.debug.save_checkpoints
            error_file = sprintf('checkpoints/error_state_iter_%04d.mat', iter);
            if ~exist('checkpoints', 'dir')
                mkdir('checkpoints');
            end
            save(error_file, 'lsf', 'theta_e', 'iter', 'ME');
            log_message('INFO', params, '错误状态已保存: %s', error_file);
        end
        
        rethrow(ME);
    end

    % 3.3 计算柔度（统一使用FE返回的F，确保C = U'*F = U'*K*U > 0）
    compliance = full(U' * F);
    compliance_history(iter) = compliance;  % 优化1.2：直接索引赋值

    % === 能量一致性自检（方案A）===
    % 验证 U'*F = U'*K*U（能量原理）
    if iter == 1 || mod(iter, 10) == 0
        compliance_Ck = full(U' * K * U);
        energy_error = abs(compliance - compliance_Ck);
        relative_error = energy_error / max(1e-12, abs(compliance));
        
        if relative_error > 1e-6
            warning('迭代%d: U''*F 与 U''*K*U 不一致 (误差=%.3e, 相对误差=%.2e)', ...
                iter, energy_error, relative_error);
        end
        
        fprintf('  [能量自检] U''*F = %.6e, U''*K*U = %.6e, 差异 = %.2e\n', ...
            compliance, compliance_Ck, energy_error);
    end

    if iter == 1
        fprintf('\n=== 优化诊断 ===\n');
        fprintf('初始柔度: %e\n', compliance);
        fprintf('施加载荷大小: %e\n', F_mag);
    end

    if iter > 1
        fprintf('  角度变化：最大=%.2f 度，平均=%.2f 度\n', ...
            max(abs(theta_e(:)-theta_old(:))) * 180/pi, ...
            mean(abs(theta_e(:)-theta_old(:))) * 180/pi);
    end

    % 3.4 应变能密度
    strain_energy = compute_strain_energy(nelx, nely, U, theta_e, E_L, E_T, nu_LT, G_LT, thickness, dx, dy);

    % 3.5 灵敏度与速度场
    sensitivity = compute_sensitivity_adjoint(nelx, nely, U, theta_e, E_L, E_T, nu_LT, G_LT, thickness, dx, dy);

    fprintf('迭代 %d - 灵敏度统计：最大=%e，最小=%e，平均=%e\n', iter, ...
        max(sensitivity(:)), min(sensitivity(:)), mean(abs(sensitivity(:))));

    node_sensitivity = aggregate_node_sensitivity(sensitivity, theta_e, lsf, nelx, nely, dx, dy);
    
    % === 修改20-A2：改进lambda_fid计算（基于dev95） ===
    % 参考：解决方案-分点总结.txt A2
    % dev95比均值更抗噪，基于v_shape_target计算更合理
    band_est = bands.narrow_10h;  % 优化1.3：使用预计算掩码
    vshape_est = -node_sensitivity;
    if any(band_est(:))
        dev = abs(lsf(band_est) - lsf_target_global(band_est));
        dev95 = prctile(dev, 95);  % 改用95%分位，更抗噪
        v_shape_target = 3.0;  % 前期目标幅值（可改为动态）
        factor = 2.5;  % 目标：Vfid95 >= 2.5*Vshape95
        lambda_fid = factor * v_shape_target / max(dev95, 1e-12);
        lambda_fid = min(lambda_fid, 5000);  % 上限放宽到5000
    else
        lambda_fid = 50.0;  % 回退固定值
    end
    gamma_curv = 0.5 * h_grid;  % 优化1.3：使用预计算的h_grid
    
    % === 强力稳住4：动态v_target（前期慢稳）===
    % build_velocity_field函数内部根据iter动态调整v_target
    % 前100步v_target=3，后期v_target=5
    [velocity, velocity_stats] = build_velocity_field(node_sensitivity, lsf, dx, dy, 1.5*h_grid, true, lsf_target_global, lambda_fid, [], gamma_curv, iter);
    
    % 诊断：检查梯度模和偏离量
    if iter == 1 || mod(iter, 10) == 0
        [grad_y, grad_x] = gradient(lsf, dy, dx);
        grad_mag = hypot(grad_x, grad_y);
        deviation = abs(lsf - lsf_target_global);
        fprintf('  [诊断] 梯度模：最小=%.3f，最大=%.3f，平均=%.3f\n', ...
            min(grad_mag(:)), max(grad_mag(:)), mean(grad_mag(:)));
        fprintf('  [诊断] 偏离量：最小=%.4f，最大=%.4f，平均=%.4f\n', ...
            min(deviation(:)), max(deviation(:)), mean(deviation(:)));
        
        % === 修改14.3：添加node_sensitivity和速度场诊断 ===
        fprintf('  [诊断] node_sensitivity：最小=%.2e，最大=%.2e，平均=%.2e\n', ...
            min(node_sensitivity(:)), max(node_sensitivity(:)), mean(abs(node_sensitivity(:))));
        fprintf('  [诊断] 速度场max_band=%.2e（包含v_fid）\n', velocity_stats.max_band);
        
        % 估算v_fid幅值（用于验证约束强度）
        if any(bands.narrow_15h(:))  % 优化1.3：使用预计算掩码
            deviation_band = deviation(bands.narrow_15h);
            v_fid_estimated = lambda_fid * mean(deviation_band);
            fprintf('  [诊断] 估算v_fid平均幅值=%.2e（lambda_fid=%.1f）\n', v_fid_estimated, lambda_fid);
        end
        
        % 路径一致性统计（在φ=0附近±0.5h范围）
        zero_band = bands.narrow_05h;  % 优化1.3：使用预计算掩码
        if any(zero_band(:))
            deviation_zero = deviation(zero_band);
            consistency_ratio = sum(deviation_zero < h) / numel(deviation_zero) * 100;
            fprintf('  [诊断] 路径一致性：零线附近%.1f%%节点偏离<1网格（平均偏离=%.4f）\n', ...
                consistency_ratio, mean(deviation_zero));
        end
        
        % === Vfid/Vshape驱动对比诊断（清单四-1）===
        % 参考：修改与保留清单.txt 四-1)
        % 验证约束驱动强度是否>=1.5
        band_chk = bands.narrow_10h;  % 优化1.3：使用预计算掩码
        if any(band_chk(:))
            vshape_band = -node_sensitivity(band_chk);
            deviation_band_chk = lsf(band_chk) - lsf_target_global(band_chk);
            vfid_band = -lambda_fid * deviation_band_chk;
            
            Vshape = prctile(abs(vshape_band), 95);  % 与lambda_fid计算一致
            Vfid = prctile(abs(vfid_band), 95);
            ratio = Vfid / max(Vshape, 1e-12);
            
            fprintf('  [驱动对比] Vfid/Vshape = %.2f (目标≥2.5, lambda_fid=%.1f)\n', ratio, lambda_fid);
            
            if ratio < 1.5
                fprintf('  ⚠️  约束驱动严重不足！\n');
            elseif ratio < 2.5
                fprintf('  ⚠  约束驱动偏弱，建议提升lambda_fid系数\n');
            else
                fprintf('  ✓  约束驱动充足（强力稳住）\n');
            end
        end
    end

    dt_adaptive = compute_adaptive_timestep(velocity, dx, dy);
    if velocity_stats.max_band > 1e-12
        dt_angle = delta_theta_max / velocity_stats.max_band;
        dt_adaptive = min(dt_adaptive, dt_angle);
    end
    max_node_sens = max(abs(node_sensitivity(:)));
    if max_node_sens > 1e-12
        dt_angle_sens = delta_theta_max / max_node_sens;
        dt_adaptive = min(dt_adaptive, dt_angle_sens);
    end

    if iter == 1 || mod(iter, 10) == 0
        fprintf('  节点灵敏度统计：最大=%.3e，最小=%.3e，平均=%.3e\n', ...
            max(node_sensitivity(:)), min(node_sensitivity(:)), mean(abs(node_sensitivity(:))));
        fprintf('  速度场统计：最大=%.3e，最小=%.3e，平均=%.3e\n', ...
            max(velocity(:)), min(velocity(:)), mean(abs(velocity(:))));
        fprintf('  窄带速度：平均=%.3e，加权平均=%.3e，最大=%.3e，正向=%.1f%%，负向=%.1f%%\n', ...
            velocity_stats.mean_band, velocity_stats.weighted_mean, velocity_stats.max_band, ...
            velocity_stats.pos_ratio*100, velocity_stats.neg_ratio*100);
        fprintf('  界面长度=%.3f，质心=(%.3f, %.3f)\n', ...
            velocity_stats.length, velocity_stats.centroid(1), velocity_stats.centroid(2));
        predicted_change = velocity_stats.max_band * dt_adaptive;
        fprintf('  预测最大角度变化 %.2f 度\n', predicted_change * 180/pi);
        fprintf('  自适应时间步长: %.6f\n', dt_adaptive);
        diag_report(iter, node_sensitivity, velocity_stats);
    end

    % 3.6 更新水平集
    lsf_before = lsf;
    lsf = update_levelset_HJ(lsf, velocity, dt_adaptive, dx, dy);
    
    % === 阶段3：水平集更新验证（检查NaN/Inf） ===
    if any(~isfinite(lsf(:)))
        log_message('ERROR', params, '水平集出现NaN/Inf（迭代 %d）', iter);
        lsf = lsf_before;  % 回退到更新前状态
        log_message('WARN', params, '已回退到更新前状态');
    end

    % === 强力稳住2：投影更硬（前期加强）===
    % 参考：强力稳住-总结清单.txt 一-2)
    % 前100步强力拉回，后期温和约束
    ENABLE_HARD_PROJECTION = true;
    
    if iter <= 100
        omega_proj = 0.7;  % 前100步：更强投影（从0.5提升）
        proj_band = bands.narrow_15h;  % 优化1.3：前100步更宽带
    else
        omega_proj = 0.5;  % 后期：正常投影
        proj_band = bands.narrow_10h;  % 优化1.3：后期标准窄带
    end
    
    if ENABLE_HARD_PROJECTION
        deviation_proj = lsf(proj_band) - lsf_target_global(proj_band);
        lsf(proj_band) = lsf(proj_band) - omega_proj * deviation_proj;
        
        if iter == 1 || mod(iter, 10) == 0
            max_correction = max(abs(omega_proj * deviation_proj));
            phase_str = '是';
            if iter > 100
                phase_str = '否';
            end
            fprintf('  [硬约束投影] omega=%.2f, 最大修正=%.4f (前期=%s)\n', ...
                omega_proj, max_correction, phase_str);
        end
    end

    lsf_change = max(abs(lsf(:) - lsf_before(:)));
    if mod(iter, 10) == 0
        fprintf('  水平集变化: %.3e\n', lsf_change);
    end

    % === 阶段2：自适应重初始化（智能判断） ===
    % 初始化计数器（第一次迭代）
    if iter == 1
        iter_since_last_reinit = 0;
    end
    
    iter_since_last_reinit = iter_since_last_reinit + 1;
    [do_reinit, reinit_reason] = should_reinitialize(lsf, lsf_before, ...
        iter_since_last_reinit, iter, params);
    
    if do_reinit
        log_message('INFO', params, '触发重初始化: %s', reinit_reason);
        zero_mask_dynamic = compute_zero_mask_from_lsf(lsf, h_grid);
        
        % === 重初始化域选择（刻意设计） ===
        % 设计意图：
        %   1. 速度场演化仅在材料边界环带更新零等值线（主纤维路径）
        %   2. FMM重初始化需将距离信息传播到全域（包括void区域）
        %   3. 全域传播确保FE分析所需的纤维角度场θ(x,y)连续无奇异
        % 默认：不传material_mask（即全域），避免后续误改
        % 可选：通过config.levelset.reinit_domain='masked'限制到材料域（仅供试验）
        if strcmp(params.levelset.reinit_domain, 'masked')
            reinit_mask = material_mask;
        else
            reinit_mask = [];  % 全域传播（默认）
        end
        lsf = fmm_reinitialize(lsf, dx, dy, zero_mask_dynamic, reinit_mask);
        iter_since_last_reinit = 0;  % 重置计数
    end

    % 3.8 纤维连续性评分
    FCS = compute_fiber_continuity(theta_e);
    FCS_history(iter) = FCS;  % 优化1.2：直接索引赋值
    theta_old = theta_e;
    
    % === 阶段3：检查点保存 ===
    save_checkpoint(iter, lsf, theta_e, compliance_history, FCS_history, params);
    
    % === 修改20验收诊断（每10步打印） ===
    if mod(iter, 10) == 0
        band = bands.narrow_10h;  % 优化1.3：使用预计算掩码
        if any(band(:))
            vshape_band = vshape_est(band);  % 使用原始node_sensitivity
            vfid_band = -lambda_fid * (lsf(band) - lsf_target_global(band));
            Vshape95 = prctile(abs(vshape_band), 95);
            Vfid95 = prctile(abs(vfid_band), 95);
            dev95 = prctile(abs(lsf(band) - lsf_target_global(band)), 95);
            mean_off = mean(abs(lsf(band) - lsf_target_global(band)));
            
            fprintf('  ✅ [验收] Vfid95/Vshape95=%.2f  dev95=%.3e  mean_off=%.3e  FCS=%.1f%%\n', ...
                Vfid95/max(Vshape95,1e-12), dev95, mean_off, FCS*100);
            
            % 达标检查
            if Vfid95/max(Vshape95,1e-12) >= 2.5 && dev95 < 1.0*h_grid && FCS >= 0.80  % 优化1.3
                fprintf('  🎉 已达标！（Vfid95/Vshape95≥2.5, dev95<h, FCS≥80%%）\n');
            end
        end
    end

    % 3.9 收敛判据（柔度滚动标准差）
    if iter > 50
        recent = compliance_history(max(1, end-19):end);
        rel_change = std(recent)/mean(recent);
        if rel_change < tol
            fprintf('优化在第 %d 次迭代时收敛\n', iter);
            break;
        end
    end

    if iter == 1
        fprintf('\n=== 系统诊断 ===\n');
        fprintf('材料各向异性 E_L/E_T = %.1f\n', E_L / E_T);
        fprintf('网格: %dx%d，单元尺寸 %.3fx%.3f\n', nelx, nely, dx, dy);
        zero_angle_elements = sum(abs(theta_e(:)) < 0.01);
        fprintf('初始接近零角度单元数: %d\n', zero_angle_elements);
        test_angle = 0;
        c = cos(test_angle);
        s = sin(test_angle);
        fprintf('当 θ=0 时: cos=%.3f，sin=%.3f，-2cs=%.3f\n', c, s, -2*c*s);
    end

    if mod(iter, 5) == 0
        fprintf('\n迭代 %d: 柔度 = %.4e，FCS = %.2f%%\n', iter, compliance, FCS*100);
        fprintf('  角度统计：最小=%.1f 度，最大=%.1f 度，平均=%.1f 度\n', ...
            min(theta_e(:))*180/pi, max(theta_e(:))*180/pi, mean(theta_e(:))*180/pi);
        fprintf('  接近 0 度单元: %d，接近 90 度单元: %d\n', ...
            sum(abs(theta_e(:)) < 0.1), sum(abs(theta_e(:) - pi/2) < 0.1));
    end
end

% === 优化1.2：裁剪历史数组到实际迭代次数 ===
final_iter = iter;
compliance_history = compliance_history(1:final_iter);
FCS_history = FCS_history(1:final_iter);

%% 4. 可视化与结果汇报
plot_if_enabled(params, @() visualize_results_article(lsf, theta_e, strain_energy, ...
    compliance_history, FCS_history, nelx, nely, Lx, Ly, dx, dy));

log_message('INFO', params, '\n优化完成!');
log_message('INFO', params, '最终柔度: %.4e', compliance);
log_message('INFO', params, '纤维连续性评分: %.2f%%', FCS*100);

% 柔度降低比例（方案A修复后应恒为正向下降）
baseC = compliance_history(1);
currC = compliance;
if baseC > 0
    improve_ratio = (baseC - currC) / baseC * 100;
    log_message('INFO', params, '柔度降低比例: %.2f%% (初始=%.4e, 最终=%.4e)', ...
        improve_ratio, baseC, currC);
else
    log_message('WARN', params, '警告：初始柔度<=0，无法计算降低比例');
end

% === 构建结果结构体 ===
results = struct();
results.lsf = lsf;
results.theta_e = theta_e;
results.compliance_history = compliance_history;
results.FCS_history = FCS_history;
results.final_compliance = compliance;
results.final_FCS = FCS;
results.final_iter = final_iter;
results.strain_energy = strain_energy;
results.params = params;

if baseC > 0
    results.improvement_ratio = improve_ratio;
else
    results.improvement_ratio = NaN;
end

log_message('INFO', params, '结果已保存至输出结构体');

end
