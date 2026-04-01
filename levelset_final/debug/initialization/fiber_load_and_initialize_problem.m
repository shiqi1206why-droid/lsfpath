function problem_ctx = fiber_load_and_initialize_problem(runtime_ctx)
%FIBER_LOAD_AND_INITIALIZE_PROBLEM Load topology and build initial optimization state.

    params = runtime_ctx.params;
    paths = runtime_ctx.paths;

    nelx = runtime_ctx.nelx;
    nely = runtime_ctx.nely;
    Lx = runtime_ctx.Lx;
    Ly = runtime_ctx.Ly;
    dx = runtime_ctx.dx;
    dy = runtime_ctx.dy;

    E_L = runtime_ctx.E_L;
    E_T = runtime_ctx.E_T;
    nu_LT = runtime_ctx.nu_LT;
    G_LT = runtime_ctx.G_LT;
    thickness = runtime_ctx.thickness;
    F_mag = runtime_ctx.F_mag;
    delta_theta_max = runtime_ctx.delta_theta_max;

    delta_phi = runtime_ctx.delta_phi;
    init_smooth_opts = runtime_ctx.init_smooth_opts;

    global DIAG; %#ok<GVMIS>
    DIAG = struct();
    diag_reset();

    topo_file = paths.topology_file;
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

    log_message('INFO', params, '开始执行边界偏移初始化...');
    enable_init_plots = isfield(params, 'debug') && isfield(params.debug, 'enable_plots') && params.debug.enable_plots;
    if enable_init_plots
        figure('Name', '拓扑与初始化检查', 'Position', [100, 100, 1200, 400]);

        subplot(1,3,1);
        imagesc(struc);
        colormap(gray);
        axis equal; axis tight;
        title('原始拓扑');
        xlabel('x方向单元索引');
        ylabel('y方向单元索引');
    end

    log_message('INFO', params, '正在清理材料掩膜...');
    [material_mask_core, mask_info] = clean_material_mask(struc, params.init.min_component_size, init_smooth_opts.morph_radius);
    material_mask_full = expand_material_mask_to_full(material_mask_core);

    if enable_init_plots
        subplot(1,3,2);
        imagesc(material_mask_core);
        colormap(gray);
        axis equal; axis tight;
        title('清理后的材料掩膜');
        xlabel('x方向单元索引');
        ylabel('y方向单元索引');
        hold on;
        boundary_mask = bwperim(material_mask_core);
        [boundary_y, boundary_x] = find(boundary_mask);
        plot(boundary_x, boundary_y, 'r.', 'MarkerSize', 2, 'DisplayName', '材料边界');
        if ~isempty(boundary_x)
            legend('材料边界', 'Location', 'best');
        end
    end
    log_message('DEBUG', params, '  连通区域数: %d，总像素数: %d', mask_info.num_components, mask_info.total_area);

    log_message('INFO', params, '正在构建基于边界偏移的符号距离场...');
    [lsf, parallel_paths, init_info] = construct_boundary_offset_levelset_with_parallel( ...
        material_mask_core, nelx, nely, dx, dy, delta_phi, init_smooth_opts);

    initial_zero_mask = init_info.zero_mask;
    lsf_initial = lsf;

    log_message('INFO', params, '正在构建边界等距目标场...');
    phi_boundary_global = init_info.phi_boundary_full;

    h_grid = min(dx, dy);
    h = h_grid;
    boundary_guard_band = abs(init_info.phi_boundary_full) <= params.levelset.boundary_guard_band_factor * h_grid;
    boundary_guard_band = boundary_guard_band & material_mask_full;
    lsf_target_global = phi_boundary_global + init_info.delta_phi_used;
    log_message('INFO', params, '  边界等距目标场已构建（用于优化约束，Δφ_used=%.4f）', init_info.delta_phi_used);

    if isfield(init_info, 'contour') && ~isempty(init_info.contour.x)
        fprintf('  [调试] init_info.contour采样点数量: %d\n', length(init_info.contour.x));
        fprintf('  [调试] X范围: [%.3f, %.3f] m\n', min(init_info.contour.x), max(init_info.contour.x));
        fprintf('  [调试] Y范围: [%.3f, %.3f] m\n', min(init_info.contour.y), max(init_info.contour.y));
    else
        fprintf('  [警告] init_info.contour为空或不存在！\n');
    end

    if enable_init_plots
        subplot(1,3,3);
        [X_idx, Y_idx] = meshgrid(0:nelx+1, 0:nely+1);
        lsf_plot = lsf;
        lsf_plot(~material_mask_full) = NaN;
        contour(X_idx, Y_idx, lsf_plot, 20, 'LineWidth', 0.5, 'DisplayName', '等值线');
        hold on;

        C_zero = contourc(0:nelx+1, 0:nely+1, lsf_plot, [0 0]);
        [zero_x, zero_y] = contourc_to_points(C_zero);
        if ~isempty(zero_x)
            fprintf('  [调试] 直接提取的零等值线采样点数量: %d\n', length(zero_x));
            plot(zero_x, zero_y, 'r-', 'LineWidth', 2, 'DisplayName', '主路径 \phi=0');
            plot(zero_x, zero_y, 'mo', 'MarkerSize', 3, 'MarkerFaceColor', 'm', ...
                'LineStyle', 'none', 'DisplayName', '主路径采样点');
        else
            fprintf('  [警告] 直接提取的零等值线也为空！\n');
        end

        axis equal; axis tight;
        set(gca,'YDir','reverse');
        title('初始化水平集等值线');
        xlabel('x方向单元索引');
        ylabel('y方向单元索引');
        legend('Location','best');
        colorbar;
    end

    fprintf('  目标Δφ = %.4f m，抽样均值 = %.4f m (标准差 = %.4f m，样本数 = %d)\n', ...
        delta_phi, init_info.mean_offset, init_info.std_offset, init_info.num_samples);
    fprintf('  最大内部距离 = %.4f m\n', init_info.max_inner_distance);
    if isfield(init_info, 'thin_ratio') && init_info.thin_ratio > 0
        fprintf('  薄壁警告：%.2f%% 的材料单元到边界距离小于 Δφ\n', 100 * init_info.thin_ratio);
    end

    enhanced_visualization_check(lsf, material_mask_core, parallel_paths, struc, nelx, nely, Lx, Ly, dx, dy, delta_phi, init_info, enable_init_plots);

    max_iter = runtime_ctx.max_iter;
    history_capacity = max_iter + 1;
    compliance_history = nan(history_capacity, 1);
    FCS_history = nan(history_capacity, 1);
    theta_only_compliance_history = nan(max_iter, 1);
    hj_raw_compliance_history = nan(max_iter, 1);
    hj_trial_compliance_history = nan(max_iter, 1);
    hj_local_reinit_compliance_history = nan(max_iter, 1);
    reinit_trial_compliance_history = nan(max_iter, 1);
    accepted_source_history = cell(max_iter, 1);
    accepted_source_detail_history = strings(max_iter, 1);
    accepted_steps = 0;
    rejected_steps = 0;
    theta_only_accept_count = 0;
    theta_only_reject_count = 0;
    theta_only_consecutive_count = 0;
    reject_due_next_guard = 0;
    reject_due_current_guard = 0;
    reinit_trigger_count = 0;
    reinit_skip_due_reject = 0;
    reinit_skip_due_objective = 0;
    reinit_skip_due_next_guard = 0;
    reinit_skip_due_current_guard = 0;
    iter_since_last_reinit = 0;
    best_state = struct('compliance', inf, 'iter', 0, 'lsf', [], 'theta', [], ...
        'strain_energy', [], 'FCS', NaN, 'theta_target', []);
    no_improve_counter = 0;
    early_stop_triggered = false;
    early_stop_reason = '';
    history_count = 0;
    current_state_recorded = false;
    prev_theta_for_diag = [];
    candidate_select_tol = max(runtime_ctx.best_state_rel_tol, 1e-12);
    loop_iter_count = 0;
    raw_final_compliance = NaN;
    raw_final_FCS = NaN;
    raw_final_improvement_ratio = NaN;
    final_to_best_gap_percent = NaN;
    quality_eval_stride = 5;
    raw_turn_history = nan(history_capacity, 1);
    raw_kappa_history = nan(history_capacity, 1);
    raw_spacing_error_history = nan(history_capacity, 1);
    raw_grad_dev_history = nan(history_capacity, 1);
    raw_near_zero_outlier_history = nan(history_capacity, 1);
    hj_fallback_history = nan(max_iter, 1);
    hj_second_order_history = nan(max_iter, 1);
    hj_frozen_incomplete_history = nan(max_iter, 1);
    hj_first_order_complete_history = nan(max_iter, 1);
    reinit_method_history = strings(max_iter, 1);
    reinit_fallback_history = false(max_iter, 1);
    reinit_reason_history = strings(max_iter, 1);
    hj_update_diagnostics_history = cell(max_iter, 1);
    reinit_diagnostics_history = cell(max_iter, 1);
    boundary_guard_ratio_history = nan(max_iter, 1);
    frozen_boundary_point_history = nan(max_iter, 1);
    accepted_hj_reinit_count = 0;
    local_reinit_shell_size_history = nan(max_iter, 1);
    post_reinit_grad_dev_mean_history = nan(max_iter, 1);
    post_reinit_grad_outlier_ratio_history = nan(max_iter, 1);
    refresh_shell_size_history = nan(max_iter, 1);
    refresh_count = 0;
    last_velocity_field = zeros(size(lsf));
    last_propagation_mask = false(size(lsf));
    last_hj_info = struct();
    last_reinit_info = struct();
    gradient_chain_cosine_history = nan(max_iter, 1);
    gradient_chain_norm_ratio_history = nan(max_iter, 1);
    gradient_chain_topk_sign_history = nan(max_iter, 1);
    gradient_chain_saturation_history = nan(max_iter, 1);
    gradient_chain_degenerate_history = nan(max_iter, 1);
    gradient_chain_active_band_coverage_history = nan(max_iter, 1);
    gradient_chain_zero_limiter_history = nan(max_iter, 1);
    gradient_chain_exact_nonzero_history = nan(max_iter, 1);
    gradient_chain_support_overlap_history = nan(max_iter, 1);
    gradient_chain_theta_raw_guard_history = nan(max_iter, 1);
    gradient_chain_theta_raw_guard_overlap_history = nan(max_iter, 1);
    gradient_chain_full_vs_opt_overlap_history = nan(max_iter, 1);
    gradient_chain_selected_source_history = strings(max_iter, 1);
    gradient_chain_audit_mode = string(params.gradient.audit_support_mode);
    last_gradient_chain_diag = struct();
    manufacturing_grad_norm_history = nan(max_iter, 1);
    manufacturing_curvature_norm_history = nan(max_iter, 1);
    manufacturing_gap_overlap_norm_history = nan(max_iter, 1);
    theta_only_vs_current_history = nan(max_iter, 1);
    hj_raw_vs_theta_only_history = nan(max_iter, 1);
    reinit_vs_theta_only_history = nan(max_iter, 1);
    last_manufacturing_diag = struct();

    current_state = evaluate_candidate_state(lsf, [], delta_theta_max, dx, dy, ...
        nelx, nely, material_mask_core, E_L, E_T, nu_LT, G_LT, thickness, F_mag, ...
        params.smooth.eta, params.smooth.iterations);

    problem_ctx = struct();
    problem_ctx.struc = struc;
    problem_ctx.material_mask_core = material_mask_core;
    problem_ctx.material_mask_full = material_mask_full;
    problem_ctx.mask_info = mask_info;
    problem_ctx.lsf = lsf;
    problem_ctx.parallel_paths = parallel_paths;
    problem_ctx.init_info = init_info;
    problem_ctx.initial_zero_mask = initial_zero_mask;
    problem_ctx.lsf_initial = lsf_initial;
    problem_ctx.phi_boundary_global = phi_boundary_global;
    problem_ctx.h_grid = h_grid;
    problem_ctx.h = h;
    problem_ctx.boundary_guard_band = boundary_guard_band;
    problem_ctx.lsf_target_global = lsf_target_global;

    problem_ctx.history_capacity = history_capacity;
    problem_ctx.compliance_history = compliance_history;
    problem_ctx.FCS_history = FCS_history;
    problem_ctx.theta_only_compliance_history = theta_only_compliance_history;
    problem_ctx.hj_raw_compliance_history = hj_raw_compliance_history;
    problem_ctx.hj_trial_compliance_history = hj_trial_compliance_history;
    problem_ctx.hj_local_reinit_compliance_history = hj_local_reinit_compliance_history;
    problem_ctx.reinit_trial_compliance_history = reinit_trial_compliance_history;
    problem_ctx.accepted_source_history = accepted_source_history;
    problem_ctx.accepted_source_detail_history = accepted_source_detail_history;
    problem_ctx.accepted_steps = accepted_steps;
    problem_ctx.rejected_steps = rejected_steps;
    problem_ctx.theta_only_accept_count = theta_only_accept_count;
    problem_ctx.theta_only_reject_count = theta_only_reject_count;
    problem_ctx.theta_only_consecutive_count = theta_only_consecutive_count;
    problem_ctx.reject_due_next_guard = reject_due_next_guard;
    problem_ctx.reject_due_current_guard = reject_due_current_guard;
    problem_ctx.reinit_trigger_count = reinit_trigger_count;
    problem_ctx.reinit_skip_due_reject = reinit_skip_due_reject;
    problem_ctx.reinit_skip_due_objective = reinit_skip_due_objective;
    problem_ctx.reinit_skip_due_next_guard = reinit_skip_due_next_guard;
    problem_ctx.reinit_skip_due_current_guard = reinit_skip_due_current_guard;
    problem_ctx.iter_since_last_reinit = iter_since_last_reinit;
    problem_ctx.best_state = best_state;
    problem_ctx.no_improve_counter = no_improve_counter;
    problem_ctx.early_stop_triggered = early_stop_triggered;
    problem_ctx.early_stop_reason = early_stop_reason;
    problem_ctx.history_count = history_count;
    problem_ctx.current_state_recorded = current_state_recorded;
    problem_ctx.prev_theta_for_diag = prev_theta_for_diag;
    problem_ctx.candidate_select_tol = candidate_select_tol;
    problem_ctx.loop_iter_count = loop_iter_count;
    problem_ctx.raw_final_compliance = raw_final_compliance;
    problem_ctx.raw_final_FCS = raw_final_FCS;
    problem_ctx.raw_final_improvement_ratio = raw_final_improvement_ratio;
    problem_ctx.final_to_best_gap_percent = final_to_best_gap_percent;
    problem_ctx.quality_eval_stride = quality_eval_stride;
    problem_ctx.raw_turn_history = raw_turn_history;
    problem_ctx.raw_kappa_history = raw_kappa_history;
    problem_ctx.raw_spacing_error_history = raw_spacing_error_history;
    problem_ctx.raw_grad_dev_history = raw_grad_dev_history;
    problem_ctx.raw_near_zero_outlier_history = raw_near_zero_outlier_history;
    problem_ctx.hj_fallback_history = hj_fallback_history;
    problem_ctx.hj_second_order_history = hj_second_order_history;
    problem_ctx.hj_frozen_incomplete_history = hj_frozen_incomplete_history;
    problem_ctx.hj_first_order_complete_history = hj_first_order_complete_history;
    problem_ctx.reinit_method_history = reinit_method_history;
    problem_ctx.reinit_fallback_history = reinit_fallback_history;
    problem_ctx.reinit_reason_history = reinit_reason_history;
    problem_ctx.hj_update_diagnostics_history = hj_update_diagnostics_history;
    problem_ctx.reinit_diagnostics_history = reinit_diagnostics_history;
    problem_ctx.boundary_guard_ratio_history = boundary_guard_ratio_history;
    problem_ctx.frozen_boundary_point_history = frozen_boundary_point_history;
    problem_ctx.accepted_hj_reinit_count = accepted_hj_reinit_count;
    problem_ctx.local_reinit_shell_size_history = local_reinit_shell_size_history;
    problem_ctx.post_reinit_grad_dev_mean_history = post_reinit_grad_dev_mean_history;
    problem_ctx.post_reinit_grad_outlier_ratio_history = post_reinit_grad_outlier_ratio_history;
    problem_ctx.refresh_shell_size_history = refresh_shell_size_history;
    problem_ctx.refresh_count = refresh_count;
    problem_ctx.last_velocity_field = last_velocity_field;
    problem_ctx.last_propagation_mask = last_propagation_mask;
    problem_ctx.last_hj_info = last_hj_info;
    problem_ctx.last_reinit_info = last_reinit_info;
    problem_ctx.gradient_chain_cosine_history = gradient_chain_cosine_history;
    problem_ctx.gradient_chain_norm_ratio_history = gradient_chain_norm_ratio_history;
    problem_ctx.gradient_chain_topk_sign_history = gradient_chain_topk_sign_history;
    problem_ctx.gradient_chain_saturation_history = gradient_chain_saturation_history;
    problem_ctx.gradient_chain_degenerate_history = gradient_chain_degenerate_history;
    problem_ctx.gradient_chain_active_band_coverage_history = gradient_chain_active_band_coverage_history;
    problem_ctx.gradient_chain_zero_limiter_history = gradient_chain_zero_limiter_history;
    problem_ctx.gradient_chain_exact_nonzero_history = gradient_chain_exact_nonzero_history;
    problem_ctx.gradient_chain_support_overlap_history = gradient_chain_support_overlap_history;
    problem_ctx.gradient_chain_theta_raw_guard_history = gradient_chain_theta_raw_guard_history;
    problem_ctx.gradient_chain_theta_raw_guard_overlap_history = gradient_chain_theta_raw_guard_overlap_history;
    problem_ctx.gradient_chain_full_vs_opt_overlap_history = gradient_chain_full_vs_opt_overlap_history;
    problem_ctx.gradient_chain_selected_source_history = gradient_chain_selected_source_history;
    problem_ctx.gradient_chain_audit_mode = gradient_chain_audit_mode;
    problem_ctx.last_gradient_chain_diag = last_gradient_chain_diag;
    problem_ctx.manufacturing_grad_norm_history = manufacturing_grad_norm_history;
    problem_ctx.manufacturing_curvature_norm_history = manufacturing_curvature_norm_history;
    problem_ctx.manufacturing_gap_overlap_norm_history = manufacturing_gap_overlap_norm_history;
    problem_ctx.theta_only_vs_current_history = theta_only_vs_current_history;
    problem_ctx.hj_raw_vs_theta_only_history = hj_raw_vs_theta_only_history;
    problem_ctx.reinit_vs_theta_only_history = reinit_vs_theta_only_history;
    problem_ctx.last_manufacturing_diag = last_manufacturing_diag;
    problem_ctx.current_state = current_state;
end
