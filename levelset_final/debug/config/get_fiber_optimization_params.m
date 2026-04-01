function params = get_fiber_optimization_params(config_name)
    % 纤维路径优化参数配置
    % 
    % 输入：
    %   config_name - 配置名称（可选，默认'default'）
    % 
    % 输出：
    %   params - 参数结构体
    %
    % 支持的配置：
    %   'default' - 标准配置（平衡精度和速度）
    %   'fast'    - 快速模式（降低迭代次数，关闭可视化）
    %   'precise' - 精确模式（更严格的收敛，更多迭代）
    %   'debug'   - 调试模式（最大化诊断信息）
    
    if nargin < 1
        config_name = 'default';
    end
    
    % ========== 基础配置 ==========
    base = struct();
    
    % 网格参数
    base.grid.nelx = 80;
    base.grid.nely = 50;
    base.grid.Lx = 1.6;
    base.grid.Ly = 1.0;
    
    % 材料参数（支持多材料切换）
    base.material = get_material_params('carbon_fiber');
    
    % 优化控制
    base.opt.max_iter = 1000;
    base.opt.tol = 1e-5;
    base.opt.alpha = 0.5;
    base.opt.dt = 0.02;
    base.opt.delta_theta_max_deg = 0.3;
    base.opt.fidelity_weight = 0.02;
    base.opt.normalize_sensitivity = false;      % 控制是否对伴随灵敏度归一化
    base.opt.enable_curvature = false;           % 是否启用曲率正则
    base.opt.enable_step_acceptance = true;      % 启用目标下降保护（单步接受/拒绝）
    base.opt.acceptance_tol = 0;                 % 允许的单步柔度相对上浮容差（0=严格不恶化）
    base.opt.enable_current_state_guard = true;  % 额外约束：候选步不应恶化当前步柔度
    base.opt.current_state_tol = 0;              % 相对当前步柔度容差（0=严格不恶化）
    base.opt.current_guard_start_iter = 1;       % 从第几步开始启用当前步门禁
    base.opt.enable_reinit_current_guard = true; % 重初始化候选也需满足当前步柔度门禁
    base.opt.reinit_current_tol = 0;             % 重初始化相对当前步柔度容差
    base.opt.reinit_guard_start_iter = 1;        % 重初始化当前步门禁启用迭代
    base.opt.backtrack_factor = 0.5;             % 回溯缩放因子
    base.opt.max_backtrack = 3;                  % 最大回溯次数
    base.opt.min_backtrack_dt = 1e-6;            % 回溯最小步长
    base.opt.enable_best_state_guard = true;     % 启用历史最优状态守护
    base.opt.best_state_rel_tol = 1e-6;          % 历史最优更新相对阈值
    base.opt.best_state_patience = 20;           % 连续无改进耐心步数（用于早停）
    base.opt.theta_only_fuse_limit = 8;          % theta_only连续接管熔断阈值
    
    % 水平集参数
    base.levelset.delta_phi_factor = 0.8;       % 边界偏移因子
    base.levelset.bandwidth_factor = 1.5;       % 窄带宽度因子
    base.levelset.transition_iter = 100;        % 前期/后期分界点
    base.levelset.reinit_freq_early = 5;        % 前期重初始化频率
    base.levelset.reinit_freq_late = 10;        % 后期重初始化频率
    base.levelset.reinit_threshold = 0.75;      % 自适应重初始化阈值（相对h）
    base.levelset.reinit_max_interval = 15;     % 最大重初始化间隔
    base.levelset.gradient_deviation_tol = 0.15;% 梯度偏差容差
    base.levelset.reinit_bandwidth_factor = 1.0;% 重初始化判据带宽（|φ|<=band*h）
    base.levelset.reinit_domain = 'masked';     % 重初始化域：'masked'=限制到材料域（默认）,'full'=兼容旧试验
    base.levelset.advection_order = 2;          % HJ空间离散阶数：1/2
    base.levelset.time_integrator = 'ssprk2';   % HJ时间推进：'euler'/'ssprk2'
    base.levelset.reinit_method = 'subcell_signed_distance'; % 重初始化方法
    base.levelset.fallback_first_order = true;  % 二阶格式在退化处回退到一阶
    base.levelset.freeze_on_incomplete_godunov = true; % Godunov一阶 stencil 不完整时冻结该点
    base.levelset.eno_smoothness_factor = 2.5;  % 二阶单边导数光滑性阈值
    base.levelset.stencil_buffer_cells = 2;     % HJ 求导读取壳层厚度（单元数）
    base.levelset.hj_rhs_mode = 'legacy';       % HJ右端计算模式：legacy/indexed/vectorized_first_order
    base.levelset.boundary_guard_band_factor = 1.5; % 材料边界护带宽度（相对h）
    base.levelset.local_reinit_zero_band_factor = 2.0; % 局部重初始化零线带宽（相对h）
    base.levelset.local_reinit_buffer_cells = 3; % 局部重初始化壳层扩展单元数（square=2r+1）
    base.levelset.refresh_interval = 5;         % 带外刷新触发的已接受HJ步间隔
    base.levelset.refresh_buffer_cells = 4;     % 带外刷新壳层扩展单元数（square=2r+1）
    base.levelset.zero_geometry_min_points = 8; % 几何重初始化所需最小折线点数
    base.levelset.zero_geometry_min_length = 0.5; % 相对h的最小零线长度阈值
    
    % 投影参数
    base.projection.enable = false;
    base.projection.omega_early = 0.7;          % 前期投影强度
    base.projection.omega_late = 0.5;           % 后期投影强度
    base.projection.band_factor_early = 1.5;    % 前期投影带宽
    base.projection.band_factor_late = 1.0;     % 后期投影带宽
    
    % 平滑参数
    base.smooth.eta = 0.10;                     % 角度平滑系数
    base.smooth.iterations = 2;                 % 平滑迭代次数

    % 梯度链参数
    base.gradient.chain_mode = 'legacy';        % legacy/shadow/exact
    base.gradient.limiter_mode = 'hard';        % hard/soft_experiment
    base.gradient.audit_support_mode = 'full';  % 审计口径：full/opt
    base.gradient.shadow_topk = 32;             % shadow比较的top-k样本数
    base.gradient.shadow_num_directional_checks = 4; % 方向导数审计向量数
    base.gradient.shadow_fd_eps_factor = 1e-4;  % 有限差分步长系数（乘以h）
    base.gradient.soft_limiter_beta = 20.0;     % soft limiter实验强度
    base.gradient.theta_raw_grad_floor = 0.05;  % exact pullback中theta_raw链的最小可微梯度阈值
    base.gradient.theta_raw_branch_cut_tol = 1e-6; % 审计时theta_raw分支切口不可微保护阈值

    % 制造约束参数（优化态附加项，不影响评估态）
    base.manufacturing.enable = true;
    base.manufacturing.grad_norm_weight = 0.0;
    base.manufacturing.curvature_weight = 0.0;
    base.manufacturing.gap_overlap_weight = 0.0;
    base.manufacturing.curvature_radius_min = 2.0;
    base.manufacturing.gap_overlap_target = 1.0;
    base.manufacturing.penalty_band_factor = 1.5;

    % 速度场参数
    base.velocity.enable_bias_removal = true;   % 是否执行形状项去偏
    base.velocity.bias_beta = 0.10;             % 最终速度场净平移抑制强度
    base.velocity.scale_quantile = 95;          % 节点灵敏度分位缩放
    base.velocity.clip_abs = inf;               % 节点灵敏度裁剪上限（inf=不裁剪）
    
    % 载荷参数
    base.load.F_mag = -1;                       % 载荷大小 (N)
    
    % 初始化参数
    base.init.morph_radius = 1;                 % 形态学半径
    base.init.min_component_size = 10;          % 最小连通区域大小
    base.init.boundary_reconstruction = 'marching_squares_linear'; % 初始化边界重建方法

    % 严格材料域约束策略（论文口径）
    base.constraint.material_constraint_mode = 'strict';
    base.constraint.boundary_contact_policy = 'no-crossing';
    
    % 调试与诊断
    base.debug.verbose = true;                  % 详细输出
    base.debug.log_level = 'INFO';              % 日志级别：DEBUG/INFO/WARN/ERROR
    base.debug.log_interval = 10;               % 日志输出间隔
    base.debug.enable_plots = true;             % 启用绘图
    base.debug.enable_diagnostics = true;       % 启用详细诊断
    base.debug.save_checkpoints = false;        % 保存检查点
    base.debug.checkpoint_interval = 50;        % 检查点间隔
    
    % ========== 预设配置 ==========
    switch lower(config_name)
        case 'default'
            params = base;
            params.gradient.chain_mode = 'legacy';
            params.opt.acceptance_tol = 0.01;      % default放宽next-guard：允许最多1%相对上浮
            params.opt.current_state_tol = 0.01;   % default放宽current-guard：允许最多1%相对上浮
            params.opt.reinit_current_tol = 0.01;  % default放宽reinit current-guard：允许最多1%相对上浮
            
        case 'fast'
            % 快速模式（牺牲精度换速度）
            params = base;
            params.gradient.chain_mode = 'legacy';
            params.opt.max_iter = 50;
            params.opt.delta_theta_max_deg = 0.8;   % 新灵敏度公式下，0.8度在柔度收益与稳定性之间最好
            params.opt.max_backtrack = 3;          % fast模式适度放宽回溯，减少坏步通过
            params.opt.acceptance_tol = 0.01;      % 放宽next-guard：允许最多1%相对上浮
            params.opt.current_state_tol = 0.01;   % 放宽current-guard：允许最多1%相对上浮
            params.opt.current_guard_start_iter = 1; % 从首轮起启用当前步门禁，阻止中前期坏步通过
            params.opt.reinit_current_tol = 0.01;  % 放宽reinit current-guard：允许最多1%相对上浮
            params.opt.reinit_guard_start_iter = 1;
            params.opt.best_state_patience = 12;   % fast模式提前收敛，避免后段回升
            params.velocity.scale_quantile = 95;
            params.velocity.clip_abs = inf;
            params.velocity.bias_beta = 0.10;
            params.levelset.transition_iter = 50;  % 修复：与max_iter保持一致
            params.debug.enable_plots = false;
            params.debug.enable_diagnostics = false;
            params.debug.log_interval = 20;
            params.levelset.reinit_freq_early = 10;
            params.levelset.reinit_freq_late = 15;
            
        case 'precise'
            % 精确模式（更严格的收敛）
            params = base;
            params.gradient.chain_mode = 'exact';
            params.opt.max_iter = 200;
            params.opt.tol = 1e-6;
            params.opt.best_state_patience = 40;
            params.smooth.iterations = 3;
            params.debug.save_checkpoints = true;
            params.debug.checkpoint_interval = 25;
            params.levelset.reinit_freq_early = 3;
            params.levelset.reinit_threshold = 0.5;
            
        case 'debug'
            % 调试模式（最大化诊断信息）
            params = base;
            params.gradient.chain_mode = 'shadow';
            params.opt.max_iter = 20;
            params.opt.best_state_patience = 8;
            params.opt.acceptance_tol = 0.01;
            params.opt.current_state_tol = 0.01;
            params.opt.reinit_current_tol = 0.01;
            params.levelset.transition_iter = 20;  % 修复：与max_iter保持一致
            params.debug.log_level = 'DEBUG';
            params.debug.log_interval = 1;
            params.debug.enable_diagnostics = true;
            params.debug.verbose = true;
            
        otherwise
            error('未知配置: %s\n支持的配置: default, fast, precise, debug', config_name);
    end
    
    % 自动计算派生参数
    params.grid.dx = params.grid.Lx / params.grid.nelx;
    params.grid.dy = params.grid.Ly / params.grid.nely;
    params.grid.h = min(params.grid.dx, params.grid.dy);
    params.opt.delta_theta_max = params.opt.delta_theta_max_deg * pi/180;

    % 添加配置名称用于日志
    params.config_name = config_name;

    % 可选：环境变量覆盖（仅用于调参试验，不改变默认行为）
    params = apply_env_overrides(params);
end

function params = apply_env_overrides(params)
    % 通用调试覆盖：仅改变运行时诊断/可视化，不改变默认配置文件内容
    v = get_env_numeric('FIBER_FORCE_ENABLE_PLOTS');
    if ~isnan(v)
        params.debug.enable_plots = (v ~= 0);
    end

    v = get_env_numeric('FIBER_FORCE_ENABLE_DIAGNOSTICS');
    if ~isnan(v)
        params.debug.enable_diagnostics = (v ~= 0);
    end

    % 仅对fast配置启用其余调参覆盖
    if ~isfield(params, 'config_name') || ~strcmpi(params.config_name, 'fast')
        return;
    end

    v = get_env_numeric('FIBER_FAST_MAX_BACKTRACK');
    if ~isnan(v)
        params.opt.max_backtrack = round(v);
    end

    v = get_env_numeric('FIBER_FAST_CURRENT_GUARD_START');
    if ~isnan(v)
        params.opt.current_guard_start_iter = round(v);
    end

    v = get_env_numeric('FIBER_FAST_REINIT_GUARD_START');
    if ~isnan(v)
        params.opt.reinit_guard_start_iter = round(v);
    end

    v = get_env_numeric('FIBER_FAST_CURRENT_TOL');
    if ~isnan(v)
        params.opt.current_state_tol = v;
    end

    v = get_env_numeric('FIBER_FAST_REINIT_CURRENT_TOL');
    if ~isnan(v)
        params.opt.reinit_current_tol = v;
    end

    v = get_env_numeric('FIBER_FAST_ENABLE_CURRENT_GUARD');
    if ~isnan(v)
        params.opt.enable_current_state_guard = (v ~= 0);
    end

    v = get_env_numeric('FIBER_FAST_ENABLE_REINIT_CURRENT_GUARD');
    if ~isnan(v)
        params.opt.enable_reinit_current_guard = (v ~= 0);
    end

    v = get_env_numeric('FIBER_FAST_PATIENCE');
    if ~isnan(v)
        params.opt.best_state_patience = round(v);
    end

    v = get_env_numeric('FIBER_FAST_ACCEPTANCE_TOL');
    if ~isnan(v)
        params.opt.acceptance_tol = v;
    end

    v = get_env_numeric('FIBER_FAST_DELTA_THETA_MAX_DEG');
    if ~isnan(v)
        params.opt.delta_theta_max_deg = v;
        params.opt.delta_theta_max = v * pi / 180;
    end

    raw_mode = getenv('FIBER_HJ_RHS_MODE');
    if ~isempty(raw_mode)
        params.levelset.hj_rhs_mode = string(raw_mode);
    end

    raw_mode = getenv('FIBER_GRADIENT_CHAIN_MODE');
    if ~isempty(raw_mode)
        params.gradient.chain_mode = string(raw_mode);
    end
end

function v = get_env_numeric(name)
    raw = getenv(name);
    if isempty(raw)
        v = NaN;
        return;
    end
    v = str2double(raw);
    if ~isfinite(v)
        v = NaN;
    end
end
