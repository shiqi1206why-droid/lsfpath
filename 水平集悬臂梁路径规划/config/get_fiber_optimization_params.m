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
    base.opt.max_iter = 100;
    base.opt.tol = 1e-5;
    base.opt.alpha = 0.5;
    base.opt.dt = 0.05;
    base.opt.delta_theta_max_deg = 3;
    base.opt.fidelity_weight = 0.02;
    base.opt.normalize_sensitivity = false;      % 控制灵敏度归一化
    base.opt.enable_curvature = true;           % 是否启用曲率正则
    base.opt.normalize_sensitivity = false;      % 控制是否对伴随灵敏度归一化
    
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
    base.levelset.reinit_domain = 'full';       % 重初始化域：'full'=全域传播（默认）,'masked'=限制到材料域（试验）
    
    % 投影参数
    base.projection.enable = true;
    base.projection.omega_early = 0.7;          % 前期投影强度
    base.projection.omega_late = 0.5;           % 后期投影强度
    base.projection.band_factor_early = 1.5;    % 前期投影带宽
    base.projection.band_factor_late = 1.0;     % 后期投影带宽
    
    % 平滑参数
    base.smooth.eta = 0.10;                     % 角度平滑系数
    base.smooth.iterations = 2;                 % 平滑迭代次数
    
    % 载荷参数
    base.load.F_mag = -1;                       % 载荷大小 (N)
    
    % 初始化参数
    base.init.morph_radius = 1;                 % 形态学半径
    base.init.min_component_size = 10;          % 最小连通区域大小
    
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
            
        case 'fast'
            % 快速模式（牺牲精度换速度）
            params = base;
            params.opt.max_iter = 50;
            params.levelset.transition_iter = 50;  % 修复：与max_iter保持一致
            params.debug.enable_plots = false;
            params.debug.enable_diagnostics = false;
            params.debug.log_interval = 20;
            params.levelset.reinit_freq_early = 10;
            params.levelset.reinit_freq_late = 15;
            
        case 'precise'
            % 精确模式（更严格的收敛）
            params = base;
            params.opt.max_iter = 200;
            params.opt.tol = 1e-6;
            params.smooth.iterations = 3;
            params.debug.save_checkpoints = true;
            params.debug.checkpoint_interval = 25;
            params.levelset.reinit_freq_early = 3;
            params.levelset.reinit_threshold = 0.5;
            
        case 'debug'
            % 调试模式（最大化诊断信息）
            params = base;
            params.opt.max_iter = 20;
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
end

