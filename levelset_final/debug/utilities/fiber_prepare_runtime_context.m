function runtime_ctx = fiber_prepare_runtime_context(config_name, bootstrap_root)
%FIBER_PREPARE_RUNTIME_CONTEXT Prepare runtime paths, params, and derived options.

    if nargin < 2 || isempty(bootstrap_root)
        bootstrap_root = fileparts(mfilename('fullpath'));
    end

    addpath(fullfile(bootstrap_root, 'utilities'), '-begin');
    project_root = get_project_root(bootstrap_root);
    cleanup = ensure_project_on_path(project_root);
    paths = build_project_paths(project_root);

    fprintf('=== 纤维路径优化（配置: %s） ===\n', config_name);
    params = get_fiber_optimization_params(config_name);
    validate_params(params);
    params.runtime = struct('project_root', project_root, 'paths', paths);

    runtime_ctx = struct();
    runtime_ctx.config_name = config_name;
    runtime_ctx.project_root = project_root;
    runtime_ctx.paths = paths;
    runtime_ctx.params = params;
    runtime_ctx.path_cleanup = cleanup;

    runtime_ctx.nelx = params.grid.nelx;
    runtime_ctx.nely = params.grid.nely;
    runtime_ctx.Lx = params.grid.Lx;
    runtime_ctx.Ly = params.grid.Ly;
    runtime_ctx.dx = params.grid.dx;
    runtime_ctx.dy = params.grid.dy;

    runtime_ctx.E_L = params.material.E_L;
    runtime_ctx.E_T = params.material.E_T;
    runtime_ctx.nu_LT = params.material.nu_LT;
    runtime_ctx.nu_TL = params.material.nu_TL;
    runtime_ctx.G_LT = params.material.G_LT;
    runtime_ctx.G_LW = params.material.G_LW;
    runtime_ctx.G_TW = params.material.G_TW;
    runtime_ctx.thickness = params.material.thickness;

    runtime_ctx.max_iter = params.opt.max_iter;
    runtime_ctx.tol = params.opt.tol;
    runtime_ctx.alpha = params.opt.alpha;
    runtime_ctx.dt = params.opt.dt;
    runtime_ctx.delta_theta_max = params.opt.delta_theta_max;
    runtime_ctx.fidelity_weight = params.opt.fidelity_weight;
    runtime_ctx.enable_step_acceptance = params.opt.enable_step_acceptance;
    runtime_ctx.acceptance_tol = params.opt.acceptance_tol;
    runtime_ctx.enable_current_state_guard = logical(params.opt.enable_current_state_guard);
    runtime_ctx.current_state_tol = params.opt.current_state_tol;
    runtime_ctx.current_guard_start_iter = params.opt.current_guard_start_iter;
    runtime_ctx.enable_reinit_current_guard = logical(params.opt.enable_reinit_current_guard);
    runtime_ctx.reinit_current_tol = params.opt.reinit_current_tol;
    runtime_ctx.reinit_guard_start_iter = params.opt.reinit_guard_start_iter;
    runtime_ctx.backtrack_factor = params.opt.backtrack_factor;
    runtime_ctx.max_backtrack = params.opt.max_backtrack;
    runtime_ctx.min_backtrack_dt = params.opt.min_backtrack_dt;
    runtime_ctx.enable_best_state_guard = logical(params.opt.enable_best_state_guard);
    runtime_ctx.best_state_rel_tol = params.opt.best_state_rel_tol;
    runtime_ctx.best_state_patience = params.opt.best_state_patience;
    runtime_ctx.theta_only_fuse_limit = params.opt.theta_only_fuse_limit;
    runtime_ctx.velocity_opts = params.velocity;
    runtime_ctx.gradient_opts = params.gradient;
    runtime_ctx.manufacturing_opts = params.manufacturing;

    runtime_ctx.delta_phi = params.levelset.delta_phi_factor * params.grid.h;
    runtime_ctx.init_smooth_opts = struct( ...
        'morph_radius', params.init.morph_radius, ...
        'boundary_reconstruction', params.init.boundary_reconstruction, ...
        'reinit_method', params.levelset.reinit_method, ...
        'zero_geometry_min_points', params.levelset.zero_geometry_min_points, ...
        'zero_geometry_min_length', params.levelset.zero_geometry_min_length);
    runtime_ctx.hj_update_opts = struct( ...
        'advection_order', params.levelset.advection_order, ...
        'time_integrator', params.levelset.time_integrator, ...
        'fallback_first_order', logical(params.levelset.fallback_first_order), ...
        'freeze_on_incomplete_godunov', logical(params.levelset.freeze_on_incomplete_godunov), ...
        'stencil_buffer_cells', params.levelset.stencil_buffer_cells, ...
        'eno_smoothness_factor', params.levelset.eno_smoothness_factor, ...
        'rhs_mode', char(params.levelset.hj_rhs_mode));
    runtime_ctx.reinit_opts = struct( ...
        'method', params.levelset.reinit_method, ...
        'zero_geometry_min_points', params.levelset.zero_geometry_min_points, ...
        'zero_geometry_min_length', params.levelset.zero_geometry_min_length, ...
        'local_shell_mask', [], ...
        'preserve_outside_shell', false);
    runtime_ctx.path_quality_opts = struct('target_ds', params.grid.h / 4, ...
        'resample_ds', params.grid.h / 4, ...
        'parallel_spacing', params.grid.h, ...
        'zero_bandwidth', 0.5 * params.grid.h, ...
        'grad_outlier_bounds', [0.5, 1.5], ...
        'boundary_overlap_bandwidth', 2.0 * params.grid.h);

    runtime_ctx.F_mag = params.load.F_mag;
end
