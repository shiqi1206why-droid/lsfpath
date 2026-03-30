function metrics = compute_raw_path_quality_metrics(lsf, dx, dy, material_mask_core, opts)
%COMPUTE_RAW_PATH_QUALITY_METRICS Adapter around compute_lsf_path_quality.

    if nargin < 5 || isempty(opts)
        opts = struct();
    end
    if nargin < 4 || isempty(material_mask_core)
        material_mask_core = true(size(lsf, 1) - 2, size(lsf, 2) - 2);
    end
    h = min(dx, dy);

    path_opts = struct();
    if isfield(opts, 'resample_ds') && ~isempty(opts.resample_ds)
        path_opts.resample_ds = opts.resample_ds;
    elseif isfield(opts, 'target_ds') && ~isempty(opts.target_ds)
        path_opts.resample_ds = opts.target_ds;
    else
        path_opts.resample_ds = h / 4;
    end

    if isfield(opts, 'zero_bandwidth') && ~isempty(opts.zero_bandwidth)
        path_opts.grad_bandwidth = opts.zero_bandwidth;
    else
        path_opts.grad_bandwidth = 0.75 * h;
    end

    if isfield(opts, 'parallel_spacing') && ~isempty(opts.parallel_spacing)
        target_spacing = opts.parallel_spacing;
    else
        target_spacing = h;
    end
    path_opts.parallel_levels = target_spacing * (1:3);

    base = compute_lsf_path_quality(lsf, dx, dy, material_mask_core, path_opts);

    metrics = base;
    metrics.mean_abs_turn_deg = base.zero_mean_abs_turn_deg;
    metrics.max_abs_turn_deg = base.zero_max_abs_turn_deg;
    metrics.max_abs_kappa = base.zero_max_curvature;
    metrics.parallel_spacing_error = base.parallel_spacing_mean_error;
    metrics.parallel_spacing_error_percent = 100 * base.parallel_spacing_mean_error;
    metrics.parallel_spacing_mean = mean(base.parallel_spacing_measured(isfinite(base.parallel_spacing_measured)), 'omitnan');
    metrics.parallel_spacing_target = target_spacing;
    metrics.parallel_spacing_measured = metrics.parallel_spacing_mean;
    metrics.grad_dev_mean = base.grad_abs_deviation_mean;
    metrics.grad_dev_p95 = base.grad_abs_deviation_p95;
    metrics.grad_dev_max = base.grad_abs_deviation_max;
    metrics.grad_outlier_ratio_0p5_1p5 = base.grad_outlier_ratio_0p5_1p5;
    metrics.near_zero_grad_outlier_ratio = base.near_zero_grad_outlier_ratio;
    metrics.high_grad_outlier_ratio = base.high_grad_outlier_ratio;
    metrics.low_grad_outlier_ratio = base.low_grad_outlier_ratio;
    metrics.high_grad_boundary_overlap_ratio = base.high_grad_boundary_overlap_ratio;
    metrics.near_zero_grad_median = base.near_zero_grad_median;
    metrics.near_zero_grad_p95 = base.near_zero_grad_p95;
    metrics.gradient_scopes = base.gradient_scopes;
    metrics.gradient_primary_scope = base.gradient_primary_scope;
    metrics.zero_segments = base.zero_contour.segments;
    metrics.segment_count = base.zero_segment_count;
    metrics.outside_nonpositive_count = nnz(lsf(~expand_material_mask_to_full(material_mask_core)) <= 0);
end
