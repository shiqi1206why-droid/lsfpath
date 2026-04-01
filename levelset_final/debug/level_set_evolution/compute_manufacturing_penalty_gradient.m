function [penalty_gradient, diagnostics] = compute_manufacturing_penalty_gradient(lsf, dx, dy, material_mask_full, lsf_target_global, active_mask, manufacturing_opts)
%COMPUTE_MANUFACTURING_PENALTY_GRADIENT Add optimization-side manufacturing penalties.

    if nargin < 7 || isempty(manufacturing_opts) || ~isfield(manufacturing_opts, 'enable') || ...
            ~logical(manufacturing_opts.enable)
        penalty_gradient = zeros(size(lsf));
        diagnostics = init_diagnostics();
        return;
    end

    band_mask = logical(active_mask) & logical(material_mask_full);
    penalty_gradient = zeros(size(lsf));
    diagnostics = init_diagnostics();
    diagnostics.active_count = nnz(band_mask);

    if diagnostics.active_count == 0
        return;
    end

    grad_norm_weight = manufacturing_opts.grad_norm_weight;
    curvature_weight = manufacturing_opts.curvature_weight;
    gap_overlap_weight = manufacturing_opts.gap_overlap_weight;

    [gy, gx] = gradient(lsf, dy, dx);
    grad_mag = hypot(gx, gy);
    grad_mag_safe = max(grad_mag, 1e-12);

    if grad_norm_weight > 0
        residual = grad_mag - 1;
        flux_x = residual .* (gx ./ grad_mag_safe);
        flux_y = residual .* (gy ./ grad_mag_safe);
        grad_norm_grad = -grad_norm_weight * divergence_like(flux_x, flux_y, dx, dy);
        grad_norm_grad(~band_mask) = 0;
        penalty_gradient = penalty_gradient + grad_norm_grad;
        diagnostics.grad_norm_weight = grad_norm_weight;
        diagnostics.grad_norm_grad_norm = norm(grad_norm_grad(:));
        diagnostics.grad_norm_energy = 0.5 * grad_norm_weight * sum((residual(band_mask)).^2, 'omitnan');
    end

    if curvature_weight > 0
        nx = gx ./ grad_mag_safe;
        ny = gy ./ grad_mag_safe;
        [~, nx_x] = gradient(nx, dy, dx);
        [ny_y, ~] = gradient(ny, dy, dx);
        kappa = nx_x + ny_y;
        kappa_max = 1 / max(manufacturing_opts.curvature_radius_min, 1e-12);
        excess = sign(kappa) .* max(abs(kappa) - kappa_max, 0);
        curvature_grad = -curvature_weight * laplacian_like(excess, dx, dy);
        curvature_grad(~band_mask) = 0;
        penalty_gradient = penalty_gradient + curvature_grad;
        diagnostics.curvature_weight = curvature_weight;
        diagnostics.curvature_grad_norm = norm(curvature_grad(:));
        diagnostics.curvature_energy = 0.5 * curvature_weight * sum((excess(band_mask)).^2, 'omitnan');
    end

    if gap_overlap_weight > 0 && isequal(size(lsf_target_global), size(lsf))
        gap_overlap_grad = gap_overlap_weight * (lsf - lsf_target_global);
        gap_overlap_grad(~band_mask) = 0;
        penalty_gradient = penalty_gradient + gap_overlap_grad;
        diagnostics.gap_overlap_weight = gap_overlap_weight;
        diagnostics.gap_overlap_grad_norm = norm(gap_overlap_grad(:));
        diagnostics.gap_overlap_energy = 0.5 * gap_overlap_weight * sum((lsf(band_mask) - lsf_target_global(band_mask)).^2, 'omitnan');
    end

    penalty_gradient(~isfinite(penalty_gradient)) = 0;
end

function diagnostics = init_diagnostics()
    diagnostics = struct();
    diagnostics.active_count = 0;
    diagnostics.grad_norm_weight = 0;
    diagnostics.curvature_weight = 0;
    diagnostics.gap_overlap_weight = 0;
    diagnostics.grad_norm_grad_norm = 0;
    diagnostics.curvature_grad_norm = 0;
    diagnostics.gap_overlap_grad_norm = 0;
    diagnostics.grad_norm_energy = 0;
    diagnostics.curvature_energy = 0;
    diagnostics.gap_overlap_energy = 0;
end

function div_field = divergence_like(fx, fy, dx, dy)
    [~, fx_x] = gradient(fx, dy, dx);
    [fy_y, ~] = gradient(fy, dy, dx);
    div_field = fx_x + fy_y;
end

function lap = laplacian_like(field, dx, dy)
    [gy, gx] = gradient(field, dy, dx);
    lap = divergence_like(gx, gy, dx, dy);
end
