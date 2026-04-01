function [velocity_field, stats] = build_velocity_exact(node_sensitivity, lsf, dx, dy, bandwidth, remove_bias, gamma_curv, iter, propagation_mask, velocity_opts)
%BUILD_VELOCITY_EXACT Adapt dC/dphi to HJ normal velocity before shaping.

    if nargin < 5 || isempty(bandwidth)
        bandwidth = 2 * min(dx, dy);
    end
    if nargin < 6 || isempty(remove_bias)
        remove_bias = true;
    end
    if nargin < 7 || isempty(gamma_curv)
        gamma_curv = 0;
    end
    if nargin < 8 || isempty(iter)
        iter = 1;
    end
    if nargin < 9 || isempty(propagation_mask)
        propagation_mask = true(size(lsf));
    elseif isstruct(propagation_mask) && (nargin < 10 || isempty(velocity_opts))
        velocity_opts = propagation_mask;
        propagation_mask = true(size(lsf));
    end
    if nargin < 10 || isempty(velocity_opts)
        velocity_opts = struct();
    end

    if ~isfield(velocity_opts, 'exact_grad_floor') || isempty(velocity_opts.exact_grad_floor)
        velocity_opts.exact_grad_floor = 1e-12;
    end

    if ~isequal(size(node_sensitivity), size(lsf))
        sensitivity_expanded = zeros(size(lsf));
        sensitivity_expanded(2:end-1, 2:end-1) = node_sensitivity;
        sensitivity_expanded(1, :) = sensitivity_expanded(2, :);
        sensitivity_expanded(end, :) = sensitivity_expanded(end-1, :);
        sensitivity_expanded(:, 1) = sensitivity_expanded(:, 2);
        sensitivity_expanded(:, end) = sensitivity_expanded(:, end-1);
        node_sensitivity = sensitivity_expanded;
    end

    propagation_mask = normalize_mask_to_lsf_grid(propagation_mask, size(lsf), 'propagation_mask');
    [grad_y, grad_x] = gradient(lsf, dy, dx);
    grad_magnitude = hypot(grad_x, grad_y);
    grad_floor = velocity_opts.exact_grad_floor;
    grad_magnitude = max(grad_magnitude, grad_floor);

    adapted_sensitivity = -node_sensitivity ./ grad_magnitude;
    adapted_sensitivity(~isfinite(adapted_sensitivity)) = 0;
    adapted_sensitivity(~propagation_mask) = 0;

    [velocity_field, stats] = build_velocity_field( ...
        adapted_sensitivity, lsf, dx, dy, bandwidth, remove_bias, gamma_curv, iter, propagation_mask, velocity_opts);

    grad_values = grad_magnitude(propagation_mask);
    if isempty(grad_values)
        mean_grad_magnitude = NaN;
        max_grad_magnitude = NaN;
    else
        mean_grad_magnitude = mean(grad_values, 'omitnan');
        max_grad_magnitude = max(grad_values);
    end
    stats.exact_velocity_adapter = struct( ...
        'grad_floor', grad_floor, ...
        'mean_grad_magnitude', mean_grad_magnitude, ...
        'max_grad_magnitude', max_grad_magnitude);
end
