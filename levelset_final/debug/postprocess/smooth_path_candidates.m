function out = smooth_path_candidates(raw_xy, material_mask_core, dx, dy, opts)
%SMOOTH_PATH_CANDIDATES Build and select smoothing candidates for one path segment.

    xr = raw_xy(:, 1);
    yr = raw_xy(:, 2);
    closed = hypot(xr(1) - xr(end), yr(1) - yr(end)) < max(dx, dy);
    [x_res, y_res] = resample_polyline_local(xr, yr, opts.target_ds, closed);
    mraw = polyline_metrics_local(x_res, y_res);

    candidates = repmat(struct( ...
        'x', [], 'y', [], ...
        'metrics', struct(), ...
        'meta', struct('method', "", 'param1', NaN, 'param2', NaN)), 0, 1);

    if opts.auto_optimize
        methods = opts.candidate_methods;
        if ischar(methods) || isstring(methods)
            methods = cellstr(methods);
        end
        for mi = 1:numel(methods)
            method = char(methods{mi});
            switch lower(method)
                case 'moving_average'
                    for win = opts.candidate_windows
                        if mod(win, 2) == 0
                            win = win + 1;
                        end
                        for niter = opts.candidate_iters
                            [xs, ys] = smooth_polyline_local(x_res, y_res, win, niter, closed);
                            cand = build_candidate(xs, ys, x_res, y_res, material_mask_core, dx, dy, ...
                                "moving_average", win, niter);
                            candidates(end+1) = cand; %#ok<AGROW>
                        end
                    end
                case 'chaikin'
                    for niter = opts.chaikin_iters
                        [xs, ys] = chaikin_smooth_polyline_local(x_res, y_res, niter, closed);
                        cand = build_candidate(xs, ys, x_res, y_res, material_mask_core, dx, dy, ...
                            "chaikin", niter, 0);
                        candidates(end+1) = cand; %#ok<AGROW>
                    end
                otherwise
                    % 未知方法直接跳过，避免影响主流程
            end
        end
    else
        [xs, ys] = smooth_polyline_local(x_res, y_res, opts.smooth_window, opts.smooth_iters, closed);
        cand = build_candidate(xs, ys, x_res, y_res, material_mask_core, dx, dy, ...
            "moving_average", opts.smooth_window, opts.smooth_iters);
        candidates(end+1) = cand; %#ok<AGROW>
    end

    if isempty(candidates)
        [xs, ys] = smooth_polyline_local(x_res, y_res, opts.smooth_window, opts.smooth_iters, closed);
        cand = build_candidate(xs, ys, x_res, y_res, material_mask_core, dx, dy, ...
            "moving_average", opts.smooth_window, opts.smooth_iters);
        candidates(end+1) = cand; %#ok<AGROW>
    end

    best = select_best_smooth_path(candidates, opts.max_deviation_limit);

    out = struct();
    out.raw_x = x_res;
    out.raw_y = y_res;
    out.smooth_x = best.x;
    out.smooth_y = best.y;
    out.raw_metrics = mraw;
    out.smooth_metrics = best.metrics;
    out.meta = best.meta;
end

function cand = build_candidate(xs, ys, x_ref, y_ref, material_mask_core, dx, dy, method, p1, p2)
    m = polyline_metrics_local(xs, ys);
    m.max_deviation_from_raw = max_min_distance_local(xs, ys, x_ref, y_ref);
    inside = is_points_in_material_local(xs, ys, material_mask_core, dx, dy);
    if any(~inside)
        m.max_abs_kappa = inf;
        m.mean_abs_turn_deg = inf;
        m.max_deviation_from_raw = inf;
    end

    cand = struct();
    cand.x = xs;
    cand.y = ys;
    cand.metrics = m;
    cand.meta = struct('method', method, 'param1', p1, 'param2', p2);
end

function inside = is_points_in_material_local(x, y, material_mask_core, dx, dy)
    nely = size(material_mask_core, 1);
    nelx = size(material_mask_core, 2);
    col = round(x / dx + 1.5);
    row = round(y / dy + 1.5);
    col = max(1, min(nelx, col));
    row = max(1, min(nely, row));
    idx = sub2ind([nely, nelx], row, col);
    inside = material_mask_core(idx);
end

function [xo, yo] = resample_polyline_local(x, y, target_ds, closed)
    if numel(x) < 2
        xo = x;
        yo = y;
        return;
    end
    dseg = hypot(diff(x), diff(y));
    s = [0; cumsum(dseg(:))];
    total = s(end);
    if total <= eps
        xo = x;
        yo = y;
        return;
    end
    npts = max(2, round(total / max(target_ds, eps)) + 1);
    s_new = linspace(0, total, npts).';
    xo = interp1(s, x(:), s_new, 'linear');
    yo = interp1(s, y(:), s_new, 'linear');
    if closed
        xo(end) = xo(1);
        yo(end) = yo(1);
    end
end

function [xs, ys] = smooth_polyline_local(x, y, win, niter, closed)
    xs = x(:);
    ys = y(:);
    win = max(3, win);
    if mod(win, 2) == 0
        win = win + 1;
    end

    for it = 1:niter
        if closed
            xs = periodic_moving_average_local(xs, win);
            ys = periodic_moving_average_local(ys, win);
        else
            xs = smoothdata(xs, 'movmean', win);
            ys = smoothdata(ys, 'movmean', win);
            xs(1) = x(1);
            ys(1) = y(1);
            xs(end) = x(end);
            ys(end) = y(end);
        end
    end
end

function [xs, ys] = chaikin_smooth_polyline_local(x, y, niter, closed)
    xs = x(:).';
    ys = y(:).';
    if numel(xs) < 3 || niter <= 0
        xs = xs(:);
        ys = ys(:);
        return;
    end

    for it = 1:niter
        if closed
            x_in = xs;
            y_in = ys;
            if ~(xs(1) == xs(end) && ys(1) == ys(end))
                x_in(end+1) = x_in(1); %#ok<AGROW>
                y_in(end+1) = y_in(1); %#ok<AGROW>
            end
            qx = 0.75 * x_in(1:end-1) + 0.25 * x_in(2:end);
            qy = 0.75 * y_in(1:end-1) + 0.25 * y_in(2:end);
            rx = 0.25 * x_in(1:end-1) + 0.75 * x_in(2:end);
            ry = 0.25 * y_in(1:end-1) + 0.75 * y_in(2:end);
            xs = reshape([qx; rx], 1, []);
            ys = reshape([qy; ry], 1, []);
            xs(end+1) = xs(1); %#ok<AGROW>
            ys(end+1) = ys(1); %#ok<AGROW>
        else
            qx = 0.75 * xs(1:end-1) + 0.25 * xs(2:end);
            qy = 0.75 * ys(1:end-1) + 0.25 * ys(2:end);
            rx = 0.25 * xs(1:end-1) + 0.75 * xs(2:end);
            ry = 0.25 * ys(1:end-1) + 0.75 * ys(2:end);
            x_mid = reshape([qx; rx], 1, []);
            y_mid = reshape([qy; ry], 1, []);
            xs = [xs(1), x_mid, xs(end)];
            ys = [ys(1), y_mid, ys(end)];
        end
    end

    xs = xs(:);
    ys = ys(:);
end

function ys = periodic_moving_average_local(x, win)
    n = numel(x);
    r = floor(win / 2);
    idx = (1:n).';
    ys = zeros(size(x));
    for ii = 1:n
        ids = mod((idx(ii)-r):(idx(ii)+r)-1, n) + 1;
        ys(ii) = mean(x(ids));
    end
end

function m = polyline_metrics_local(x, y)
    dx = diff(x(:));
    dy = diff(y(:));
    ds = hypot(dx, dy);
    m.length = sum(ds);
    if numel(ds) < 2
        m.mean_abs_turn_deg = 0;
        m.max_abs_turn_deg = 0;
        m.max_abs_kappa = 0;
        m.min_turn_radius = inf;
        m.max_deviation_from_raw = 0;
        return;
    end

    tx = dx ./ max(ds, eps);
    ty = dy ./ max(ds, eps);
    heading = atan2(ty, tx);
    dtheta = diff(heading);
    dtheta = wrap_to_pi_local(dtheta);
    m.mean_abs_turn_deg = mean(abs(dtheta)) * 180 / pi;
    m.max_abs_turn_deg = max(abs(dtheta)) * 180 / pi;

    ds_mid = 0.5 * (ds(1:end-1) + ds(2:end));
    kappa = dtheta ./ max(ds_mid, eps);
    m.max_abs_kappa = max(abs(kappa));
    if m.max_abs_kappa > 0
        m.min_turn_radius = 1 / m.max_abs_kappa;
    else
        m.min_turn_radius = inf;
    end
    m.max_deviation_from_raw = NaN;
end

function w = wrap_to_pi_local(a)
    w = mod(a + pi, 2*pi) - pi;
end

function dmax = max_min_distance_local(xa, ya, xb, yb)
    if isempty(xa) || isempty(xb)
        dmax = inf;
        return;
    end
    d2 = pdist2([xa(:), ya(:)], [xb(:), yb(:)]);
    dmax = max(min(d2, [], 2));
end
