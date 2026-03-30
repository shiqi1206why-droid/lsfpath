function best = select_best_smooth_path(candidates, max_deviation_limit)
%SELECT_BEST_SMOOTH_PATH Select best smoothing candidate by quality metrics.

    if isempty(candidates)
        error('select_best_smooth_path:EmptyCandidates', '候选平滑路径为空。');
    end

    valid_mask = false(size(candidates));
    for i = 1:numel(candidates)
        m = candidates(i).metrics;
        valid_mask(i) = isfinite(m.max_abs_kappa) && isfinite(m.mean_abs_turn_deg) && ...
            isfinite(m.max_deviation_from_raw) && (m.max_deviation_from_raw <= max_deviation_limit);
    end
    if any(valid_mask)
        pool = candidates(valid_mask);
    else
        pool = candidates;
    end

    best = pool(1);
    for i = 2:numel(pool)
        cand = pool(i);
        if is_better_candidate(cand, best)
            best = cand;
        end
    end
end

function tf = is_better_candidate(cand_new, cand_old)
    m_new = cand_new.metrics;
    m_old = cand_old.metrics;
    if m_new.max_abs_kappa < m_old.max_abs_kappa * 0.98
        tf = true;
        return;
    end
    if m_new.max_abs_kappa > m_old.max_abs_kappa * 1.02
        tf = false;
        return;
    end
    if m_new.mean_abs_turn_deg < m_old.mean_abs_turn_deg * 0.98
        tf = true;
        return;
    end
    if m_new.mean_abs_turn_deg > m_old.mean_abs_turn_deg * 1.02
        tf = false;
        return;
    end
    if m_new.max_deviation_from_raw < m_old.max_deviation_from_raw * 0.98
        tf = true;
        return;
    end
    if m_new.max_deviation_from_raw > m_old.max_deviation_from_raw * 1.02
        tf = false;
        return;
    end
    tf = prefer_method(cand_new.meta.method, cand_old.meta.method);
end

function tf = prefer_method(method_new, method_old)
    order = {'moving_average', 'chaikin'};
    i_new = find(strcmpi(order, char(method_new)), 1, 'first');
    i_old = find(strcmpi(order, char(method_old)), 1, 'first');
    if isempty(i_new)
        i_new = numel(order) + 1;
    end
    if isempty(i_old)
        i_old = numel(order) + 1;
    end
    tf = i_new < i_old;
end
