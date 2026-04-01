function diagnostics = compute_gradient_chain_diagnostics(legacy_grad, exact_opt_grad, exact_full_grad, active_mask, material_mask_full, transition_cache, topk)
%COMPUTE_GRADIENT_CHAIN_DIAGNOSTICS Compare legacy and exact-chain gradients.

    if nargin < 7 || isempty(topk)
        topk = 32;
    end

    active_mask = logical(active_mask);
    material_mask_full = logical(material_mask_full);
    compare_mask = active_mask & material_mask_full;

    legacy_vec = legacy_grad(compare_mask);
    exact_vec = exact_opt_grad(compare_mask);

    legacy_vec(~isfinite(legacy_vec)) = 0;
    exact_vec(~isfinite(exact_vec)) = 0;

    legacy_norm = norm(legacy_vec);
    exact_norm = norm(exact_vec);
    if legacy_norm > 0 && exact_norm > 0
        cosine_similarity = dot(legacy_vec, exact_vec) / (legacy_norm * exact_norm);
    else
        cosine_similarity = NaN;
    end

    support_legacy = abs(legacy_vec) > 1e-12;
    support_exact = abs(exact_vec) > 1e-12;
    support_union = support_legacy | support_exact;
    if any(support_union)
        band_support_overlap = nnz(support_legacy & support_exact) / nnz(support_union);
    else
        band_support_overlap = 1;
    end

    combined_mag = max(abs(legacy_vec), abs(exact_vec));
    [~, order] = sort(combined_mag, 'descend');
    topk = min(topk, numel(order));
    if topk > 0
        idx = order(1:topk);
        sign_agree = sign_match_fraction(legacy_vec(idx), exact_vec(idx));
    else
        sign_agree = NaN;
    end

    diagnostics = struct();
    diagnostics.cosine_similarity = cosine_similarity;
    diagnostics.norm_ratio = exact_norm / max(legacy_norm, 1e-12);
    diagnostics.topk_sign_agreement = sign_agree;
    diagnostics.saturation_ratio = safe_mask_fraction(transition_cache.limiter_saturated_mask, transition_cache.material_mask_core);
    diagnostics.degenerate_ratio = safe_mask_fraction(transition_cache.degenerate_grad_mask, transition_cache.material_mask_core);
    diagnostics.active_band_coverage = safe_fraction(nnz(compare_mask), nnz(material_mask_full));
    diagnostics.zero_gradient_due_to_limiter_ratio = diagnostics.saturation_ratio;
    diagnostics.exact_gradient_active_nonzero_ratio = safe_fraction(nnz(support_exact), numel(exact_vec));
    diagnostics.exact_vs_legacy_band_support_overlap = band_support_overlap;
    diagnostics.full_vs_opt_support_overlap = compute_support_overlap( ...
        abs(exact_opt_grad(material_mask_full)) > 1e-12, ...
        abs(exact_full_grad(material_mask_full)) > 1e-12);
end

function overlap = compute_support_overlap(a, b)
    union_mask = a | b;
    if any(union_mask)
        overlap = nnz(a & b) / nnz(union_mask);
    else
        overlap = 1;
    end
end

function frac = sign_match_fraction(a, b)
    a = a(:);
    b = b(:);
    both_zero = (abs(a) <= 1e-12) & (abs(b) <= 1e-12);
    same_sign = sign(a) == sign(b);
    frac = mean(both_zero | same_sign);
end

function ratio = safe_mask_fraction(mask, base_mask)
    mask = logical(mask);
    base_mask = logical(base_mask);
    ratio = safe_fraction(nnz(mask & base_mask), nnz(base_mask));
end

function ratio = safe_fraction(num, den)
    if den <= 0
        ratio = 0;
    else
        ratio = num / den;
    end
end
