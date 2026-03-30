function zero_mask = compute_zero_mask_from_lsf(lsf, bandwidth)
    % 从当前lsf提取“细零集”掩膜，供FMM重初始化使用
    % 目标：保留主路径位置，避免把整条宽窄带都当作零种子

    if nargin < 2 || isempty(bandwidth)
        bandwidth = 2;
    end
    bandwidth = max(bandwidth, eps);

    % 1) 首选符号变号边（更接近phi=0）
    cross_mask = false(size(lsf));
    lr_cross = lsf(:, 1:end-1) .* lsf(:, 2:end) <= 0;
    ud_cross = lsf(1:end-1, :) .* lsf(2:end, :) <= 0;

    cross_mask(:, 1:end-1) = cross_mask(:, 1:end-1) | lr_cross;
    cross_mask(:, 2:end) = cross_mask(:, 2:end) | lr_cross;
    cross_mask(1:end-1, :) = cross_mask(1:end-1, :) | ud_cross;
    cross_mask(2:end, :) = cross_mask(2:end, :) | ud_cross;

    % 2) 用更窄的近零带补点，避免断裂
    tol_zero = 0.25 * bandwidth;
    near_zero = abs(lsf) <= tol_zero;
    narrow = abs(lsf) <= bandwidth;

    zero_mask = (cross_mask | near_zero) & narrow;

    % 3) 细化为近似单像素零集，不再扩张
    if any(zero_mask(:))
        zero_mask = bwmorph(zero_mask, 'clean');
        zero_mask = bwmorph(zero_mask, 'thin', Inf);
    end

    % 4) 兜底：逐步放宽，但仍保持为窄零集
    if ~any(zero_mask(:))
        zero_mask = abs(lsf) <= 0.5 * bandwidth;
    end
    if ~any(zero_mask(:))
        [~, idx_min] = min(abs(lsf(:)));
        zero_mask = false(size(lsf));
        zero_mask(idx_min) = true;
    end

    zero_mask = logical(zero_mask);
    zero_mask(1,:) = zero_mask(2,:);
    zero_mask(end,:) = zero_mask(end-1,:);
    zero_mask(:,1) = zero_mask(:,2);
    zero_mask(:,end) = zero_mask(:,end-1);
end

