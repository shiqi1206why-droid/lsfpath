function zero_mask = compute_zero_mask_from_lsf(lsf, bandwidth)
    % 通过细化窄带构造零水平集掩膜

    if nargin < 2 || isempty(bandwidth)
        bandwidth = 2;
    end
    narrow = abs(lsf) <= bandwidth;
    zero_mask = bwmorph(narrow, 'thin', Inf);
    zero_mask = imdilate(zero_mask, strel('disk', 1));
    zero_mask = logical(zero_mask);
    if nnz(zero_mask) == 0
        zero_mask = narrow;
    end
    zero_mask(1,:) = zero_mask(2,:);
    zero_mask(end,:) = zero_mask(end-1,:);
    zero_mask(:,1) = zero_mask(:,2);
    zero_mask(:,end) = zero_mask(:,end-1);
end

