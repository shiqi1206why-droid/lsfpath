function mask_full = normalize_mask_to_lsf_grid(mask_in, lsf_size, context_name)
%NORMALIZE_MASK_TO_LSF_GRID Normalize core/full masks to the LSF grid size.

    if nargin < 3 || isempty(context_name)
        context_name = 'mask';
    end

    if isequal(size(mask_in), lsf_size)
        mask_full = logical(mask_in);
        return;
    end

    core_size = [lsf_size(1) - 2, lsf_size(2) - 2];
    if all(size(mask_in) == core_size)
        mask_full = false(lsf_size);
        mask_full(2:end-1, 2:end-1) = logical(mask_in);
        mask_full(1, :) = mask_full(2, :);
        mask_full(end, :) = mask_full(end-1, :);
        mask_full(:, 1) = mask_full(:, 2);
        mask_full(:, end) = mask_full(:, end-1);
        return;
    end

    error('%s尺寸不匹配：期望[%d,%d]或[%d,%d]，实际[%d,%d]。', ...
        context_name, lsf_size(1), lsf_size(2), core_size(1), core_size(2), ...
        size(mask_in, 1), size(mask_in, 2));
end
