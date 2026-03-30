function material_mask_full = expand_material_mask_to_full(material_mask_core)
%EXPAND_MATERIAL_MASK_TO_FULL Expand a core material mask to the LSF grid.

    if ndims(material_mask_core) ~= 2
        error('material_mask_core must be a 2-D array.');
    end

    material_mask_core = logical(material_mask_core);
    material_mask_full = false(size(material_mask_core, 1) + 2, size(material_mask_core, 2) + 2);
    material_mask_full(2:end-1, 2:end-1) = material_mask_core;
    material_mask_full(1, :) = material_mask_full(2, :);
    material_mask_full(end, :) = material_mask_full(end-1, :);
    material_mask_full(:, 1) = material_mask_full(:, 2);
    material_mask_full(:, end) = material_mask_full(:, end-1);
end
