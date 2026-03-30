function [material_mask, mask_info] = clean_material_mask(struc, min_area, morph_radius)
    % 在水平集初始化前对拓扑掩膜进行清洗和去噪

    if nargin < 2 || isempty(min_area)
        min_area = 10;
    end
    if nargin < 3 || isempty(morph_radius)
        morph_radius = 1;
    end

    material_mask = struc > 0.5;
    material_mask = bwareaopen(material_mask, max(1, round(min_area)));
    material_mask = imfill(material_mask, 'holes');

    if morph_radius > 0
        se = strel('disk', morph_radius);
        material_mask = imopen(material_mask, se);
    end

    material_mask = logical(material_mask);

    cc = bwconncomp(material_mask);
    areas = cellfun(@numel, cc.PixelIdxList);
    mask_info = struct();
    mask_info.num_components = cc.NumObjects;
    mask_info.total_area = sum(areas);
    mask_info.min_area = min_area;
    mask_info.morph_radius = morph_radius;
    if mask_info.num_components > 0
        [~, max_idx] = max(areas);
        mask_info.largest_component_area = areas(max_idx);
    else
        mask_info.largest_component_area = 0;
    end
end

