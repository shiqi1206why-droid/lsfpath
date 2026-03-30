function mask_out = dilate_binary_mask(mask_in, radius_cells)
%DILATE_BINARY_MASK Dilate a binary mask with a square structuring element.

    mask_in = logical(mask_in);
    if radius_cells <= 0
        mask_out = mask_in;
        return;
    end
    window_size = 2 * radius_cells + 1;
    mask_out = imdilate(mask_in, strel('square', window_size));
end
