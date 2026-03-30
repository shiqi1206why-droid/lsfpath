function contour_entries = extract_path_contours(lsf, levels, x_coords, y_coords, material_mask_core, dx, dy, opts)
%EXTRACT_PATH_CONTOURS Extract and material-clip contour segments from LSF.

    contour_entries = {};
    for li = 1:numel(levels)
        lv = levels(li);
        C = contourc(x_coords, y_coords, lsf, [lv, lv]);
        segments = parse_contourc_segments_local(C);
        for si = 1:numel(segments)
            clipped_segments = clip_segment_to_material_local(segments{si}, material_mask_core, dx, dy, opts);
            for cj = 1:numel(clipped_segments)
                entry = struct();
                entry.level = lv;
                entry.xy = clipped_segments{cj};
                contour_entries{end+1} = entry; %#ok<AGROW>
            end
        end
    end
end

function segments = parse_contourc_segments_local(C)
    segments = {};
    if isempty(C)
        return;
    end

    idx = 1;
    ncol = size(C, 2);
    while idx <= ncol
        if idx + 1 > ncol
            break;
        end
        npts = C(2, idx);
        j0 = idx + 1;
        j1 = min(idx + npts, ncol);
        if j0 <= j1
            xy = C(:, j0:j1).';
            if size(xy, 1) >= 2
                segments{end+1} = xy; %#ok<AGROW>
            end
        end
        idx = idx + npts + 1;
    end
end

function clipped_segments = clip_segment_to_material_local(raw_xy, material_mask_core, dx, dy, opts)
    clipped_segments = {};
    if isempty(raw_xy)
        return;
    end
    inside = is_points_in_material_local(raw_xy(:, 1), raw_xy(:, 2), material_mask_core, dx, dy);
    if all(inside)
        clipped_segments = {raw_xy};
        return;
    end
    if ~any(inside)
        return;
    end

    idx_inside = find(inside);
    breaks = [1; find(diff(idx_inside) > 1) + 1; numel(idx_inside) + 1];
    for bi = 1:numel(breaks)-1
        ids = idx_inside(breaks(bi):breaks(bi+1)-1);
        seg = raw_xy(ids, :);
        if size(seg, 1) >= max(2, floor(opts.min_points / 2))
            clipped_segments{end+1} = seg; %#ok<AGROW>
        end
    end
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
