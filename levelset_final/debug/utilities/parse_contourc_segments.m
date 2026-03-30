function segments = parse_contourc_segments(C)
%PARSE_CONTOURC_SEGMENTS Convert contourc output to a cell array of Nx2 paths.

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

        level_npts = C(:, idx);
        npts = round(level_npts(2));
        j0 = idx + 1;
        j1 = idx + npts;
        if npts < 2 || j1 > ncol
            break;
        end

        pts = C(:, j0:j1)';
        finite_rows = all(isfinite(pts), 2);
        pts = pts(finite_rows, :);
        if size(pts, 1) >= 2
            segments{end+1} = pts; %#ok<AGROW>
        end
        idx = j1 + 1;
    end
end
