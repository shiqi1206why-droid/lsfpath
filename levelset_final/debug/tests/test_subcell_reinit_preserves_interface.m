clear; clc; close all;
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))), '-begin');

nelx = 28;
nely = 20;
dx = 0.05;
dy = 0.05;

material_mask_core = false(nely, nelx);
material_mask_core(3:18, 4:25) = true;
material_mask_full = expand_material_mask_to_full(material_mask_core);

[x_full, y_full] = get_lsf_grid_coordinates([nely + 2, nelx + 2], dx, dy);
[X, Y] = meshgrid(x_full, y_full);
cx = mean(x_full(6:end-5));
cy = mean(y_full(6:end-5));
r = 0.18;
lsf_ref = sqrt((X - cx).^2 + (Y - cy).^2) - r;
lsf_ref(~material_mask_full) = 0.5;
lsf_ref = impose_neumann(lsf_ref);

lsf_distorted = lsf_ref + 0.12 * min(dx, dy) * sin(4 * pi * X) .* cos(3 * pi * Y);
lsf_distorted(~material_mask_full) = 0.5;
lsf_distorted = impose_neumann(lsf_distorted);

[lsf_reinit, reinit_diag] = fmm_reinitialize(lsf_distorted, dx, dy, [], material_mask_core, ...
    struct('method', 'subcell_signed_distance'));

geom_ref = extract_levelset_geometry(lsf_ref, dx, dy, 0, material_mask_full);
geom_reinit = extract_levelset_geometry(lsf_reinit, dx, dy, 0, material_mask_full);
assert(geom_ref.segment_count > 0 && geom_reinit.segment_count > 0, '零等值线提取失败。');

dist = contour_distance(geom_reinit.segments, geom_ref.segments);
mean_dist = mean(dist);
assert(mean_dist < 0.35 * min(dx, dy), '重初始化前后零等值线位置偏移过大。');
assert(nnz(lsf_reinit(~material_mask_full) <= 0) == 0, '材料域外不应产生 phi<=0。');
assert(strcmpi(reinit_diag.method_used, 'subcell_signed_distance') || reinit_diag.fallback_used, ...
    '重初始化诊断字段异常。');

fprintf('method_used=%s\n', reinit_diag.method_used);
fprintf('mean_contour_shift=%.6e\n', mean_dist);
fprintf('PASS: subcell reinitialization preserves interface test.\n');

function field = impose_neumann(field)
field(1, :) = field(2, :);
field(end, :) = field(end-1, :);
field(:, 1) = field(:, 2);
field(:, end) = field(:, end-1);
end

function distances = contour_distance(source_segments, target_segments)
    distances = [];
    for i = 1:numel(source_segments)
        src = source_segments{i};
        if size(src, 1) < 2
            continue;
        end
        min_dist = inf(size(src, 1), 1);
        for j = 1:numel(target_segments)
            tgt = target_segments{j};
            if size(tgt, 1) < 2
                continue;
            end
            min_dist = min(min_dist, point_to_polyline(src, tgt));
        end
        distances = [distances; min_dist(isfinite(min_dist))]; %#ok<AGROW>
    end
end

function dist = point_to_polyline(points, polyline)
    dist = inf(size(points, 1), 1);
    for s = 1:size(polyline, 1) - 1
        a = polyline(s, :);
        b = polyline(s + 1, :);
        ab = b - a;
        denom = dot(ab, ab);
        if denom <= eps
            proj = repmat(a, size(points, 1), 1);
        else
            t = ((points(:, 1) - a(1)) * ab(1) + (points(:, 2) - a(2)) * ab(2)) / denom;
            t = min(max(t, 0), 1);
            proj = a + t .* ab;
        end
        d = hypot(points(:, 1) - proj(:, 1), points(:, 2) - proj(:, 2));
        dist = min(dist, d);
    end
end
