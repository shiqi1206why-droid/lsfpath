clear; clc; close all;
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))), '-begin');

dx = 0.1;
dy = 0.1;
mask = false(14, 18);
mask(4:11, 6:15) = true;

geometry = reconstruct_material_boundary_subpixel(mask, dx, dy);
assert(geometry.success, '子像素边界重建失败。');
assert(geometry.segment_count >= 1, '至少应提取出一条边界折线。');

sign_reference = ones(size(mask));
sign_reference(mask) = -1;
[phi_boundary, dist_info] = build_signed_distance_from_segments( ...
    geometry.x_centers, geometry.y_centers, geometry.segments, sign_reference, true(size(mask)), 0);

assert(dist_info.success, '无法根据重建边界恢复符号距离场。');
assert(all(phi_boundary(mask) < 0), '材料域内的符号距离应为负。');
assert(all(phi_boundary(~mask) > 0), '材料域外的符号距离应为正。');
assert(all(isfinite(phi_boundary(:))), '符号距离场不应包含非有限值。');

fprintf('segment_count=%d\n', geometry.segment_count);
fprintf('inside_max_abs=%.6e\n', max(abs(phi_boundary(mask))));
fprintf('outside_min=%.6e\n', min(phi_boundary(~mask)));
fprintf('PASS: subcell boundary reconstruction test.\n');
