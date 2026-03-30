clear; clc; close all;
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))), '-begin');

nelx = 24;
nely = 18;
dx = 0.05;
dy = 0.05;
h = min(dx, dy);

material_mask_core = true(nely, nelx);
material_mask_full = expand_material_mask_to_full(material_mask_core);
[x_full, y_full] = get_lsf_grid_coordinates([nely + 2, nelx + 2], dx, dy);
[X, Y] = meshgrid(x_full, y_full);
lsf = (X - mean(x_full)) + 0.02 * sin(2 * pi * Y / max(y_full));
velocity = ones(size(lsf));

phi_boundary_full = min(cat(3, X - min(x_full), max(x_full) - X, Y - min(y_full), max(y_full) - Y), [], 3);
boundary_guard_band = abs(phi_boundary_full) <= 1.5 * h;
narrow_band = abs(lsf) <= 1.5 * h;
primary_update_mask = narrow_band & material_mask_full & ~boundary_guard_band;
stencil_mask = imdilate(primary_update_mask, strel('square', 5)) & material_mask_full;

[lsf_new, ~] = update_levelset_HJ(lsf, velocity, 0.02, dx, dy, primary_update_mask, ...
    struct('advection_order', 2, 'time_integrator', 'ssprk2', 'freeze_on_incomplete_godunov', true, ...
           'stencil_mask', stencil_mask));

guard_point = find(narrow_band & boundary_guard_band, 1, 'first');
active_point = find(primary_update_mask, 1, 'first');
assert(~isempty(guard_point) && ~isempty(active_point), '测试场未形成足够的护带/活动带样本。');
assert(abs(lsf_new(guard_point) - lsf(guard_point)) < 1e-14, '边界护带内点不应更新。');
assert(abs(lsf_new(active_point) - lsf(active_point)) > 1e-12, '活动带内至少应有一点发生更新。');

fprintf('guard_point_change=%.6e\n', abs(lsf_new(guard_point) - lsf(guard_point)));
fprintf('active_point_change=%.6e\n', abs(lsf_new(active_point) - lsf(active_point)));
fprintf('PASS: boundary guard band blocks updates test.\n');
