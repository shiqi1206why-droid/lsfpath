clear; clc; close all;
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
addpath(fullfile(project_root, 'utilities'), '-begin');
project_root = get_project_root(project_root);
cleanup_path = ensure_project_on_path(project_root); %#ok<NASGU>

nelx = 24;
nely = 18;
dx = 0.04;
dy = 0.04;

material_mask_core = false(nely, nelx);
material_mask_core(4:15, 5:21) = true;
material_mask_full = expand_material_mask_to_full(material_mask_core);

[x_full, y_full] = get_lsf_grid_coordinates([nely + 2, nelx + 2], dx, dy);
[X, Y] = meshgrid(x_full, y_full);
x_core = x_full(2:end-1);
lsf = (X - mean(x_core)) + 0.1 * sin(2 * pi * Y / max(y_full));

velocity = ones(size(lsf));
propagation_mask = (abs(lsf) <= 2 * min(dx, dy)) & material_mask_full;
hj_opts = struct( ...
    'advection_order', 2, ...
    'time_integrator', 'ssprk2', ...
    'fallback_first_order', true, ...
    'stencil_mask', material_mask_full);

[lsf_new, hj_diag] = update_levelset_HJ(lsf, velocity, 0.01, dx, dy, propagation_mask, hj_opts);

core_before = lsf(2:end-1, 2:end-1);
core_after = lsf_new(2:end-1, 2:end-1);
outside_change = max(abs(core_after(~material_mask_core) - core_before(~material_mask_core)));
assert(outside_change < 1e-14, '核心材料域外不应被 HJ 更新修改。');
assert(hj_diag.second_order_count > 0, '二阶格式应至少在部分点上生效。');
assert(hj_diag.fallback_count > 0, '靠近活动域边缘时应触发一阶回退。');

fprintf('second_order_count=%d\n', hj_diag.second_order_count);
fprintf('fallback_count=%d\n', hj_diag.fallback_count);
fprintf('outside_change=%.6e\n', outside_change);
fprintf('PASS: second-order HJ material-mask test.\n');
