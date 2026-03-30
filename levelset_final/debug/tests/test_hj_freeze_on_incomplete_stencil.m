clear; clc; close all;
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))), '-begin');

ny = 9;
nx = 9;
dx = 0.05;
dy = 0.05;
[x_full, y_full] = get_lsf_grid_coordinates([ny, nx], dx, dy);
[X, Y] = meshgrid(x_full, y_full);
lsf = X + 0.2 * Y;
velocity = ones(size(lsf));

active_mask = false(size(lsf));
active_mask(5, 5) = true;
stencil_mask = true(size(lsf));
stencil_mask(5, 4) = false;  % 缺失 Godunov 所需左邻点

hj_opts = struct( ...
    'advection_order', 2, ...
    'time_integrator', 'ssprk2', ...
    'fallback_first_order', true, ...
    'freeze_on_incomplete_godunov', true, ...
    'stencil_mask', stencil_mask);

[lsf_new, hj_diag] = update_levelset_HJ(lsf, velocity, 0.02, dx, dy, active_mask, hj_opts);

assert(abs(lsf_new(5, 5) - lsf(5, 5)) < 1e-14, '不完整 stencil 点应被冻结。');
assert(hj_diag.frozen_incomplete_godunov_count >= 1, '应记录冻结的 Godunov 点数。');

fprintf('frozen_incomplete_godunov_count=%d\n', hj_diag.frozen_incomplete_godunov_count);
fprintf('PASS: HJ freeze on incomplete stencil test.\n');
