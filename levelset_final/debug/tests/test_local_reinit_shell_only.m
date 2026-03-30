clear; clc; close all;
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))), '-begin');

nelx = 30;
nely = 20;
dx = 0.05;
dy = 0.05;

material_mask_core = false(nely, nelx);
material_mask_core(3:18, 4:27) = true;
material_mask_full = expand_material_mask_to_full(material_mask_core);

[x_full, y_full] = get_lsf_grid_coordinates([nely + 2, nelx + 2], dx, dy);
[X, Y] = meshgrid(x_full, y_full);
cx = mean(x_full(6:end-5));
cy = mean(y_full(6:end-5));
lsf = sqrt((X - cx).^2 + (Y - cy).^2) - 0.22;
lsf(~material_mask_full) = 0.8;
lsf = impose_neumann(lsf);
lsf_distorted = lsf + 0.08 * min(dx, dy) * sin(5 * pi * X) .* cos(4 * pi * Y);
lsf_distorted(~material_mask_full) = 0.8;
lsf_distorted = impose_neumann(lsf_distorted);

zero_mask = abs(lsf_distorted) <= 0.5 * min(dx, dy);
local_shell_mask = imdilate(abs(lsf_distorted) <= 2 * min(dx, dy), strel('square', 7)) & material_mask_full;

[lsf_reinit, reinit_diag] = fmm_reinitialize(lsf_distorted, dx, dy, zero_mask, material_mask_core, ...
    struct('method', 'subcell_signed_distance', 'local_shell_mask', local_shell_mask, ...
           'preserve_outside_shell', true));

interior_outside_shell = false(size(lsf_reinit));
interior_outside_shell(2:end-1, 2:end-1) = true;
interior_outside_shell = interior_outside_shell & material_mask_full & ~local_shell_mask;
outside_diff = max(abs(lsf_reinit(interior_outside_shell) - lsf_distorted(interior_outside_shell)));

assert(outside_diff < 1e-12, '局部重初始化不应修改 shell 外部的内部材料域。');
assert(reinit_diag.local_shell_applied, '应记录 local shell 已应用。');

fprintf('local_shell_size=%d\n', reinit_diag.local_shell_size);
fprintf('outside_diff=%.6e\n', outside_diff);
fprintf('PASS: local reinit shell only test.\n');

function field = impose_neumann(field)
field(1, :) = field(2, :);
field(end, :) = field(end-1, :);
field(:, 1) = field(:, 2);
field(:, end) = field(:, end-1);
end
