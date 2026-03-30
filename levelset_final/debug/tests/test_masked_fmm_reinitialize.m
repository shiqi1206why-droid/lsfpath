clear; clc; close all;
addpath(genpath(fileparts(fileparts(mfilename('fullpath')))), '-begin');

nelx = 16;
nely = 12;
dx = 0.1;
dy = 0.1;

material_mask_core = false(nely, nelx);
material_mask_core(3:10, 4:14) = true;
material_mask_full = expand_material_mask_to_full(material_mask_core);

[X, ~] = meshgrid(0:nelx+1, 0:nely+1);
lsf = (X - 8.5) * dx;

zero_mask = false(size(lsf));
zero_mask(:, 9) = true;
outside_seed_count = nnz(zero_mask & ~material_mask_full);

[lsf_reinit, diag_info] = fmm_reinitialize(lsf, dx, dy, zero_mask, material_mask_core, ...
    struct('method', 'masked_fmm'));

zero_outside = nnz(abs(lsf_reinit) <= 1e-12 & ~material_mask_full);
nonpositive_outside = nnz(lsf_reinit(~material_mask_full) <= 0);
zero_inside = nnz(abs(lsf_reinit) <= 1e-12 & material_mask_full);

fprintf('outside_seed_count=%d\n', outside_seed_count);
fprintf('zero_outside=%d\n', zero_outside);
fprintf('nonpositive_outside=%d\n', nonpositive_outside);
fprintf('zero_inside=%d\n', zero_inside);
fprintf('used_method=%s\n', diag_info.method_used);

assert(outside_seed_count > 0, '测试场景无效：需要先在void域布置种子。');
assert(zero_outside == 0, 'masked FMM 不应在void域保留零种子。');
assert(nonpositive_outside == 0, 'masked FMM 后，void域必须全部为严格正值。');
assert(zero_inside > 0, 'masked FMM 后，材料域内应保留零水平集。');
assert(strcmpi(diag_info.method_used, 'masked_fmm_fallback'), '测试应显式走 masked_fmm 路径。');

fprintf('PASS: masked FMM reinitialization test.\n');

function mask_full = expand_material_mask_to_full(mask_core)
mask_full = false(size(mask_core, 1) + 2, size(mask_core, 2) + 2);
mask_full(2:end-1, 2:end-1) = logical(mask_core);
mask_full(1, :) = mask_full(2, :);
mask_full(end, :) = mask_full(end-1, :);
mask_full(:, 1) = mask_full(:, 2);
mask_full(:, end) = mask_full(:, end-1);
end
