% Step4 narrow-band aggregation test
clc;
addpath(genpath(fileparts(fileparts(fileparts(mfilename('fullpath'))))), '-begin');

nelx = 12;
nely = 8;
dx = 1.0;
dy = 1.0;

rng(42);
element_sensitivity = randn(nely, nelx);
theta_field = zeros(nely, nelx); %#ok<NASGU>

% Build a simple lsf and a central narrow-band update mask
[X, Y] = meshgrid(0:nelx+1, 0:nely+1);
lsf = (X - (nelx+2)/2) * 0.05 + (Y - (nely+2)/2) * 0.02;
material_mask_core = false(nely, nelx);
material_mask_core(2:end-1, 3:end-2) = true;
material_mask_full = expand_material_mask_to_full(material_mask_core);
primary_update_mask = false(size(lsf));
primary_update_mask(3:end-2, 4:end-3) = true;
primary_update_mask = primary_update_mask & material_mask_full;

node_sens = aggregate_node_sensitivity(element_sensitivity, theta_field, lsf, nelx, nely, dx, dy, primary_update_mask);

outside_nonzero = nnz(abs(node_sens(~primary_update_mask)) > 1e-14);
inside_nonzero = nnz(abs(node_sens(primary_update_mask)) > 1e-14);
void_nonzero = nnz(abs(node_sens(~material_mask_full)) > 1e-14);

fprintf('outside_nonzero=%d\n', outside_nonzero);
fprintf('inside_nonzero=%d\n', inside_nonzero);
fprintf('void_nonzero=%d\n', void_nonzero);

assert(outside_nonzero == 0, 'Sensitivity outside update mask must be zero.');
assert(inside_nonzero > 0, 'Sensitivity inside update mask should be non-zero.');
assert(void_nonzero == 0, 'Sensitivity in void cells must stay zero.');

fprintf('PASS: Step4 narrow-band sensitivity aggregation test.\n');

function mask_full = expand_material_mask_to_full(mask_core)
mask_full = false(size(mask_core, 1) + 2, size(mask_core, 2) + 2);
mask_full(2:end-1, 2:end-1) = logical(mask_core);
mask_full(1, :) = mask_full(2, :);
mask_full(end, :) = mask_full(end-1, :);
mask_full(:, 1) = mask_full(:, 2);
mask_full(:, end) = mask_full(:, end-1);
end
