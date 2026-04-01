clear; clc; close all;
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(script_dir));
addpath(fullfile(project_root, 'utilities'), '-begin');
project_root = get_project_root(project_root);
cleanup_path = ensure_project_on_path(project_root); %#ok<NASGU>

rng(11);
nelx = 10;
nely = 9;
dx = 0.08;
dy = 0.06;
delta_theta_max = deg2rad(0.4);
smooth_eta = 0.12;
smooth_iterations = 2;

[X, Y] = meshgrid(0:nelx+1, 0:nely+1);
lsf = 0.03 * X + 0.05 * Y + 0.02 * sin(0.6 * X) - 0.04 * cos(0.3 * Y);
material_mask_core = true(nely, nelx);
material_mask_core(1, :) = false;
theta_prev = mod(rand(nely, nelx) * pi, pi);
theta_prev(~material_mask_core) = NaN;

transition = compute_theta_transition_forward(lsf, theta_prev, delta_theta_max, dx, dy, ...
    material_mask_core, smooth_eta, smooth_iterations);
cache = transition.cache;

assert(isequal(size(cache.phi_base), size(lsf)), 'phi_base 尺寸错误。');
assert(isequal(size(cache.theta_prev), [nely, nelx]), 'theta_prev 尺寸错误。');
assert(isequal(size(cache.theta_raw), [nely, nelx]), 'theta_raw 尺寸错误。');
assert(isequal(size(cache.theta_target), [nely, nelx]), 'theta_target 尺寸错误。');
assert(isequal(size(cache.theta_next), [nely, nelx]), 'theta_next 尺寸错误。');
assert(isequal(size(cache.limiter_linear_mask), [nely, nelx]), 'limiter_linear_mask 尺寸错误。');
assert(isequal(size(cache.limiter_saturated_mask), [nely, nelx]), 'limiter_saturated_mask 尺寸错误。');
assert(isequal(size(cache.degenerate_grad_mask), [nely, nelx]), 'degenerate_grad_mask 尺寸错误。');
assert(numel(cache.smooth_cache.iterations) == smooth_iterations, '平滑缓存迭代数错误。');

for k = 1:smooth_iterations
    iter_cache = cache.smooth_cache.iterations(k);
    assert(isequal(size(iter_cache.z_before_conv), [nely, nelx]), 'z_before_conv 尺寸错误。');
    assert(isequal(size(iter_cache.z_after_conv), [nely, nelx]), 'z_after_conv 尺寸错误。');
    assert(isequal(size(iter_cache.z_after_bc), [nely, nelx]), 'z_after_bc 尺寸错误。');
    assert(isequal(size(iter_cache.norm_before_clamp), [nely, nelx]), 'norm_before_clamp 尺寸错误。');
    assert(isequal(size(iter_cache.clamped_norm), [nely, nelx]), 'clamped_norm 尺寸错误。');
    assert(isequal(size(iter_cache.clamped_mask), [nely, nelx]), 'clamped_mask 尺寸错误。');
    assert(isequal(size(iter_cache.z_after_normalize), [nely, nelx]), 'z_after_normalize 尺寸错误。');

    z_reconstructed = iter_cache.z_after_bc ./ iter_cache.clamped_norm;
    assert(max(abs(z_reconstructed(:) - iter_cache.z_after_normalize(:)), [], 'omitnan') < 1e-12, ...
        '归一化缓存重建失败。');
end

theta_step_ref = sign(cache.theta_diff_wrapped) .* min(abs(cache.theta_diff_wrapped), cache.delta_theta_max);
theta_next_ref = mod(cache.theta_prev + theta_step_ref, pi);
theta_next_ref(~material_mask_core) = NaN;
assert(max(abs(theta_next_ref(:) - cache.theta_next(:)), [], 'omitnan') < 1e-12, ...
    'limiter 缓存重建 theta_next 失败。');

partition_mask = cache.limiter_linear_mask | cache.limiter_saturated_mask | cache.limiter_kink_mask;
assert(nnz(partition_mask(material_mask_core)) == nnz(material_mask_core), ...
    'limiter 区域划分未覆盖全部材料域。');

fprintf('PASS: theta transition cache self-check test.\n');

