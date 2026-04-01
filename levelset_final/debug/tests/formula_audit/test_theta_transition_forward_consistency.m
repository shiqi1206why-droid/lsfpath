clear; clc; close all;
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(script_dir));
addpath(fullfile(project_root, 'utilities'), '-begin');
project_root = get_project_root(project_root);
cleanup_path = ensure_project_on_path(project_root); %#ok<NASGU>

rng(7);
nelx = 18;
nely = 12;
dx = 0.04;
dy = 0.05;
delta_theta_max = deg2rad(0.6);
smooth_eta = 0.10;
smooth_iterations = 3;

[X, Y] = meshgrid(0:nelx+1, 0:nely+1);
lsf = 0.12 * sin(0.31 * X) + 0.09 * cos(0.27 * Y) + 0.04 * sin(0.11 * X .* Y);
material_mask_core = true(nely, nelx);
material_mask_core(2, 2) = false;
material_mask_core(end-1, end-1) = false;
theta_prev = mod(rand(nely, nelx) * pi, pi);
theta_prev(~material_mask_core) = NaN;

transition = compute_theta_transition_forward(lsf, theta_prev, delta_theta_max, dx, dy, ...
    material_mask_core, smooth_eta, smooth_iterations);
[theta_next, theta_target, theta_cache] = advance_theta_state(lsf, theta_prev, delta_theta_max, ...
    dx, dy, material_mask_core, smooth_eta, smooth_iterations);
[theta_raw_ref, angle_cache] = compute_fiber_angles_from_lsf(lsf, dx, dy, material_mask_core);
theta_for_smooth = theta_raw_ref;
theta_for_smooth(~material_mask_core) = 0;
[z_ref, smooth_cache] = angle_smooth_vectorized(theta_for_smooth, smooth_eta, smooth_iterations);

assert(max(abs(theta_next(:) - transition.theta_next(:)), [], 'omitnan') < 1e-12, ...
    'advance_theta_state 与共享前向核心的 theta_next 不一致。');
assert(max(abs(theta_target(:) - transition.theta_target(:)), [], 'omitnan') < 1e-12, ...
    'advance_theta_state 与共享前向核心的 theta_target 不一致。');
assert(max(abs(theta_cache.theta_raw(:) - theta_raw_ref(:)), [], 'omitnan') < 1e-12, ...
    '共享缓存中的 theta_raw 与直接计算不一致。');
assert(max(abs(theta_cache.angle_cache.dphi_dx(:) - angle_cache.dphi_dx(:)), [], 'omitnan') < 1e-12, ...
    '共享缓存中的 dphi_dx 与直接计算不一致。');
assert(max(abs(theta_cache.angle_cache.dphi_dy(:) - angle_cache.dphi_dy(:)), [], 'omitnan') < 1e-12, ...
    '共享缓存中的 dphi_dy 与直接计算不一致。');
assert(max(abs(theta_cache.z_smooth(:) - z_ref(:)), [], 'omitnan') < 1e-12, ...
    '共享缓存中的 z_smooth 与平滑算子不一致。');
assert(max(abs(theta_cache.smooth_cache.z_output(:) - smooth_cache.z_output(:)), [], 'omitnan') < 1e-12, ...
    '共享缓存中的 smooth_cache.z_output 与平滑算子缓存不一致。');

target_ref = mod(0.5 * angle(z_ref), pi);
target_ref(~material_mask_core) = NaN;
assert(max(abs(theta_target(:) - target_ref(:)), [], 'omitnan') < 1e-12, ...
    'theta_target 与 z_smooth 派生值不一致。');

transition_no_prev = compute_theta_transition_forward(lsf, [], delta_theta_max, dx, dy, ...
    material_mask_core, smooth_eta, smooth_iterations);
assert(max(abs(transition_no_prev.theta_next(:) - transition_no_prev.theta_target(:)), [], 'omitnan') < 1e-12, ...
    '未提供 theta_prev 时 theta_next 应与 theta_target 一致。');

fprintf('PASS: theta transition forward consistency test.\n');
