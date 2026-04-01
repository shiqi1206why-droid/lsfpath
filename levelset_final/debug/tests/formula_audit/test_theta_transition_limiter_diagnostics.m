clear; clc; close all;
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(script_dir));
addpath(fullfile(project_root, 'utilities'), '-begin');
project_root = get_project_root(project_root);
cleanup_path = ensure_project_on_path(project_root); %#ok<NASGU>

nelx = 14;
nely = 10;
dx = 0.04;
dy = 0.05;
delta_theta_max = deg2rad(0.5);

[X, Y] = meshgrid(0:nelx+1, 0:nely+1);
lsf = 0.2 * sin(0.2 * X) + 0.15 * cos(0.3 * Y);
material_mask_core = true(nely, nelx);
[~, theta_target] = advance_theta_state(lsf, [], delta_theta_max, dx, dy, material_mask_core, 0.1, 2);
theta_prev = mod(theta_target + deg2rad(8), pi);

transition = compute_theta_transition_forward(lsf, theta_prev, delta_theta_max, dx, dy, ...
    material_mask_core, 0.1, 2);
dtheta = ones(nely, nelx);
[~, diagnostics] = differentiate_theta_transition_exact(dtheta, transition.cache, ...
    struct('limiter_mode', 'hard', 'soft_limiter_beta', 20.0));

fprintf('saturation_ratio=%.4f\n', diagnostics.saturation_ratio);
fprintf('zero_gradient_due_to_limiter_ratio=%.4f\n', diagnostics.zero_gradient_due_to_limiter_ratio);

assert(diagnostics.saturation_ratio > 0.8, '测试场景未触发足够大的 limiter 饱和比例。');
assert(abs(diagnostics.zero_gradient_due_to_limiter_ratio - diagnostics.saturation_ratio) < 1e-12, ...
    'hard limiter 零梯度诊断应与 saturation_ratio 一致。');

fprintf('PASS: theta transition limiter diagnostics test.\n');
