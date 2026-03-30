clc;
close all;
set(0, 'DefaultFigureVisible', 'off');

project_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(project_dir, 'utilities'), '-begin');
project_root = get_project_root(project_dir);
cleanup_path = ensure_project_on_path(project_root); %#ok<NASGU>
paths = build_project_paths(project_root);

stamp = datestr(now, 'yyyymmdd_HHMMSS');
out_dir = fullfile(paths.baseline_dir, ['default_run_' stamp]);
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

results = fiber_levelset('default');

save(fullfile(out_dir, 'results.mat'), 'results', '-v7.3');

fig1 = figure('Visible', 'off', 'Position', [80, 80, 1200, 700]);
yyaxis left;
plot(results.compliance_history, 'b-', 'LineWidth', 1.6, 'DisplayName', 'Compliance');
hold on;
ylabel('Compliance');
xlabel('Iteration');
grid on;
yyaxis right;
plot(results.FCS_history * 100, 'r-', 'LineWidth', 1.2, 'DisplayName', 'FCS');
ylabel('FCS (%)');
title('Default Run: Compliance and FCS');
legend('Location', 'best');
print(fig1, fullfile(out_dir, 'convergence.png'), '-dpng', '-r140');
close(fig1);

fig2 = figure('Visible', 'off', 'Position', [100, 100, 900, 600]);
material_mask_full = results.material_mask_full;
lsf_plot = results.lsf;
lsf_plot(~material_mask_full) = NaN;
contour(lsf_plot, 25, 'LineWidth', 0.8);
hold on;
contour(lsf_plot, [0 0], 'r', 'LineWidth', 2);
axis equal tight;
set(gca, 'YDir', 'reverse');
title('Default Run: Final Level Set');
xlabel('x index');
ylabel('y index');
print(fig2, fullfile(out_dir, 'final_lsf.png'), '-dpng', '-r140');
close(fig2);

levels = (-2:2) * results.params.grid.h;
opts = struct();
opts.target_ds = results.params.grid.h / 4;
opts.smooth_window = 9;
opts.smooth_iters = 3;
opts.min_points = 20;
export_printable_paths_from_lsf(results.lsf, results.params.grid.dx, results.params.grid.dy, ...
    levels, fullfile(out_dir, 'printable_paths'), results.material_mask_core, opts);

report_path = fullfile(out_dir, 'summary.txt');
fid = fopen(report_path, 'w');
fprintf(fid, 'out_dir=%s\n', out_dir);
fprintf(fid, 'final_compliance=%.12e\n', results.final_compliance);
fprintf(fid, 'final_FCS=%.6f\n', results.final_FCS);
fprintf(fid, 'final_iter=%d\n', results.final_iter);
fclose(fid);

fprintf('OUT_DIR=%s\n', out_dir);
