clc; close all;

script_path = mfilename('fullpath');
if isempty(script_path)
    script_path = which('run_default_precise_rerun_and_stitch.m');
end
script_dir = fileparts(script_path);
project_dir = fileparts(fileparts(script_dir));
addpath(fullfile(project_dir, 'utilities'), '-begin');
project_root = get_project_root(project_dir);
cleanup_path = ensure_project_on_path(project_root); %#ok<NASGU>
set(0, 'DefaultFigureVisible', 'off');

paths = build_project_paths(project_root);
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
artifact_dir = fullfile(paths.refactor_dir, ['rerun_default_precise_compare_' timestamp]);
if ~exist(artifact_dir, 'dir')
    mkdir(artifact_dir);
end

configs = {'default', 'precise'};
template = struct( ...
    'config', '', ...
    'chain_mode', '', ...
    'initial_compliance', NaN, ...
    'final_compliance', NaN, ...
    'best_compliance', NaN, ...
    'final_FCS', NaN, ...
    'best_iter', NaN, ...
    'final_iter', NaN, ...
    'executed_iter', NaN, ...
    'accepted_steps', NaN, ...
    'rejected_steps', NaN, ...
    'rollback_to_best', false, ...
    'early_stop_triggered', false, ...
    'early_stop_reason', '', ...
    'compliance_history', [], ...
    'FCS_history', [], ...
    'lsf', [], ...
    'material_mask_full', []);
summaries = repmat(template, numel(configs), 1);

for i = 1:numel(configs)
    cfg = configs{i};
    fprintf('=== Running %s ===\n', cfg);
    addpath(genpath(project_root), '-begin');
    rehash;
    log_path = fullfile(artifact_dir, sprintf('rerun_%s_log.txt', cfg));
    diary(log_path);
    diary on;
    try
        results = fiber_levelset(cfg); %#ok<NASGU>
    catch ME
        diary off;
        rethrow(ME);
    end
    diary off;

    results_path = fullfile(artifact_dir, sprintf('rerun_%s_results.mat', cfg));
    save(results_path, 'results');

    entry = template;
    entry.config = cfg;
    entry.chain_mode = char(results.params.gradient.chain_mode);
    entry.initial_compliance = results.compliance_history(1);
    entry.final_compliance = results.final_compliance;
    entry.best_compliance = results.best_compliance;
    entry.final_FCS = results.final_FCS;
    entry.best_iter = results.best_iter;
    entry.final_iter = results.final_iter;
    entry.executed_iter = results.executed_iter;
    entry.accepted_steps = results.accepted_steps;
    entry.rejected_steps = results.rejected_steps;
    entry.rollback_to_best = logical(results.rollback_to_best);
    entry.early_stop_triggered = logical(results.early_stop_triggered);
    if isfield(results, 'early_stop_reason') && ~isempty(results.early_stop_reason)
        entry.early_stop_reason = results.early_stop_reason;
    end
    entry.compliance_history = results.compliance_history(:);
    entry.FCS_history = results.FCS_history(:);
    entry.lsf = results.lsf;
    entry.material_mask_full = results.material_mask_full;
    summaries(i) = entry;
end

% === Per-run figures ===
for i = 1:numel(summaries)
    s = summaries(i);

    fig_conv = figure('Color', 'w', 'Position', [100, 100, 820, 460]);
    yyaxis left;
    plot(0:(numel(s.compliance_history)-1), s.compliance_history, 'b-', 'LineWidth', 1.8);
    ylabel('Compliance');
    yyaxis right;
    plot(0:(numel(s.FCS_history)-1), s.FCS_history * 100, 'r-', 'LineWidth', 1.6);
    ylabel('FCS (%)');
    xlabel('Iteration');
    grid on;
    title(sprintf('%s convergence (chain=%s)', s.config, s.chain_mode));
    exportgraphics(fig_conv, fullfile(artifact_dir, sprintf('%s_convergence.png', s.config)), 'Resolution', 180);
    close(fig_conv);

    fig_lsf = figure('Color', 'w', 'Position', [120, 120, 820, 460]);
    lsf_plot = s.lsf;
    lsf_plot(~s.material_mask_full) = NaN;
    contour(lsf_plot, 25, 'LineWidth', 0.6, 'LineColor', [0.2 0.4 0.9]);
    hold on;
    contour(lsf_plot, [0 0], 'r', 'LineWidth', 2.0);
    axis equal tight;
    set(gca, 'YDir', 'reverse');
    grid on;
    title(sprintf('%s final lsf (chain=%s)', s.config, s.chain_mode));
    xlabel('X index');
    ylabel('Y index');
    exportgraphics(fig_lsf, fullfile(artifact_dir, sprintf('%s_final_lsf.png', s.config)), 'Resolution', 180);
    close(fig_lsf);
end

% === Stitched comparison ===
fig_cmp = figure('Color', 'w', 'Position', [60, 60, 1500, 1000]);
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

for i = 1:numel(summaries)
    s = summaries(i);
    nexttile(i);
    yyaxis left;
    plot(0:(numel(s.compliance_history)-1), s.compliance_history, 'b-', 'LineWidth', 1.8);
    ylabel('Compliance');
    yyaxis right;
    plot(0:(numel(s.FCS_history)-1), s.FCS_history * 100, 'r-', 'LineWidth', 1.6);
    ylabel('FCS (%)');
    xlabel('Iteration');
    grid on;
    title(sprintf('%s convergence', s.config));
end

for i = 1:numel(summaries)
    s = summaries(i);
    nexttile(numel(summaries) + i);
    lsf_plot = s.lsf;
    lsf_plot(~s.material_mask_full) = NaN;
    contour(lsf_plot, 25, 'LineWidth', 0.6, 'LineColor', [0.2 0.4 0.9]);
    hold on;
    contour(lsf_plot, [0 0], 'r', 'LineWidth', 2.0);
    axis equal tight;
    set(gca, 'YDir', 'reverse');
    grid on;
    title(sprintf('%s final lsf', s.config));
    xlabel('X index');
    ylabel('Y index');
end

sgtitle('Default vs Precise Comparison');
stitched_path = fullfile(artifact_dir, 'default_vs_precise_stitched.png');
exportgraphics(fig_cmp, stitched_path, 'Resolution', 220);
close(fig_cmp);

summary_txt = fullfile(artifact_dir, 'summary.txt');
fid = fopen(summary_txt, 'w');
if fid < 0
    error('无法写入摘要文件: %s', summary_txt);
end
cleanup_summary = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, 'run_tag=rerun_default_precise_compare_%s\n', timestamp);
fprintf(fid, 'artifact_dir=%s\n', artifact_dir);
fprintf(fid, 'stitched_image=%s\n', stitched_path);
for i = 1:numel(summaries)
    s = summaries(i);
    fprintf(fid, '\n[%s]\n', s.config);
    fprintf(fid, 'chain_mode=%s\n', s.chain_mode);
    fprintf(fid, 'initial_compliance=%.12e\n', s.initial_compliance);
    fprintf(fid, 'final_compliance=%.12e\n', s.final_compliance);
    fprintf(fid, 'best_compliance=%.12e\n', s.best_compliance);
    fprintf(fid, 'final_FCS=%.12f\n', s.final_FCS);
    fprintf(fid, 'best_iter=%d\n', s.best_iter);
    fprintf(fid, 'final_iter=%d\n', s.final_iter);
    fprintf(fid, 'executed_iter=%d\n', s.executed_iter);
    fprintf(fid, 'accepted_steps=%d\n', s.accepted_steps);
    fprintf(fid, 'rejected_steps=%d\n', s.rejected_steps);
    fprintf(fid, 'rollback_to_best=%d\n', s.rollback_to_best);
    fprintf(fid, 'early_stop_triggered=%d\n', s.early_stop_triggered);
    fprintf(fid, 'early_stop_reason=%s\n', s.early_stop_reason);
end
clear cleanup_summary

summary_json = fullfile(artifact_dir, 'summary.json');
fid_json = fopen(summary_json, 'w');
if fid_json < 0
    error('无法写入JSON文件: %s', summary_json);
end
cleanup_json = onCleanup(@() fclose(fid_json)); %#ok<NASGU>
fprintf(fid_json, '%s', jsonencode(summaries));
clear cleanup_json

save(fullfile(artifact_dir, 'results.mat'), 'summaries');

fprintf('ARTIFACT_DIR=%s\n', artifact_dir);
fprintf('STITCHED=%s\n', stitched_path);
fprintf('SUMMARY_TXT=%s\n', summary_txt);
