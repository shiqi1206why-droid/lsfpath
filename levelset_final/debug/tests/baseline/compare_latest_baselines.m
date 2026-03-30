% 一键比较最近两次 baseline 指标
% 用法：直接运行本脚本（无输入）

clc;
fprintf('=== Compare Latest Baselines ===\n');

script_path = mfilename('fullpath');
if isempty(script_path)
    script_path = which('compare_latest_baselines.m');
end
if isempty(script_path)
    error('无法定位 compare_latest_baselines.m');
end

script_dir = fileparts(script_path);          % .../tests/baseline
project_dir = fileparts(fileparts(script_dir)); % .../debug
addpath(fullfile(project_dir, 'utilities'), '-begin');
project_root = get_project_root(project_dir);
cleanup_path = ensure_project_on_path(project_root); %#ok<NASGU>
paths = build_project_paths(project_root);
output_root = paths.baseline_dir;

if ~exist(output_root, 'dir')
    error('未找到 baseline_artifacts: %s', output_root);
end

d = dir(fullfile(output_root, 'baseline_*'));
d = d([d.isdir]);
d = d(~ismember({d.name}, {'.', '..'}));

if numel(d) < 2
    error('可比较的 baseline 目录不足 2 个。');
end

[~, idx] = sort({d.name});
d = d(idx);
prev_dir = fullfile(output_root, d(end-1).name);
curr_dir = fullfile(output_root, d(end).name);

prev_json = fullfile(prev_dir, 'baseline_metrics.json');
curr_json = fullfile(curr_dir, 'baseline_metrics.json');

if ~exist(prev_json, 'file') || ~exist(curr_json, 'file')
    error('缺少 baseline_metrics.json：\n  prev=%s\n  curr=%s', prev_json, curr_json);
end

prev = jsondecode(fileread(prev_json));
curr = jsondecode(fileread(curr_json));

delta_final_compliance = curr.optimization.final_compliance - prev.optimization.final_compliance;
delta_improvement_ratio = curr.optimization.improvement_ratio - prev.optimization.improvement_ratio;
delta_fcs = curr.optimization.final_FCS - prev.optimization.final_FCS;
delta_runtime = curr.runtime_seconds - prev.runtime_seconds;

comparison = struct();
comparison.prev_run = d(end-1).name;
comparison.curr_run = d(end).name;
comparison.prev = prev.optimization;
comparison.curr = curr.optimization;
comparison.delta = struct( ...
    'final_compliance', delta_final_compliance, ...
    'improvement_ratio_percent', delta_improvement_ratio, ...
    'final_FCS', delta_fcs, ...
    'runtime_seconds', delta_runtime);

cmp_json_path = fullfile(curr_dir, 'compare_latest_metrics.json');
fid = fopen(cmp_json_path, 'w');
if fid == -1
    error('无法写入: %s', cmp_json_path);
end
fprintf(fid, '%s', jsonencode(comparison));
fclose(fid);

cmp_md_path = fullfile(curr_dir, 'compare_latest_report.md');
fid = fopen(cmp_md_path, 'w');
if fid == -1
    error('无法写入: %s', cmp_md_path);
end
fprintf(fid, '# Baseline Compare (Latest Two Runs)\n\n');
fprintf(fid, '- Previous run: `%s`\n', d(end-1).name);
fprintf(fid, '- Current run: `%s`\n\n', d(end).name);
fprintf(fid, '## Core Metrics\n\n');
fprintf(fid, '| Metric | Previous | Current | Delta |\n');
fprintf(fid, '|---|---:|---:|---:|\n');
fprintf(fid, '| Final compliance | %.10e | %.10e | %.10e |\n', ...
    prev.optimization.final_compliance, curr.optimization.final_compliance, delta_final_compliance);
fprintf(fid, '| Improvement ratio (%%) | %.6f | %.6f | %.6f |\n', ...
    prev.optimization.improvement_ratio, curr.optimization.improvement_ratio, delta_improvement_ratio);
fprintf(fid, '| Final FCS | %.6f | %.6f | %.6f |\n', ...
    prev.optimization.final_FCS, curr.optimization.final_FCS, delta_fcs);
fprintf(fid, '| Runtime (s) | %.3f | %.3f | %.3f |\n', ...
    prev.runtime_seconds, curr.runtime_seconds, delta_runtime);
fclose(fid);

fprintf('Previous: %s\n', prev_dir);
fprintf('Current : %s\n', curr_dir);
fprintf('Delta final compliance: %.10e\n', delta_final_compliance);
fprintf('Delta improvement ratio: %.6f %%\n', delta_improvement_ratio);
fprintf('Delta final FCS: %.6f\n', delta_fcs);
fprintf('Delta runtime: %.3f s\n', delta_runtime);
fprintf('Report: %s\n', cmp_md_path);
fprintf('JSON  : %s\n', cmp_json_path);
