clear; clc;
load('\\wsl.localhost\Ubuntu\home\again\projects\recover\debug\refactor_artifacts\rerun_default_20260331_170739\rerun_default_results.mat','results');
base = results.compliance_history(:);
valid_hist = base(isfinite(base));
fprintf('initial=%.12e\n', valid_hist(1));
fprintf('best=%.12e best_iter=%d\n', results.best_compliance, results.best_iter);
fprintf('final=%.12e final_iter=%d\n', results.final_compliance, results.final_iter);
if isfield(results,'raw_final_compliance')
  fprintf('raw_final=%.12e\n', results.raw_final_compliance);
end
if isfield(results,'final_to_best_gap_percent')
  fprintf('final_to_best_gap_percent=%.6f\n', results.final_to_best_gap_percent);
end
src = string(results.accepted_source_detail_history(:));
nonempty = src(src ~= "");
if ~isempty(nonempty)
  u = unique(nonempty);
  fprintf('accepted_source_detail unique:\n');
  disp(u);
end
fprintf('accepted=%d rejected=%d\n', results.accepted_steps, results.rejected_steps);
for name = {'theta_only_compliance_history','hj_trial_compliance_history','hj_local_reinit_compliance_history','reinit_trial_compliance_history'}
  f = name{1};
  if isfield(results,f)
    v = results.(f); v = v(:); v = v(isfinite(v));
    if ~isempty(v)
      fprintf('%s: min=%.12e max=%.12e first=%.12e last=%.12e\n', f, min(v), max(v), v(1), v(end));
    end
  end
end
hist = results.compliance_history(:);
N = numel(hist);
[best_val, best_idx] = min(hist);
[last_worse, idx] = max(hist(1:min(end, results.final_iter)));
fprintf('hist_best=%.12e hist_best_idx=%d\n', best_val, best_idx);
fprintf('hist_max=%.12e hist_max_idx=%d\n', last_worse, idx);
for i = 1:min(results.final_iter, numel(results.accepted_source_detail_history))
  src_i = results.accepted_source_detail_history(i);
  if strlength(src_i) == 0
    continue;
  end
  cur = results.compliance_history(min(i+1,numel(results.compliance_history)));
  theta = results.theta_only_compliance_history(i);
  hj = results.hj_trial_compliance_history(i);
  fprintf('iter=%d src=%s current_after=%.12e theta=%.12e hj=%.12e\n', i, src_i, cur, theta, hj);
end
exit;
