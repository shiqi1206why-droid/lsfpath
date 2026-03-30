function step_out = fiber_run_hj_backtracking_step(step_in)
%FIBER_RUN_HJ_BACKTRACKING_STEP Evaluate HJ candidate with guard/backtracking policy.

    step_out = struct();
    step_out.hj_trial_compliance = NaN;
    step_out.hj_trial_state = [];
    step_out.step_accepted = true;
    step_out.backtrack_used = 0;
    step_out.fail_next_guard = false;
    step_out.fail_current_guard = false;
    step_out.hj_info_trial = struct();

    accepted_steps = step_in.accepted_steps;
    rejected_steps = step_in.rejected_steps;
    reject_due_next_guard = step_in.reject_due_next_guard;
    reject_due_current_guard = step_in.reject_due_current_guard;

    if step_in.enable_step_acceptance
        step_out.step_accepted = false;
        for bt = 0:step_in.max_backtrack
            dt_trial = step_in.dt_adaptive * (step_in.backtrack_factor ^ bt);
            if dt_trial < step_in.min_backtrack_dt
                break;
            end

            hj_update_opts = step_in.hj_update_opts;
            hj_update_opts.stencil_mask = step_in.stencil_mask;
            [lsf_trial, hj_info_trial] = update_levelset_HJ(step_in.lsf_before, step_in.velocity, dt_trial, ...
                step_in.dx, step_in.dy, step_in.primary_update_mask, hj_update_opts);
            if any(~isfinite(lsf_trial(:)))
                continue;
            end

            hj_candidate = evaluate_candidate_state(lsf_trial, step_in.theta_e, step_in.delta_theta_max, ...
                step_in.dx, step_in.dy, step_in.nelx, step_in.nely, step_in.material_mask_core, ...
                step_in.E_L, step_in.E_T, step_in.nu_LT, step_in.G_LT, step_in.thickness, ...
                step_in.F_mag, step_in.smooth_eta, step_in.smooth_iterations);
            hj_trial_compliance = hj_candidate.compliance;

            trial_ok_next = isfinite(hj_trial_compliance) && ...
                hj_trial_compliance <= step_in.ref_next_compliance * (1 + step_in.acceptance_tol);
            if step_in.current_guard_active
                trial_ok_current = isfinite(hj_trial_compliance) && ...
                    hj_trial_compliance <= step_in.compliance * (1 + step_in.current_state_tol);
            else
                trial_ok_current = true;
            end

            if trial_ok_next && trial_ok_current
                step_out.step_accepted = true;
                step_out.backtrack_used = bt;
                step_out.hj_trial_compliance = hj_trial_compliance;
                step_out.hj_trial_state = hj_candidate;
                step_out.hj_info_trial = hj_info_trial;
                step_out.dt_trial = dt_trial;
                break;
            else
                step_out.fail_next_guard = step_out.fail_next_guard || ~trial_ok_next;
                step_out.fail_current_guard = step_out.fail_current_guard || ~trial_ok_current;
                step_out.hj_trial_compliance = hj_trial_compliance;
                step_out.hj_trial_state = hj_candidate;
                step_out.hj_info_trial = hj_info_trial;
                step_out.dt_trial = dt_trial;
            end
        end

        if step_out.step_accepted
            accepted_steps = accepted_steps + 1;
            if step_out.backtrack_used > 0 && (step_in.iter == 1 || mod(step_in.iter, 10) == 0)
                fprintf('  [步长回溯] 接受更新：回溯%d次，dt=%.3e，C_trial=%.4e，C_ref_next=%.4e\n', ...
                    step_out.backtrack_used, step_out.dt_trial, step_out.hj_trial_compliance, step_in.ref_next_compliance);
            end
        else
            rejected_steps = rejected_steps + 1;
            if step_out.fail_next_guard
                reject_due_next_guard = reject_due_next_guard + 1;
            end
            if step_out.fail_current_guard
                reject_due_current_guard = reject_due_current_guard + 1;
            end
            if step_in.iter == 1 || mod(step_in.iter, 10) == 0
                fprintf('  [步长回溯] 拒绝更新：最小dt=%.1e，C_trial=%.4e，C_ref_next=%.4e，C_current=%.4e，guard_active=%d，fail_next=%d，fail_current=%d\n', ...
                    step_in.min_backtrack_dt, step_out.hj_trial_compliance, step_in.ref_next_compliance, step_in.compliance, ...
                    step_in.current_guard_active, step_out.fail_next_guard, step_out.fail_current_guard);
            end
        end
    else
        hj_update_opts = step_in.hj_update_opts;
        hj_update_opts.stencil_mask = step_in.stencil_mask;
        [lsf_trial, hj_info_trial] = update_levelset_HJ(step_in.lsf_before, step_in.velocity, step_in.dt_adaptive, ...
            step_in.dx, step_in.dy, step_in.primary_update_mask, hj_update_opts);

        hj_trial_state = evaluate_candidate_state(lsf_trial, step_in.theta_e, step_in.delta_theta_max, ...
            step_in.dx, step_in.dy, step_in.nelx, step_in.nely, step_in.material_mask_core, ...
            step_in.E_L, step_in.E_T, step_in.nu_LT, step_in.G_LT, step_in.thickness, ...
            step_in.F_mag, step_in.smooth_eta, step_in.smooth_iterations);

        step_out.hj_trial_state = hj_trial_state;
        step_out.hj_trial_compliance = hj_trial_state.compliance;
        step_out.hj_info_trial = hj_info_trial;
        step_out.step_accepted = true;
        accepted_steps = accepted_steps + 1;
    end

    step_out.accepted_steps = accepted_steps;
    step_out.rejected_steps = rejected_steps;
    step_out.reject_due_next_guard = reject_due_next_guard;
    step_out.reject_due_current_guard = reject_due_current_guard;
end
