function [next_state, accepted_source, accepted_source_detail] = fiber_select_candidate_state(sel_in)
%FIBER_SELECT_CANDIDATE_STATE Select accepted state among theta-only, HJ, and reinit candidates.

    accepted_source = 'hold';
    accepted_source_detail = 'hold';
    next_state = sel_in.current_state;

    if sel_in.theta_only_accepted && isfinite(sel_in.theta_only_compliance)
        next_state = sel_in.theta_only_state;
        accepted_source = 'theta_only';
        accepted_source_detail = 'theta_only';
    end

    hj_state_candidate = [];
    hj_candidate_detail = 'hj';
    if sel_in.step_accepted
        if isfield(sel_in, 'hj_local_reinit_state') && ~isempty(sel_in.hj_local_reinit_state)
            hj_state_candidate = sel_in.hj_local_reinit_state;
            hj_candidate_detail = 'hj_local_reinit';
        elseif isfield(sel_in, 'hj_raw_state') && ~isempty(sel_in.hj_raw_state)
            hj_state_candidate = sel_in.hj_raw_state;
            hj_candidate_detail = 'hj_raw';
        elseif isfield(sel_in, 'hj_trial_state') && ~isempty(sel_in.hj_trial_state)
            hj_state_candidate = sel_in.hj_trial_state;
            hj_candidate_detail = 'hj_trial';
        end
    end

    if ~isempty(hj_state_candidate) && ...
            (strcmp(accepted_source, 'hold') || ...
             hj_state_candidate.compliance < next_state.compliance * (1 - sel_in.candidate_select_tol))
        next_state = hj_state_candidate;
        accepted_source = 'hj';
        accepted_source_detail = hj_candidate_detail;
    end

    if sel_in.reinit_candidate_selected && ~isempty(sel_in.reinit_trial_state) && ...
            (strcmp(accepted_source, 'hold') || ...
             sel_in.reinit_trial_state.compliance < next_state.compliance * (1 - sel_in.candidate_select_tol) || ...
              (strcmp(accepted_source, 'hj') && ...
              ~isempty(hj_state_candidate) && ...
              sel_in.reinit_trial_state.compliance <= hj_state_candidate.compliance * (1 + sel_in.acceptance_tol)))
        next_state = sel_in.reinit_trial_state;
        accepted_source = 'reinit';
        accepted_source_detail = 'reinit';
    end
end
