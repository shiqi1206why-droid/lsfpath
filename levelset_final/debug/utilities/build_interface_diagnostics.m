function interface_diagnostics = build_interface_diagnostics(input_data)
%BUILD_INTERFACE_DIAGNOSTICS Assemble interface diagnostics for results.

    interface_diagnostics = struct();
    interface_diagnostics.path_quality_history = input_data.path_quality_history;
    interface_diagnostics.hj_scheme = sprintf('order%d_%s', ...
        input_data.params.levelset.advection_order, input_data.params.levelset.time_integrator);
    interface_diagnostics.hj_fallback_count_history = input_data.hj_fallback_history;
    interface_diagnostics.hj_second_order_count_history = input_data.hj_second_order_history;
    interface_diagnostics.hj_frozen_incomplete_history = input_data.hj_frozen_incomplete_history;
    interface_diagnostics.hj_first_order_complete_history = input_data.hj_first_order_complete_history;
    interface_diagnostics.reinit_method_history = cellstr(input_data.reinit_method_history);
    interface_diagnostics.reinit_fallback_history = input_data.reinit_fallback_history;
    interface_diagnostics.reinit_reason_history = cellstr(input_data.reinit_reason_history);
    interface_diagnostics.hj_update_diagnostics_history = input_data.hj_update_diagnostics_history;
    interface_diagnostics.reinit_diagnostics_history = input_data.reinit_diagnostics_history;
    interface_diagnostics.boundary_guard_ratio_history = input_data.boundary_guard_ratio_history;
    interface_diagnostics.frozen_boundary_point_history = input_data.frozen_boundary_point_history;
    interface_diagnostics.accepted_hj_reinit_count = input_data.accepted_hj_reinit_count;
    interface_diagnostics.local_reinit_shell_size_history = input_data.local_reinit_shell_size_history;
    interface_diagnostics.post_reinit_grad_dev_mean_history = input_data.post_reinit_grad_dev_mean_history;
    interface_diagnostics.post_reinit_grad_outlier_ratio_history = input_data.post_reinit_grad_outlier_ratio_history;
    interface_diagnostics.refresh_shell_size_history = input_data.refresh_shell_size_history;
    interface_diagnostics.refresh_count = input_data.refresh_count;
    interface_diagnostics.last_hj_info = input_data.last_hj_info;
    interface_diagnostics.last_reinit_info = input_data.last_reinit_info;
    interface_diagnostics.last_velocity_field = input_data.last_velocity_field;
    interface_diagnostics.last_propagation_mask = input_data.last_propagation_mask;
    interface_diagnostics.boundary_guard_band = input_data.boundary_guard_band;
    interface_diagnostics.outside_velocity_nonzero_count = nnz(abs( ...
        input_data.last_velocity_field(~input_data.material_mask_full)) > 1e-14);
    interface_diagnostics.outside_propagation_nonzero_count = nnz(abs( ...
        input_data.last_velocity_field(~input_data.last_propagation_mask)) > 1e-14);
    interface_diagnostics.outside_phi_nonpositive_count = nnz( ...
        input_data.lsf(~input_data.material_mask_full) <= 0);
end
