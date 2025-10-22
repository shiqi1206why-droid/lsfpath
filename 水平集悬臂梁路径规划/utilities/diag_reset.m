function diag_reset()
    % 重置全局诊断计数器

    global DIAG;
    fields = {'theta_grad_samples','theta_grad_small','den_small','contrib_total','contrib_nonzero'};
    for k = 1:numel(fields)
        DIAG.(fields{k}) = 0;
    end
end

