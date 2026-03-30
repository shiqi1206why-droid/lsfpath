function results = build_results_struct(input_data)
%BUILD_RESULTS_STRUCT Convert a flat input struct into the final results struct.

    results = struct();
    field_names = fieldnames(input_data);
    for field_idx = 1:numel(field_names)
        key = field_names{field_idx};
        results.(key) = input_data.(key);
    end
end
