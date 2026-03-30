function value = get_struct_field_or_default(data, field_name, default_value)
%GET_STRUCT_FIELD_OR_DEFAULT Return struct field value or a fallback default.

    value = default_value;
    if isstruct(data) && isfield(data, field_name) && ~isempty(data.(field_name))
        value = data.(field_name);
    end
end
