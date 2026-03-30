function paths = resolve_runtime_paths(context)
%RESOLVE_RUNTIME_PATHS Resolve canonical project paths from params or fallback.

    if nargin >= 1 && isstruct(context)
        if isfield(context, 'runtime') && isstruct(context.runtime) && ...
                isfield(context.runtime, 'paths') && ~isempty(context.runtime.paths)
            paths = context.runtime.paths;
            return;
        end
        if isfield(context, 'paths') && isstruct(context.paths) && ...
                isfield(context.paths, 'project_root')
            paths = context.paths;
            return;
        end
        if isfield(context, 'project_root') && ~isempty(context.project_root)
            paths = build_project_paths(context.project_root);
            return;
        end
    end

    paths = build_project_paths(get_project_root());
end
