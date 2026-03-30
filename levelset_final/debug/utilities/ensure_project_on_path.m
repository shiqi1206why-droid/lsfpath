function cleanup = ensure_project_on_path(project_root)
%ENSURE_PROJECT_ON_PATH Add the debug project tree to the MATLAB path.

    if nargin < 1 || isempty(project_root)
        project_root = get_project_root();
    end

    project_path = genpath(project_root);
    addpath(project_path, '-begin');
    cleanup = onCleanup(@() rmpath(project_path));
end
