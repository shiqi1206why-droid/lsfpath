function paths = build_project_paths(project_root)
%BUILD_PROJECT_PATHS Build canonical absolute paths for the debug project.

    if nargin < 1 || isempty(project_root)
        project_root = get_project_root();
    end

    paths = struct();
    paths.project_root = project_root;
    paths.topology_file = fullfile(project_root, 'topo_result.mat');
    paths.checkpoint_dir = fullfile(project_root, 'checkpoints');
    paths.visualization_dir = fullfile(project_root, 'visualization_artifacts');
    paths.baseline_dir = fullfile(project_root, 'baseline_artifacts');
    paths.refactor_dir = fullfile(project_root, 'refactor_artifacts');
    paths.tests_dir = fullfile(project_root, 'tests');
end
