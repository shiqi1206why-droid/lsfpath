function project_root = get_project_root(anchor_path)
%GET_PROJECT_ROOT Resolve the maintained debug project root.

    if nargin < 1 || isempty(anchor_path)
        anchor_path = mfilename('fullpath');
    end

    if isfolder(anchor_path)
        search_dir = char(anchor_path);
    else
        search_dir = fileparts(char(anchor_path));
    end

    while true
        is_project_root = exist(fullfile(search_dir, 'fiber_levelset.m'), 'file') == 2 && ...
            exist(fullfile(search_dir, 'config'), 'dir') == 7 && ...
            exist(fullfile(search_dir, 'utilities'), 'dir') == 7;
        if is_project_root
            project_root = search_dir;
            return;
        end

        parent_dir = fileparts(search_dir);
        if isempty(parent_dir) || strcmp(parent_dir, search_dir)
            error('get_project_root:NotFound', ...
                'Unable to resolve project root from anchor: %s', char(anchor_path));
        end
        search_dir = parent_dir;
    end
end
