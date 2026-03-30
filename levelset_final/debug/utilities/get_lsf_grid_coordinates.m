function [x_full, y_full, x_core, y_core] = get_lsf_grid_coordinates(lsf_or_size, dx, dy)
    % 返回与当前lsf存储方式一致的物理坐标

    if isnumeric(lsf_or_size) && isvector(lsf_or_size) && numel(lsf_or_size) == 2
        lsf_size = lsf_or_size;
    else
        lsf_size = size(lsf_or_size);
    end

    ny = lsf_size(1);
    nx = lsf_size(2);

    x_full = ((1:nx) - 1.5) * dx;
    y_full = ((1:ny) - 1.5) * dy;

    if nx >= 3
        x_core = x_full(2:end-1);
    else
        x_core = x_full;
    end
    if ny >= 3
        y_core = y_full(2:end-1);
    else
        y_core = y_full;
    end
end
