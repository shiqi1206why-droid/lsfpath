function [x_points, y_points] = contourc_to_points(C)
    % 将 contourc 的输出转换为坐标点序列
    % 多段等值线之间插入NaN分隔，避免跨段错误连线

    x_points = [];
    y_points = [];
    if isempty(C)
        x_points = x_points(:);
        y_points = y_points(:);
        return;
    end

    idx = 1;
    first_segment = true;
    while idx < size(C, 2)
        count = C(2, idx);
        segment = C(:, idx+1:idx+count);
        
        % 如果不是第一段，插入NaN分隔
        if ~first_segment
            x_points = [x_points, NaN];
            y_points = [y_points, NaN];
        end
        first_segment = false;
        
        x_points = [x_points, segment(1, :)];
        y_points = [y_points, segment(2, :)];
        idx = idx + count + 1;
    end
    x_points = x_points(:);
    y_points = y_points(:);
end

