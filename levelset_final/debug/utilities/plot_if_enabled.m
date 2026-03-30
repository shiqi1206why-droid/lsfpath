function plot_if_enabled(params, plot_function)
    % 可选可视化控制
    % 
    % 输入：
    %   params - 参数结构（包含debug.enable_plots字段）
    %   plot_function - 绘图函数句柄
    %
    % 示例：
    %   plot_if_enabled(params, @() plot_results(data));
    %   plot_if_enabled(params, @() my_plot_function(x, y, z));
    
    if isfield(params, 'debug') && isfield(params.debug, 'enable_plots')
        enable = params.debug.enable_plots;
    else
        enable = true;  % 默认启用
    end
    
    if enable
        try
            plot_function();
        catch ME
            warning(ME.identifier, '绘图失败: %s', ME.message);
        end
    end
end

