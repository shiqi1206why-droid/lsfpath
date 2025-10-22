function save_checkpoint(iter, lsf, theta_e, compliance_history, FCS_history, params)
    % 保存优化检查点
    % 
    % 输入：
    %   iter - 当前迭代次数
    %   lsf - 水平集函数
    %   theta_e - 角度场
    %   compliance_history - 柔度历史
    %   FCS_history - FCS历史
    %   params - 参数结构
    
    % 检查是否启用检查点保存
    if ~params.debug.save_checkpoints
        return;
    end
    
    % 检查是否到达保存间隔
    if mod(iter, params.debug.checkpoint_interval) ~= 0
        return;
    end
    
    % 创建检查点目录
    checkpoint_dir = 'checkpoints';
    if ~exist(checkpoint_dir, 'dir')
        mkdir(checkpoint_dir);
    end
    
    % 生成文件名
    filename = sprintf('%s/checkpoint_iter_%04d.mat', checkpoint_dir, iter);
    
    % 保存数据
    try
        save(filename, 'iter', 'lsf', 'theta_e', 'compliance_history', ...
             'FCS_history', 'params', '-v7.3');
        
        fprintf('  [检查点] 已保存: %s\n', filename);
    catch ME
        warning('检查点保存失败: %s', ME.message);
    end
end

