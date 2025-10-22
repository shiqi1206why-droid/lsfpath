function [iter, lsf, theta_e, compliance_history, FCS_history, params] = ...
    load_checkpoint(checkpoint_file)
    % 从检查点恢复优化
    % 
    % 输入：
    %   checkpoint_file - 检查点文件路径
    % 
    % 输出：
    %   iter - 保存时的迭代次数
    %   lsf - 水平集函数
    %   theta_e - 角度场
    %   compliance_history - 柔度历史
    %   FCS_history - FCS历史
    %   params - 参数结构
    %
    % 示例：
    %   [iter, lsf, theta_e, ch, fh, params] = ...
    %       load_checkpoint('checkpoints/checkpoint_iter_0050.mat');
    
    if ~exist(checkpoint_file, 'file')
        error('检查点文件不存在: %s', checkpoint_file);
    end
    
    try
        data = load(checkpoint_file);
        iter = data.iter;
        lsf = data.lsf;
        theta_e = data.theta_e;
        compliance_history = data.compliance_history;
        FCS_history = data.FCS_history;
        params = data.params;
        
        fprintf('✅ 从检查点恢复: iter=%d, 文件=%s\n', iter, checkpoint_file);
    catch ME
        error('检查点加载失败: %s', ME.message);
    end
end

