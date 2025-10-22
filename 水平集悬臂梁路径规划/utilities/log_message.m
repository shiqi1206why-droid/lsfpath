function log_message(level, params, format, varargin)
    % 分级日志输出
    % 
    % 输入：
    %   level - 日志级别（'DEBUG', 'INFO', 'WARN', 'ERROR'）
    %   params - 参数结构（包含debug.log_level字段）
    %   format - printf格式字符串
    %   varargin - 格式化参数
    %
    % 日志级别优先级（从低到高）：
    %   DEBUG - 详细调试信息
    %   INFO  - 一般信息
    %   WARN  - 警告
    %   ERROR - 错误
    %
    % 示例：
    %   log_message('INFO', params, '正在初始化...');
    %   log_message('WARN', params, '参数 %s 超出范围', param_name);
    %   log_message('DEBUG', params, '中间值: %.3e', value);
    
    % 定义日志级别优先级
    level_priority = struct('DEBUG', 4, 'INFO', 3, 'WARN', 2, 'ERROR', 1);
    
    % 获取配置的日志级别
    if isfield(params, 'debug') && isfield(params.debug, 'log_level')
        config_level = params.debug.log_level;
    else
        config_level = 'INFO';  % 默认INFO级别
    end
    
    % 检查是否应该输出
    if ~isfield(level_priority, level)
        warning('未知日志级别: %s，使用INFO', level);
        level = 'INFO';
    end
    
    if ~isfield(level_priority, config_level)
        warning('配置的日志级别无效: %s，使用INFO', config_level);
        config_level = 'INFO';
    end
    
    config_priority = level_priority.(config_level);
    msg_priority = level_priority.(level);
    
    % 只输出优先级足够高的消息
    if msg_priority <= config_priority
        % 添加级别前缀
        prefix = sprintf('[%5s]', level);
        
        % 格式化并输出消息
        if nargin > 3
            message = sprintf(format, varargin{:});
        else
            message = format;
        end
        
        fprintf('%s %s\n', prefix, message);
        
        % ERROR级别额外使用warning函数（生成堆栈跟踪）
        if strcmp(level, 'ERROR')
            warning('off', 'backtrace');  % 关闭重复的堆栈信息
            warning('Error: %s', message);
            warning('on', 'backtrace');
        end
    end
end

