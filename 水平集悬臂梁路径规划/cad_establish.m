function cad_establish(lsf, dx, dy, nelx, nely, output_filename, options)
% CAD_ESTABLISH 将优化后的纤维路径导出为DXF格式CAD文件
%
% 输入参数：
%   lsf              - 水平集函数 [nely+2, nelx+2]（包含ghost cells）
%   dx, dy           - 网格物理间距（米）
%   nelx, nely       - 网格单元数量
%   output_filename  - 输出DXF文件名（不含扩展名）
%   options          - 可选参数结构体（可选）
%                      .levels: 等值线层级数组（默认：[-2h, -h, 0, h, 2h]）
%                      .layer_prefix: DXF图层名前缀（默认：'Fiber_Path'）
%
% 输出：
%   生成 [output_filename].dxf 文件
%
% 示例：
%   cad_establish(lsf_final, 0.02, 0.02, 80, 50, 'fiber_paths');
%
% 作者：基于fiber_levelset.m项目
% 日期：2025-10-23

    %% 参数验证
    if nargin < 6
        error('至少需要6个输入参数：lsf, dx, dy, nelx, nely, output_filename');
    end
    
    if nargin < 7
        options = struct();
    end
    
    % 检查lsf尺寸
    if ~all(size(lsf) == [nely+2, nelx+2])
        error('lsf尺寸应为 [nely+2=%d, nelx+2=%d]，实际为 [%d, %d]', ...
            nely+2, nelx+2, size(lsf, 1), size(lsf, 2));
    end
    
    %% 设置默认参数
    h = min(dx, dy);  % 最小网格间距
    
    if ~isfield(options, 'levels')
        % 默认输出5条等值线：零等值线及其上下各2层
        options.levels = [-2*h, -h, 0, h, 2*h];
    end
    
    if ~isfield(options, 'layer_prefix')
        options.layer_prefix = 'Fiber_Path';
    end
    
    %% 提取等值线数据
    fprintf('正在提取纤维路径等值线...\n');
    
    % 去除ghost cells，获取核心区域
    lsf_core = lsf(2:end-1, 2:end-1);  % 尺寸: nely × nelx
    
    % 构建物理坐标网格（使用单元中心坐标）
    % lsf_core是nely × nelx，需要对应长度的坐标向量
    x_coords = dx * (0.5:1:nelx-0.5);  % 长度: nelx
    y_coords = dy * (0.5:1:nely-0.5);  % 长度: nely
    
    % 存储所有等值线
    all_contours = struct('level', {}, 'lines', {}, 'layer_name', {});
    
    for i = 1:length(options.levels)
        level = options.levels(i);
        
        % 使用contourc提取等值线（x_coords和y_coords是向量）
        C = contourc(x_coords, y_coords, lsf_core, [level level]);
        
        if isempty(C)
            fprintf('  警告：水平集层级 %.4f 无等值线\n', level);
            continue;
        end
        
        % 解析contourc返回的数据
        lines = parse_contourc_matrix(C);
        
        if isempty(lines)
            continue;
        end
        
        % 图层名称
        if level == 0
            layer_name = [options.layer_prefix '_Zero'];
        elseif level > 0
            layer_name = sprintf('%s_Plus%d', options.layer_prefix, i);
        else
            layer_name = sprintf('%s_Minus%d', options.layer_prefix, i);
        end
        
        % 存储
        contour_data = struct();
        contour_data.level = level;
        contour_data.lines = lines;
        contour_data.layer_name = layer_name;
        all_contours(end+1) = contour_data; %#ok<AGROW>
        
        fprintf('  层级 %+.4f: 提取 %d 条路径\n', level, length(lines));
    end
    
    if isempty(all_contours)
        error('未能提取到任何等值线，请检查水平集函数');
    end
    
    %% 写入DXF文件
    fprintf('正在写入DXF文件: %s.dxf\n', output_filename);
    
    dxf_filename = [output_filename '.dxf'];
    fid = fopen(dxf_filename, 'w');
    
    if fid == -1
        error('无法创建文件: %s', dxf_filename);
    end
    
    try
        % 写入各部分（传递实际坐标范围）
        x_min = min(x_coords);
        x_max = max(x_coords);
        y_min = min(y_coords);
        y_max = max(y_coords);
        write_dxf_header(fid, x_min, x_max, y_min, y_max);
        write_dxf_tables(fid, all_contours);
        write_dxf_entities(fid, all_contours);
        write_dxf_eof(fid);
        
        fclose(fid);
        fprintf('✓ 成功生成CAD文件: %s\n', dxf_filename);
        fprintf('  总路径数: %d\n', sum(arrayfun(@(c) length(c.lines), all_contours)));
        
    catch ME
        fclose(fid);
        error('写入DXF文件失败: %s', ME.message);
    end
end

%% ========== 辅助函数 ==========

function lines = parse_contourc_matrix(C)
    % 解析MATLAB contourc返回的矩阵格式
    % 输入：C - contourc返回的2×N矩阵
    % 输出：lines - 结构体数组，每个元素包含x, y坐标数组
    
    lines = struct('x', {}, 'y', {});
    
    if isempty(C)
        return;
    end
    
    idx = 1;
    while idx < size(C, 2)
        % C(:,idx) = [level; num_points]
        num_points = C(2, idx);
        
        if idx + num_points > size(C, 2)
            break;  % 数据不完整
        end
        
        % 提取坐标
        x_coords = C(1, idx+1:idx+num_points);
        y_coords = C(2, idx+1:idx+num_points);
        
        % 存储
        line = struct();
        line.x = x_coords;
        line.y = y_coords;
        lines(end+1) = line; %#ok<AGROW>
        
        % 移动到下一条曲线
        idx = idx + num_points + 1;
    end
end

function write_dxf_header(fid, x_min, x_max, y_min, y_max)
    % 写入DXF文件头部
    
    fprintf(fid, '0\nSECTION\n');
    fprintf(fid, '2\nHEADER\n');
    
    % AutoCAD版本 (R12 = AC1009，兼容性更好)
    fprintf(fid, '9\n$ACADVER\n');
    fprintf(fid, '1\nAC1009\n');
    
    % 图纸范围
    fprintf(fid, '9\n$EXTMIN\n');
    fprintf(fid, '10\n%.6f\n20\n%.6f\n30\n0.0\n', x_min, y_min);
    
    fprintf(fid, '9\n$EXTMAX\n');
    fprintf(fid, '10\n%.6f\n20\n%.6f\n30\n0.0\n', x_max, y_max);
    
    % 单位（米）
    fprintf(fid, '9\n$INSUNITS\n');
    fprintf(fid, '70\n6\n');  % 6 = 米
    
    fprintf(fid, '0\nENDSEC\n');
end

function write_dxf_tables(fid, all_contours)
    % 写入DXF图层定义（R12格式，含LTYPE表和Layer 0）
    
    fprintf(fid, '0\nSECTION\n');
    fprintf(fid, '2\nTABLES\n');
    
    % === 线型表（LTYPE）===
    fprintf(fid, '0\nTABLE\n');
    fprintf(fid, '2\nLTYPE\n');
    fprintf(fid, '70\n1\n');  % 1个线型
    
    % 定义CONTINUOUS线型（规范要求）
    fprintf(fid, '0\nLTYPE\n');
    fprintf(fid, '2\nCONTINUOUS\n');        % 线型名
    fprintf(fid, '70\n0\n');                 % 标准标志
    fprintf(fid, '3\nSolid line\n');         % 描述
    fprintf(fid, '72\n65\n');                % 对齐代码
    fprintf(fid, '73\n0\n');                 % 划线段数
    fprintf(fid, '40\n0.0\n');               % 总长度
    
    fprintf(fid, '0\nENDTAB\n');
    
    % === 图层表（LAYER）===
    fprintf(fid, '0\nTABLE\n');
    fprintf(fid, '2\nLAYER\n');
    fprintf(fid, '70\n%d\n', length(all_contours) + 1);  % 图层数（含Layer 0）
    
    % Layer 0（默认图层，规范要求）
    fprintf(fid, '0\nLAYER\n');
    fprintf(fid, '2\n0\n');                  % 图层名 = 0
    fprintf(fid, '70\n0\n');                 % 标准标志
    fprintf(fid, '62\n7\n');                 % 颜色 = 白色
    fprintf(fid, '6\nCONTINUOUS\n');         % 线型
    
    % 为每个等值线层级定义图层
    for i = 1:length(all_contours)
        layer_name = all_contours(i).layer_name;
        color_index = mod(i-1, 7) + 1;  % 颜色1-7循环
        
        fprintf(fid, '0\nLAYER\n');
        fprintf(fid, '2\n%s\n', layer_name);  % 图层名
        fprintf(fid, '70\n0\n');              % 标志
        fprintf(fid, '62\n%d\n', color_index); % 颜色
        fprintf(fid, '6\nCONTINUOUS\n');      % 线型
    end
    
    fprintf(fid, '0\nENDTAB\n');
    fprintf(fid, '0\nENDSEC\n');
end

function write_dxf_entities(fid, all_contours)
    % 写入DXF几何实体（R12格式：POLYLINE/VERTEX/SEQEND）
    
    fprintf(fid, '0\nSECTION\n');
    fprintf(fid, '2\nENTITIES\n');
    
    % 遍历所有等值线
    for i = 1:length(all_contours)
        layer_name = all_contours(i).layer_name;
        lines = all_contours(i).lines;
        
        % 遍历该层级的所有路径
        for j = 1:length(lines)
            x_coords = lines(j).x;
            y_coords = lines(j).y;
            
            % === POLYLINE 实体头（R12格式）===
            fprintf(fid, '0\nPOLYLINE\n');
            fprintf(fid, '8\n%s\n', layer_name);     % 图层
            fprintf(fid, '66\n1\n');                 % 后续顶点标志（1=有VERTEX）
            fprintf(fid, '70\n0\n');                 % 标志（0=开放，1=闭合）
            fprintf(fid, '10\n0.0\n');               % 默认X（R12要求）
            fprintf(fid, '20\n0.0\n');               % 默认Y
            fprintf(fid, '30\n0.0\n');               % 默认Z
            
            % === VERTEX 实体（每个顶点）===
            for k = 1:length(x_coords)
                fprintf(fid, '0\nVERTEX\n');
                fprintf(fid, '8\n%s\n', layer_name); % 图层
                fprintf(fid, '10\n%.6f\n', x_coords(k));  % X坐标
                fprintf(fid, '20\n%.6f\n', y_coords(k));  % Y坐标
                fprintf(fid, '30\n0.0\n');                % Z坐标
                fprintf(fid, '70\n0\n');                  % 标志
            end
            
            % === SEQEND 结束标记（R12要求）===
            fprintf(fid, '0\nSEQEND\n');
            fprintf(fid, '8\n%s\n', layer_name);     % 图层
        end
    end
    
    fprintf(fid, '0\nENDSEC\n');
end

function write_dxf_eof(fid)
    % 写入DXF文件结束标记
    fprintf(fid, '0\nEOF\n');
end

