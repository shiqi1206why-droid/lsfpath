function lsf = fmm_reinitialize(lsf, dx, dy, zero_mask, material_mask)
% FMM_REINITIALIZE 使用快速行进法（Fast Marching Method）重新初始化水平集函数
%
% 输入：
%   lsf           - 当前水平集函数（包含ghost cells，尺寸 (nely+2) × (nelx+2)）
%   dx, dy        - 网格单元尺寸（可不等）
%   zero_mask     - 零水平集种子区域的逻辑掩膜（与lsf同尺寸）
%   material_mask - 材料区域掩膜（可选，默认为全域可达）
%
% 输出：
%   lsf           - 重新初始化后的符号距离函数，满足 |∇φ| ≈ 1
%
% 方法：
%   基于Eikonal方程 |∇T| = 1/F_s，其中F_s≡1（纯几何距离）
%   1) 从zero_mask构造初始Alive集合（φ=0的种子点）
%   2) 使用FMM单调推进计算到达时间T
%   3) 用原lsf的符号恢复：φ = sgn(φ_old) · T
%   4) 维持ghost cells边界条件
%
% 参考：接入FMM实现流程.txt 第1-163行

    % 处理可选参数
    if nargin < 5 || isempty(material_mask)
        material_mask = true(size(lsf));
    end
    
    % 保存原始符号
    sign_field = sign(lsf);
    sign_field(sign_field == 0) = 1;  % 零处默认为正
    
    % 提取核心区域（去除ghost cells）
    [ny_total, nx_total] = size(lsf);
    nely = ny_total - 2;
    nelx = nx_total - 2;
    
    % 确保material_mask与lsf尺寸一致
    if all(size(material_mask) == [nely, nelx])
        % material_mask是核心区域，需要扩展到包含ghost cells
        material_mask_full = false(ny_total, nx_total);
        material_mask_full(2:end-1, 2:end-1) = material_mask;
        % 复制边界
        material_mask_full(1, :) = material_mask_full(2, :);
        material_mask_full(end, :) = material_mask_full(end-1, :);
        material_mask_full(:, 1) = material_mask_full(:, 2);
        material_mask_full(:, end) = material_mask_full(:, end-1);
        material_mask = material_mask_full;
    end
    
    % 执行FMM计算距离场
    T = fmm_compute_distance(zero_mask, dx, dy, material_mask);
    
    % === FMM结果有限化：防止Inf传播到下游计算 ===
    % 参考：代码修改交流.md (Codex方案 2025-10-17)
    T_finite = T(isfinite(T));
    if isempty(T_finite)
        % 极端兜底：全域为0，随后由零水平集恢复符号
        T(:) = 0;
    else
        % 将不可达区域（T=Inf）替换为有限大值
        Tmax = max(T_finite);
        pad = 10 * max(dx, dy);  % 远离零水平集的"标记距离"
        T(~isfinite(T)) = Tmax + pad;
    end
    
    % 恢复符号：φ = sgn(φ_old) · T
    lsf = sign_field .* T;
    
    % 确保零水平集种子点保持为零
    lsf(zero_mask) = 0;
    
    % 维持ghost cells边界条件（Neumann边界）
    lsf(1, :) = lsf(2, :);
    lsf(end, :) = lsf(end-1, :);
    lsf(:, 1) = lsf(:, 2);
    lsf(:, end) = lsf(:, end-1);
end

%% ================== FMM核心算法 ==================

function T = fmm_compute_distance(zero_mask, dx, dy, material_mask)
% FMM_COMPUTE_DISTANCE 使用快速行进法计算到达时间场（距离场）
%
% 输入：
%   zero_mask     - 种子点集合（Γ₀）的逻辑掩膜
%   dx, dy        - 网格间距
%   material_mask - 可达区域掩膜
%
% 输出：
%   T             - 到达时间场（几何距离），不可达点为 +∞

    [ny, nx] = size(zero_mask);
    
    % 初始化到达时间为无穷大
    T = inf(ny, nx);
    
    % 定义状态标签：0=Far, 1=Trial, 2=Alive
    State = zeros(ny, nx, 'uint8');
    
    % 速度场（纯几何距离，F_s ≡ 1）
    F_s = 1.0;
    
    % 初始化：将种子点设为Alive，T=0
    seed_points = find(zero_mask & material_mask);
    T(seed_points) = 0;
    State(seed_points) = 2;  % Alive
    
    % 创建最小堆（优先队列）
    heap = [];  % 每个元素: [T_value, linear_index]
    
    % 将种子点的邻居加入Trial集合
    for idx = 1:length(seed_points)
        sp = seed_points(idx);
        [i, j] = ind2sub([ny, nx], sp);
        
        % 检查四邻域
        neighbors = get_neighbors(i, j, ny, nx);
        for n = 1:size(neighbors, 1)
            ni = neighbors(n, 1);
            nj = neighbors(n, 2);
            nidx = sub2ind([ny, nx], ni, nj);
            
            if State(ni, nj) == 0 && material_mask(ni, nj)  % Far状态且可达
                % 计算初始到达时间
                T_new = compute_local_update(T, ni, nj, dx, dy, F_s, State);
                T(ni, nj) = T_new;
                State(ni, nj) = 1;  % Trial
                heap = heap_push(heap, T_new, nidx);
            end
        end
    end
    
    % FMM主循环：从堆中取最小Trial点，推进前沿
    while ~isempty(heap)
        % Pop最小的Trial点
        [heap, T_min, idx_min] = heap_pop(heap);
        [i, j] = ind2sub([ny, nx], idx_min);
        
        % 标记为Alive
        State(i, j) = 2;
        
        % 更新其邻居
        neighbors = get_neighbors(i, j, ny, nx);
        for n = 1:size(neighbors, 1)
            ni = neighbors(n, 1);
            nj = neighbors(n, 2);
            nidx = sub2ind([ny, nx], ni, nj);
            
            if State(ni, nj) ~= 2 && material_mask(ni, nj)  % 非Alive且可达
                % 计算新的到达时间
                T_new = compute_local_update(T, ni, nj, dx, dy, F_s, State);
                
                if State(ni, nj) == 0  % Far -> Trial
                    T(ni, nj) = T_new;
                    State(ni, nj) = 1;
                    heap = heap_push(heap, T_new, nidx);
                elseif T_new < T(ni, nj)  % Trial，更新值
                    T(ni, nj) = T_new;
                    heap = heap_decrease_key(heap, T_new, nidx);
                end
            end
        end
    end
    
    % 不可达点保持为 +∞
end

%% ================== 局部更新公式 ==================

function T_new = compute_local_update(T, i, j, dx, dy, F_s, State)
% COMPUTE_LOCAL_UPDATE 使用上风格式计算网格点(i,j)的新到达时间
%
% 基于Eikonal方程的离散化：
%   对于各向异性网格（dx≠dy），使用二次方程求解
%   参考：接入FMM实现流程.txt 第44-74行

    [ny, nx] = size(T);
    
    % 获取四邻居的到达时间（仅使用Alive的邻居）
    T_L = inf; T_R = inf; T_D = inf; T_U = inf;
    
    if j > 1 && State(i, j-1) == 2
        T_L = T(i, j-1);
    end
    if j < nx && State(i, j+1) == 2
        T_R = T(i, j+1);
    end
    if i > 1 && State(i-1, j) == 2
        T_D = T(i-1, j);
    end
    if i < ny && State(i+1, j) == 2
        T_U = T(i+1, j);
    end
    
    % 取每个方向的最小值（上风格式）
    a = min(T_L, T_R);  % x方向
    b = min(T_D, T_U);  % y方向
    
    % 根据可用方向数量选择更新方式
    if isinf(a) && isinf(b)
        % 无可用邻居
        T_new = inf;
    elseif isinf(a)
        % 仅y方向可用
        T_new = b + dy / F_s;
    elseif isinf(b)
        % 仅x方向可用
        T_new = a + dx / F_s;
    else
        % 两个方向都可用，使用二元更新公式
        % 解方程：((T-a)/dx)^2 + ((T-b)/dy)^2 = 1/F_s^2
        % 转化为：A*T^2 + B*T + C = 0
        
        A = 1/dx^2 + 1/dy^2;
        B = -2 * (a/dx^2 + b/dy^2);
        C = a^2/dx^2 + b^2/dy^2 - 1/F_s^2;
        
        discriminant = B^2 - 4*A*C;
        
        if discriminant >= 0
            % 取较大根（保证T >= max(a,b)）
            T_candidate = (-B + sqrt(discriminant)) / (2*A);
            
            % 验证单调性：T必须 >= max(a,b)
            if T_candidate >= max(a, b)
                T_new = T_candidate;
            else
                % 退化为一次更新
                T_new = min(a + dx/F_s, b + dy/F_s);
            end
        else
            % 判别式为负，退化为一次更新
            T_new = min(a + dx/F_s, b + dy/F_s);
        end
    end
end

%% ================== 辅助函数 ==================

function neighbors = get_neighbors(i, j, ny, nx)
% GET_NEIGHBORS 获取网格点(i,j)的四邻域坐标
%
% 输出：N×2矩阵，每行为[i, j]

    neighbors = [];
    
    if i > 1
        neighbors = [neighbors; i-1, j];  % 下
    end
    if i < ny
        neighbors = [neighbors; i+1, j];  % 上
    end
    if j > 1
        neighbors = [neighbors; i, j-1];  % 左
    end
    if j < nx
        neighbors = [neighbors; i, j+1];  % 右
    end
end

%% ================== 最小堆实现 ==================

function heap = heap_push(heap, value, index)
% HEAP_PUSH 向最小堆中插入新元素
%
% 输入：
%   heap  - 当前堆（N×2矩阵，第1列为值，第2列为索引）
%   value - 插入的值
%   index - 对应的线性索引
%
% 输出：
%   heap  - 更新后的堆

    heap = [heap; value, index];
    
    % 上浮操作
    pos = size(heap, 1);
    while pos > 1
        parent = floor(pos / 2);
        if heap(pos, 1) < heap(parent, 1)
            % 交换
            temp = heap(pos, :);
            heap(pos, :) = heap(parent, :);
            heap(parent, :) = temp;
            pos = parent;
        else
            break;
        end
    end
end

function [heap, min_value, min_index] = heap_pop(heap)
% HEAP_POP 从最小堆中取出最小元素
%
% 输出：
%   heap      - 更新后的堆
%   min_value - 最小值
%   min_index - 最小值对应的索引

    if isempty(heap)
        min_value = inf;
        min_index = 0;
        return;
    end
    
    min_value = heap(1, 1);
    min_index = heap(1, 2);
    
    % 将最后一个元素移到堆顶
    heap(1, :) = heap(end, :);
    heap(end, :) = [];
    
    if isempty(heap)
        return;
    end
    
    % 下沉操作
    pos = 1;
    n = size(heap, 1);
    
    while true
        left = 2 * pos;
        right = 2 * pos + 1;
        smallest = pos;
        
        if left <= n && heap(left, 1) < heap(smallest, 1)
            smallest = left;
        end
        if right <= n && heap(right, 1) < heap(smallest, 1)
            smallest = right;
        end
        
        if smallest ~= pos
            % 交换
            temp = heap(pos, :);
            heap(pos, :) = heap(smallest, :);
            heap(smallest, :) = temp;
            pos = smallest;
        else
            break;
        end
    end
end

function heap = heap_decrease_key(heap, new_value, target_index)
% HEAP_DECREASE_KEY 减小堆中某个元素的值
%
% 输入：
%   heap         - 当前堆
%   new_value    - 新的更小的值
%   target_index - 要更新的元素的线性索引
%
% 输出：
%   heap         - 更新后的堆

    % 查找目标元素在堆中的位置
    pos = 0;
    for k = 1:size(heap, 1)
        if heap(k, 2) == target_index
            pos = k;
            break;
        end
    end
    
    if pos == 0
        % 未找到，直接返回
        return;
    end
    
    % 更新值
    heap(pos, 1) = new_value;
    
    % 上浮操作（因为值减小了）
    while pos > 1
        parent = floor(pos / 2);
        if heap(pos, 1) < heap(parent, 1)
            % 交换
            temp = heap(pos, :);
            heap(pos, :) = heap(parent, :);
            heap(parent, :) = temp;
            pos = parent;
        else
            break;
        end
    end
end

