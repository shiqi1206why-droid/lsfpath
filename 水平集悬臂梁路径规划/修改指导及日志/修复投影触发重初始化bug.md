# 修复投影触发重初始化bug

## 问题诊断

当前问题：投影修正量(≈0.03)持续超过重初始化阈值(0.015)，导致每步都触发重初始化。

根本原因：
- `lsf_change = max(abs(lsf - lsf_before))` 包含了HJ演化量 + 投影修正量
- 投影强度 omega=0.7 导致修正量≈0.03 > 阈值0.015
- 重初始化应该监控HJ演化导致的符号距离性质退化，而非投影修正

## 修改方案

### 1. 修改 `fiber_levelset.m` 的水平集更新与监控逻辑

**位置**: `fiber_levelset.m:413-455`

**修改内容**:

在投影之前保存HJ演化后的状态：
```matlab
% 3.6 更新水平集
lsf_before = lsf;
lsf = update_levelset_HJ(lsf, velocity, dt_adaptive, dx, dy);

% 保存HJ演化后状态（用于监控符号距离退化）
lsf_after_HJ = lsf;

% [NaN/Inf检查保持不变]

% [投影代码保持不变]

% 分别计算HJ演化量和投影修正量
lsf_change_HJ = max(abs(lsf_after_HJ(:) - lsf_before(:)));
lsf_change_total = max(abs(lsf(:) - lsf_before(:)));
lsf_change_proj = max(abs(lsf(:) - lsf_after_HJ(:)));
```

添加诊断输出（每10步或重初始化时）：
```matlab
if mod(iter, 10) == 0 || do_reinit
    fprintf('  [水平集变化] HJ演化=%.3e, 投影修正=%.3e, 总计=%.3e\n', ...
        lsf_change_HJ, lsf_change_proj, lsf_change_total);
end
```

### 2. 修改重初始化判断传参

**位置**: `fiber_levelset.m:464`

**修改前**:
```matlab
[do_reinit, reinit_reason] = should_reinitialize(lsf, lsf_before, ...
    iter_since_last_reinit, iter, params);
```

**修改后**:
```matlab
% 传入HJ演化后的lsf，而不是投影后的lsf
[do_reinit, reinit_reason] = should_reinitialize(lsf_after_HJ, lsf_before, ...
    iter_since_last_reinit, iter, params);
```

**说明**：传入`lsf_after_HJ`（HJ演化后但投影前的状态）作为第一个参数，使得`should_reinitialize`内部计算`lsf_change = max(abs(lsf_after_HJ - lsf_before))`时，只监控HJ演化导致的变化，排除投影修正的影响。

### 3. 在重初始化逻辑中添加说明注释

**位置**: `fiber_levelset.m:467-485`（修改后位置可能变化）

在重初始化代码块前添加注释说明为何只监控HJ演化：
```matlab
% === 重初始化判断说明 ===
% 传入lsf_after_HJ（而非投影后的lsf）是刻意设计：
% 1. 重初始化目的：恢复|∇φ|≈1的符号距离性质
% 2. HJ演化会导致符号距离性质退化
% 3. 投影只是平移零等值线位置，不改变符号距离性质
% 4. 因此只监控HJ演化量，避免投影触发过度重初始化
```

## 预期效果

修复后：
- 重初始化频率恢复正常（前期每5-10步，后期每10-15步）
- 诊断输出清晰区分HJ演化和投影修正的贡献
- 优化收敛性能可能提升（减少不必要的重初始化开销）
- 逻辑更符合物理意义

## 验证方法

运行 `fiber_levelset('default')` 观察：
1. 重初始化触发频率是否降低
2. 诊断输出中HJ演化量是否 < 0.015h
3. 优化质量指标（FCS、柔度）是否保持或改善

## 技术细节补充

### should_reinitialize 函数接口
```matlab
function [should_reinit, reason] = should_reinitialize(lsf, lsf_before, ...
    iter_since_last_reinit, iter, params)
```

该函数内部计算：
```matlab
lsf_change = max(abs(lsf(:) - lsf_before(:)));
threshold_change = params.levelset.reinit_threshold * h;
if lsf_change > threshold_change
    should_reinit = true;
end
```

修复方案通过传入`lsf_after_HJ`而非投影后的`lsf`，确保`lsf_change`只反映HJ演化导致的变化。

