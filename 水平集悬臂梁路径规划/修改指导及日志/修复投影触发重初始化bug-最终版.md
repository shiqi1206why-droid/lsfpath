# 修复投影触发重初始化bug - 最终版

## 问题诊断

当前问题：投影修正量(≈0.03)持续超过重初始化阈值(0.015)，导致每步都触发重初始化。

根本原因：
- `lsf_change = max(abs(lsf - lsf_before))` 包含了HJ演化量 + 投影修正量
- 投影强度 omega=0.7 导致修正量≈0.03 > 阈值0.015
- 重初始化应该监控HJ演化导致的符号距离性质退化，而非投影修正

## 修改方案（已整合Codex建议）

### 1. 修改 `fiber_levelset.m` 的水平集更新与监控逻辑

**位置**: `fiber_levelset.m:413-455`

**完整修改后的代码结构**：

```matlab
% 3.6 更新水平集
lsf_before = lsf;
lsf = update_levelset_HJ(lsf, velocity, dt_adaptive, dx, dy);

% 保存HJ演化后状态（用于监控符号距离退化）
lsf_after_HJ = lsf;

% === 阶段3：水平集更新验证（检查NaN/Inf） ===
if any(~isfinite(lsf(:)))
    log_message('ERROR', params, '水平集出现NaN/Inf（迭代 %d）', iter);
    lsf = lsf_before;  % 回退到更新前状态
    lsf_after_HJ = lsf;  % 【关键修复】同步lsf_after_HJ，避免传NaN到should_reinitialize
    log_message('WARN', params, '已回退到更新前状态');
end

% === 强力稳住2：投影更硬（前期加强）===
% [投影代码保持不变，包括424-450行的所有投影逻辑]
if ENABLE_HARD_PROJECTION
    deviation_proj = lsf(proj_band) - lsf_target_global(proj_band);
    lsf(proj_band) = lsf(proj_band) - omega_proj * deviation_proj;
    
    if iter == 1 || mod(iter, 10) == 0
        max_correction = max(abs(omega_proj * deviation_proj));
        phase_str = '是';
        if iter > 100
            phase_str = '否';
        end
        fprintf('  [硬约束投影] omega=%.2f, 最大修正=%.4f (前期=%s)\n', ...
            omega_proj, max_correction, phase_str);
    end
end

% 分别计算HJ演化量和投影修正量
lsf_change_HJ = max(abs(lsf_after_HJ(:) - lsf_before(:)));
lsf_change_total = max(abs(lsf(:) - lsf_before(:)));
lsf_change_proj = max(abs(lsf(:) - lsf_after_HJ(:)));

% 诊断输出（每10步显示）
if mod(iter, 10) == 0
    fprintf('  [水平集变化] HJ演化=%.3e, 投影修正=%.3e, 总计=%.3e\n', ...
        lsf_change_HJ, lsf_change_proj, lsf_change_total);
end
```

**关键点**：
1. 在HJ更新后立即保存 `lsf_after_HJ`
2. **【Codex建议】** 在NaN回退时同步更新 `lsf_after_HJ = lsf`
3. 分别统计三种变化量：HJ演化、投影修正、总计

### 2. 修改重初始化判断传参

**位置**: `fiber_levelset.m:463-465`（原464行，位置可能因插入代码而变化）

**修改前**:
```matlab
iter_since_last_reinit = iter_since_last_reinit + 1;
[do_reinit, reinit_reason] = should_reinitialize(lsf, lsf_before, ...
    iter_since_last_reinit, iter, params);
```

**修改后**:
```matlab
iter_since_last_reinit = iter_since_last_reinit + 1;

% === 重初始化判断说明 ===
% 传入lsf_after_HJ（而非投影后的lsf）是刻意设计：
% 1. 重初始化目的：恢复|∇φ|≈1的符号距离性质
% 2. HJ演化会导致符号距离性质退化
% 3. 投影只是平移零等值线位置，不改变符号距离性质
% 4. 因此只监控HJ演化量，避免投影触发过度重初始化
% 注意：传入的lsf_after_HJ会用于梯度偏差计算，硬投影对梯度影响通常很小
[do_reinit, reinit_reason] = should_reinitialize(lsf_after_HJ, lsf_before, ...
    iter_since_last_reinit, iter, params);
```

**说明**：
- 第一个参数改为 `lsf_after_HJ`（HJ演化后、投影前）
- `should_reinitialize` 内部会计算 `lsf_change = max(abs(lsf_after_HJ - lsf_before))`
- 梯度偏差也基于 `lsf_after_HJ` 计算，通常投影对梯度影响很小

### 3. 可选：在重初始化触发时也输出变化量

**位置**: `fiber_levelset.m:467-470`（重初始化触发后）

**在重初始化日志中添加变化量信息**：
```matlab
if do_reinit
    log_message('INFO', params, '触发重初始化: %s', reinit_reason);
    fprintf('  [触发时状态] HJ演化=%.3e, 投影修正=%.3e\n', ...
        lsf_change_HJ, lsf_change_proj);
    zero_mask_dynamic = compute_zero_mask_from_lsf(lsf, h_grid);
    % [后续重初始化代码保持不变]
```

## 预期效果

修复后：
- ✅ 重初始化频率恢复正常（前期每5-10步，后期每10-15步）
- ✅ 诊断输出清晰区分HJ演化和投影修正的贡献
- ✅ 优化收敛性能提升（减少不必要的重初始化开销，约减少40-50次/500步）
- ✅ 逻辑更符合物理意义
- ✅ NaN回退时不会传递NaN到判据函数

## 验证方法

运行 `fiber_levelset('default')` 观察：

1. **重初始化频率**：应该恢复到正常模式
   - 前100步：每5-10步触发一次（约10-20次）
   - 后期：每10-15步触发一次（约25-40次）
   - 总计：约35-60次（vs 修复前500次）

2. **诊断输出检查**：
   ```
   [水平集变化] HJ演化=5.234e-03, 投影修正=3.081e-02, 总计=3.234e-02
   ```
   - HJ演化量应该 < 0.015（阈值）
   - 投影修正量约0.03（符合预期）
   - 总计应该≈HJ+投影

3. **优化质量指标**：
   - FCS（纤维连续性）≥85%（应保持或改善）
   - 柔度曲线单调下降
   - 无NaN/Inf错误

4. **重初始化触发原因分布**：
   - 前期主要是"前期固定频率"
   - 后期主要是"后期固定频率"或"达到最大间隔"
   - "水平集变化过大"应该很少出现（<5次）

## Codex评审要点总结

✅ **方案大方向正确**：把重初始化判据只针对HJ演化，不让投影掺进去

✅ **关键细节已修复**：NaN回退时同步`lsf_after_HJ = lsf`

⚠️ **后续观察点**：硬投影对梯度的影响（当前假设很小，若后续发现问题可补充统计）

## 实施顺序

1. 先修改水平集更新逻辑（保存lsf_after_HJ，NaN同步，三种变化量统计）
2. 再修改should_reinitialize调用（传参改为lsf_after_HJ）
3. 添加诊断输出
4. 运行测试验证
5. 观察500步完整优化结果

## 技术细节

### should_reinitialize 函数接口
```matlab
function [should_reinit, reason] = should_reinitialize(lsf, lsf_before, ...
    iter_since_last_reinit, iter, params)
```

函数内部两个判据：
1. **最大差值判据**：`lsf_change = max(abs(lsf(:) - lsf_before(:)))`
2. **梯度偏差判据**：基于传入的`lsf`计算`|∇φ|-1`的均值

修复方案：
- 传入 `lsf_after_HJ` 确保两个判据都只看HJ演化的影响
- 投影主要改变零等值线位置（平移），对梯度模影响通常很小

