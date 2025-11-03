# Hamilton-Jacobi 更新导致梯度偏离的分析与修复建议

## 现象

在改进后的水平集悬臂梁优化流程中，几乎每一次 Hamilton-Jacobi (HJ) 更新之后，`|∇φ|` 都会显著偏离 1。即便随后使用 FMM 重初始化，也无法完全纠正，最终导致灵敏度传播异常、主路径漂移。

## 根因

问题出在 `core_computation/update_levelset_HJ.m` 中对梯度范数的离散化：

```matlab
grad_plus = sqrt( max(dmx,0)^2 + min(dpx,0)^2 + ...
                  max(dmy,0)^2 + min(dpy,0)^2 );
```

- 上述写法把正向、反向差分的平方 **直接累加**；当局部存在交替趋势时（如最大值附近），正向和反向项同时参与，实际等价于“把两个方向都算成了正贡献”，导致 `|∇φ|` 系统性放大。
- Godunov 型格式的正确做法应是对每个方向择取 **上风项的最大平方**：即在正速度下使用 `max(dmx,0)` 和 `max(dmy,0)`，负速度下使用 `max(dpx,0)` 和 `max(dpy,0)`；若反向差分发挥作用，则比较正向、反向平方，取较大的那一个，而不是二者求和。

## 修复建议

将 HJ 更新中的梯度计算改为 Godunov 形式。例如：

```matlab
if velocity(i,j) > 0
    a_plus  = max(dmx, 0);
    a_minus = min(dpx, 0);
    b_plus  = max(dmy, 0);
    b_minus = min(dpy, 0);
    grad = sqrt( max(a_plus^2, (-a_minus)^2) + ...
                 max(b_plus^2, (-b_minus)^2) );
else
    a_plus  = max(dpx, 0);
    a_minus = min(dmx, 0);
    b_plus  = max(dpy, 0);
    b_minus = min(dmy, 0);
    grad = sqrt( max(a_plus^2, (-a_minus)^2) + ...
                 max(b_plus^2, (-b_minus)^2) );
end

lsf_new(i,j) = lsf_old(i,j) - dt * velocity(i,j) * grad;
```

要点：

1. **逐方向选最大平方**：确保只由真正的上风项主导。
2. **区分正负速度**：正向传播与负向传播使用的差分不一样。
3. 如需兼容各向异性栅格，可在平方前乘以 `1/dx`、`1/dy` 等权重。

修正后，`|∇φ|` 将稳定在 1 附近，重初始化频率也会随之降低，整个灵敏度分析链条才能保持稳健。
