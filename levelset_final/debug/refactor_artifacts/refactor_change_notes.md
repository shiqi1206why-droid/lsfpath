# 重构修改说明（按模块）

## A. 通用工具模块（`utilities`）

### 1) `normalize_mask_to_lsf_grid.m`（新增）
- 修改目的：统一“核心网格掩膜/含 ghost 掩膜”尺寸归一逻辑，避免多文件重复实现。
- 影响范围：仅掩膜尺寸转换与错误信息一致化，不改数值算法。

### 2) `apply_neumann_boundary.m`（新增）
- 修改目的：统一 ghost ring Neumann 边界复制逻辑，减少重复代码。
- 影响范围：边界复制行为保持原实现。

## B. 核心算法调用层

### 3) `core_computation/update_levelset_HJ.m`
- 修改目的：移除本地重复 helper，改为调用 `utilities` 公共函数。
- 具体变化：
  - `active_mask`、`stencil_mask` 统一走 `normalize_mask_to_lsf_grid`
  - 边界处理统一走 `apply_neumann_boundary`
- 逻辑变化：无（仅去重复）。

### 4) `level_set_evolution/build_velocity_field.m`
- 修改目的：掩膜归一逻辑去重。
- 具体变化：`propagation_mask` 改为调用公共掩膜归一函数。
- 逻辑变化：无。

### 5) `fmm_reinitialize.m`
- 修改目的：掩膜归一与 Neumann 边界逻辑去重。
- 具体变化：
  - `material_mask`、`zero_mask`、`local_shell_mask` 统一调用公共掩膜归一
  - 边界复制调用公共 `apply_neumann_boundary`
- 逻辑变化：无。

## C. 主流程文件（`fiber_levelset.m`）

### 6) 历史数据收尾函数化
- 新增局部函数：`finalize_history_data`
- 修改目的：把结尾阶段的大段数组裁剪代码收拢到单点，减少主流程噪音。
- 逻辑变化：无（切片规则保持一致）。

### 7) 路径质量与接口诊断组装函数化
- 新增局部函数：
  - `build_path_quality_history`
  - `build_interface_diagnostics`
- 修改目的：将结果打包与诊断组装从主流程中抽离，增强可读性。
- 逻辑变化：无（字段保持一致）。

### 8) 结果结构体构建函数化
- 新增局部函数：`build_results_struct`
- 修改目的：统一结果字段赋值入口，避免长段重复赋值块。
- 逻辑变化：无（结果字段保持不变）。

### 9) 局部命名规范化
- 修改项：`baseC/currC` -> `base_compliance/current_compliance`
- 修改目的：统一可读命名，避免缩写歧义。
- 逻辑变化：无。
