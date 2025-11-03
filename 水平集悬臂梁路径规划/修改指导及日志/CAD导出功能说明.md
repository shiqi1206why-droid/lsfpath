# CAD纤维路径导出功能说明

**文件**: `cad_establish.m`  
**功能**: 将优化后的水平集纤维路径导出为DXF格式CAD文件  
**创建日期**: 2025-10-23

---

## 功能特点

- ✅ 导出DXF格式CAD文件（AutoCAD R2000兼容）
- ✅ 支持多层等值线（可自定义）
- ✅ 自动图层分组（每层等值线独立图层）
- ✅ 物理坐标准确（基于dx, dy）
- ✅ 兼容主流CAD软件（AutoCAD, FreeCAD, LibreCAD等）

---

## 函数接口

```matlab
function cad_establish(lsf, dx, dy, nelx, nely, output_filename, options)
```

### 必需参数

| 参数 | 类型 | 说明 |
|------|------|------|
| `lsf` | double [nely+2, nelx+2] | 水平集函数（包含ghost cells） |
| `dx` | double | X方向网格物理间距（米） |
| `dy` | double | Y方向网格物理间距（米） |
| `nelx` | int | X方向单元数量 |
| `nely` | int | Y方向单元数量 |
| `output_filename` | string | 输出文件名（不含扩展名） |

### 可选参数（options结构体）

| 字段 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `levels` | double数组 | `[-2h, -h, 0, h, 2h]` | 要导出的等值线层级 |
| `layer_prefix` | string | `'Fiber_Path'` | DXF图层名前缀 |

其中 `h = min(dx, dy)` 为最小网格间距。

---

## 使用示例

### 示例1：从fiber_levelset结果导出

```matlab
% 运行优化
results = fiber_levelset('default');

% 导出纤维路径
cad_establish(results.lsf_final, ...
              results.params.grid.dx, ...
              results.params.grid.dy, ...
              results.params.grid.nelx, ...
              results.params.grid.nely, ...
              'my_fiber_paths');

% 生成文件: my_fiber_paths.dxf
```

### 示例2：自定义等值线层级

```matlab
% 自定义参数
opts = struct();
opts.levels = [-0.04, -0.02, 0, 0.02, 0.04];  % 自定义5层
opts.layer_prefix = 'CustomFiber';

cad_establish(lsf_final, 0.02, 0.02, 80, 50, 'custom_output', opts);
```

### 示例3：只导出零等值线（主纤维路径）

```matlab
opts = struct();
opts.levels = 0;  % 只导出零等值线

cad_establish(lsf_final, dx, dy, nelx, nely, 'main_path_only', opts);
```

---

## 输出文件格式

### DXF文件结构

```
fiber_paths.dxf
├── HEADER      % 文件头（版本、范围、单位）
├── TABLES      % 图层定义
│   ├── Fiber_Path_Zero     (零等值线)
│   ├── Fiber_Path_Minus1   (负向第1层)
│   ├── Fiber_Path_Plus1    (正向第1层)
│   └── ...
├── ENTITIES    % 几何实体（LWPOLYLINE）
│   ├── 路径1 (在Fiber_Path_Zero图层)
│   ├── 路径2 (在Fiber_Path_Zero图层)
│   └── ...
└── EOF         % 文件结束
```

### 图层命名规则

| 层级类型 | 图层名示例 |
|---------|-----------|
| 零等值线 | `Fiber_Path_Zero` |
| 正向层级 | `Fiber_Path_Plus1`, `Fiber_Path_Plus2`, ... |
| 负向层级 | `Fiber_Path_Minus1`, `Fiber_Path_Minus2`, ... |

### 坐标系统

- **原点**: (0, 0) 位于左下角
- **X轴**: 水平向右（0 → Lx）
- **Y轴**: 垂直向上（0 → Ly）
- **单位**: 米（DXF单位码=6）
- **Z坐标**: 0（2D平面）

---

## CAD软件打开方法

### AutoCAD
1. 打开AutoCAD
2. 文件 → 打开 → 选择.dxf文件
3. 在"图层"面板查看各层路径

### FreeCAD
1. 打开FreeCAD
2. 文件 → 导入 → 选择.dxf文件
3. 在树状视图中展开查看

### LibreCAD（免费开源）
1. 下载：https://librecad.org/
2. 文件 → 打开 → 选择.dxf文件
3. 右侧图层面板可切换显示

### 在线查看器
- ShareCAD: https://sharecad.org/
- A360 Viewer: https://viewer.autodesk.com/

---

## 测试验证

### 运行测试脚本

```matlab
run test_cad_export.m
```

**测试内容**：
1. 从拓扑结果导出
2. 简单圆形路径导出
3. 可视化对比

**生成文件**：
- `test_fiber_paths.dxf`
- `test_circle.dxf`

### 验证清单

- [ ] DXF文件可以在CAD软件中打开
- [ ] 路径位置与原始水平集等值线一致
- [ ] 图层正确分组（不同层级在不同图层）
- [ ] 物理尺寸正确（使用CAD测量工具验证）
- [ ] 路径连续平滑（无断裂或锯齿）

---

## 常见问题

### Q1: 为什么有些层级没有路径？

**A**: 如果水平集函数在某个层级没有等值线（例如层级超出了水平集范围），会显示警告并跳过该层级。这是正常的。

### Q2: 如何调整导出的路径数量？

**A**: 通过`options.levels`参数控制。例如：
```matlab
opts.levels = linspace(-0.06, 0.06, 7);  % 均匀分布的7层
```

### Q3: DXF文件很大怎么办？

**A**: 
- 减少等值线层级数量（只保留关键层）
- 只导出零等值线：`opts.levels = 0;`
- 使用CAD软件简化曲线

### Q4: 坐标不对齐怎么办？

**A**: 检查：
- `dx, dy`是否正确（应该是物理间距，不是单元数）
- `lsf`尺寸是否为`[nely+2, nelx+2]`
- CAD软件的单位设置是否为米

### Q5: 能否导出STL或其他格式？

**A**: 当前仅支持DXF。如需其他格式：
- DXF → STL: 使用FreeCAD转换
- DXF → SVG: 使用Inkscape转换
- DXF → STEP: 使用AutoCAD或FreeCAD转换

---

## 技术细节

### 等值线提取

使用MATLAB内置`contourc`函数：
```matlab
C = contourc(X, Y, lsf_core, [level level]);
```

返回格式：
```
C = [level1, x1, x2, ..., xn, level2, y1, y2, ...]
    [n,     y1, y2, ..., yn, m,     y1, y2, ...]
```

### DXF几何实体

使用`LWPOLYLINE`（轻量级多段线）：
- 支持度：AutoCAD R2000+
- 优点：文件小，读取快
- 特性：2D多段线，支持开放/闭合

### 颜色编码

DXF颜色索引（1-7循环）：
1. 红色
2. 黄色
3. 绿色
4. 青色
5. 蓝色
6. 品红色
7. 白色/黑色

---

## 集成到主程序

### 在fiber_levelset.m中自动调用

在`fiber_levelset.m`的结果保存部分添加：

```matlab
% 保存结果到输出结构体
results.lsf_final = lsf;
% ... 其他结果 ...

% 自动导出CAD文件
if params.debug.export_cad  % 添加配置开关
    try
        output_name = sprintf('fiber_result_%s', datestr(now, 'yyyymmdd_HHMMSS'));
        cad_establish(lsf, dx, dy, nelx, nely, output_name);
        log_message('INFO', params, 'CAD文件已导出: %s.dxf', output_name);
    catch ME
        log_message('WARN', params, 'CAD导出失败: %s', ME.message);
    end
end
```

### 配置参数添加

在`config/get_fiber_optimization_params.m`中添加：

```matlab
% CAD导出参数
base.debug.export_cad = false;  % 默认关闭，按需开启
base.cad.levels = [];           % 空则使用默认[-2h, -h, 0, h, 2h]
base.cad.layer_prefix = 'Fiber_Path';
```

---

## 性能说明

| 网格规模 | 等值线数 | 生成时间 | 文件大小 |
|---------|---------|---------|---------|
| 40×40 | 5层 | <1秒 | ~10KB |
| 80×50 | 5层 | ~1秒 | ~50KB |
| 160×100 | 5层 | ~2秒 | ~200KB |

---

## 许可与引用

本功能基于`fiber_levelset.m`项目开发，遵循相同许可协议。

如在学术工作中使用，请引用：
- 纤维路径优化论文（水平集方法用于纤维路径规划）
- 本项目GitHub仓库（如有）

---

## 更新日志

### v1.0 (2025-10-23)
- ✨ 初始版本
- ✅ DXF R2000格式导出
- ✅ 多层等值线支持
- ✅ 自动图层分组
- ✅ 物理坐标映射

