# Step1 与 Step4 测试结果

- 生成时间: 2026-04-01
- 测试目录: `/home/again/projects/recover/debug/refactor_artifacts/checklist_remaining_single_20260401_134356`

## Step1: strict fast `accepted_steps > 0`

- 判定: **FAIL**
- 结果数据:
  - `fast.accepted_steps = 0`
  - `fast.rejected_steps = 41`
- 结论: strict fast 条件下未出现可接受步。

## Step4: TestB + strict fast accepted 约束

### 4.1 TestB (`test_exact_hj_descent_smoke`)
- 判定: **PASS**
- 原始输出:
  - `C_current = 3.777978e-06`
  - `C_hj = 3.778026e-06`
  - `delta = 4.720293e-11`
  - `dt_tiny = 1.000000e-02`
  - `max_band = 4.116696e-07`

### 4.2 Step4 最终清单判定
- 判定: **FAIL**
- 原因: Step4 在清单中要求同时满足 strict fast `accepted_steps > 0`，但当前 `accepted_steps = 0`。

---

## 摘要

- Step1: **FAIL**
- Step4: **FAIL**（其中 TestB 子测试通过，但清单级判定失败）
