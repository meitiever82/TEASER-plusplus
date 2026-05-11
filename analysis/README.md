# TEASER++ 代码分析文档

本目录是对 TEASER++ `teaser_registration` 库的深度代码分析，依据 `codebase-analysis-skill` 规范生成。

**基准 commit**：`52a9c52` (master)
**生成方式**：手工逐行精读 + 数学推导 + 代码-公式对照
**所有断言**：均以 `[VERIFY: file:line]` 标签锚定到具体代码位置

---

## 阅读路线

### 新人路线（按依赖顺序）

1. **`00_SYSTEM_OVERVIEW.md`** — 先看这份。整套库的大局：注册管线 6 个 Stage、模块边界、参数全景。约 470 行。

2. **`01_DATA_STRUCTURES.md`** — 主求解器内部字段、参数结构、Strategy 模式下的 solver 类家族。约 880 行。

3. **`02_DATA_FLOW.md`** — `solve()` 一次调用过程中数据如何流转，每阶段输入输出维度、index 语义。约 720 行。

### 算法深度路线

4. **`ALGORITHM_01-Scalar_TLS.md`** — 1D 截断最小二乘（sweep-line adaptive voting）。scale + translation 共用。约 670 行。

5. **`ALGORITHM_02-GNC_TLS_Rotation.md`** — SO(3) 旋转的 Graduated Non-Convexity TLS。Black-Rangarajan 对偶展开、三段权重函数、加权 Procrustes。约 690 行。

6. **`ALGORITHM_03-MaxClique_TIM.md`** — TIM 构造 + Inlier 图建图 + Max-Clique（PMC + K-core）。约 720 行。

---

## 范围与不范围

### 已覆盖

- 完整管线：`RobustRegistrationSolver::solve()` 的 6 个 Stage
- 1D 估计器：`ScalarTLSEstimator::estimate` 与 `estimate_tiled`
- 旋转估计器：`GNCTLSRotationSolver`（默认）
- 图与团：`Graph`、`MaxCliqueSolver`、PMC interop
- 公开数据结构：`PointXYZ`、`PointCloud`、`Params`、`RegistrationSolution` 等
- 三大 abstract base 类家族
- 参数结构与三个嵌套枚举的语义

### 未覆盖（用户选择「核心算法」范围时排除）

- `FastGlobalRegistrationSolver` 实现细节（在算法 02 中作对比）
- `QuatroSolver` 实现细节（指出 static 变量隐患）
- `DRSCertifier` SDP 关系的证明与实现
- `teaser_features` 的 FPFH 描述子 + matcher
- `teaser_io` 的 PLY IO
- Python/MATLAB 绑定层

如需扩展上述，提示「请补 X 模块的分析」。

---

## 验证机制

每个断言（如「函数 X 在第 N 行做了 Y」）都有形如 `[VERIFY: relative/path:line]` 的标签。读者可以：

```bash
# 抽取所有 VERIFY tags
grep -ohrE "VERIFY:[^]]+\]" analysis/*.md | sort -u

# 检查某个引用
sed -n '770p' teaser/src/registration.cc   # 例：验证 GNC-TLS 函数起始
```

**统计**：共 194 个 [VERIFY:] 标签，覆盖核心代码段。

---

## 已知偏差与边界情况记录

阅读这些文档时请特别留意以下几点（散见于各文档「工程坑点」章节）：

| 主题 | 位置 | 简述 |
|------|------|------|
| `RegistrationSolution.valid` 不全 | `00_SYSTEM_OVERVIEW.md` §6.2 | 只有 max-clique ≤ 1 时被置 false；GNC 不收敛不触发 |
| 多次 `solve()` 状态泄漏 | `02_DATA_FLOW.md` §9 | `inlier_graph_`、`rotation_inliers_` 等不自动清空 |
| Quatro static 变量 | `ALGORITHM_02` §8.3 | `noise_bound_sq` 被 `static` 修饰，setParams 改不掉 |
| KCORE_HEU 不总走 K-core | `ALGORITHM_03` §5.5 | 触发依赖 `max_core > threshold·N`，否则 fallback |
| Scalar-TLS 罚分形式 | `ALGORITHM_01` §2.3.7 | 代码用 Σα（线性）而非论文的 Σα² |
| Stage 4 noise_bound 调整 | `02_DATA_FLOW.md` §5.5 | rotation 内部看到的 noise_bound 已被 *= 2/scale |

---

## 反馈与扩展

如发现 `[VERIFY:]` 标签与实际代码行号不符（可能因后续提交漂移），或某段代码理解有误，请：

1. 给出具体 `file:line`，附正确解读；
2. 修复对应文档段落（保持 `[VERIFY:]` 标签机制）；
3. 在本 README 「已知偏差」表里追加一条。

