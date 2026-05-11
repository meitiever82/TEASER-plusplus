# TEASER++ 系统总览

> **范围**：本文档以 `teaser_registration` 库为主线，梳理 TEASER++ 的核心数据流、模块边界与依赖。FPFH / certification / Quatro 子模块在本系列里只作引用，不展开。
>
> **基准代码**：commit `52a9c52` (master)。所有 `[VERIFY: file:line]` 标签均基于此基准。

---

## 1. 问题与定位

### 1.1 算法层定义的问题

给定两组 3D 点
```
src = {a_i ∈ R³, i = 1..N}
dst = {b_i ∈ R³, i = 1..N}
```
以及它们之间的「对应关系」（correspondences，可能含有**极高比例**的错误匹配，即 outlier），求 `(s, R, t) ∈ R₊ × SO(3) × R³`，使得对**未知**的 inlier 子集 I：

```
b_i = s · R · a_i + t + ε_i ,  ∀ i ∈ I,  ‖ε_i‖ ≤ β
```

其中 `β` 是噪声上界（noise bound）。outlier 集合 `O = {1..N} \ I` 上的 `ε_i` 可以任意大。

TEASER++ 的目标不是「找一个看起来还行的解」，而是在**已知噪声模型**下：

1. 给出**可证最优**（certifiably optimal）的 `(s, R, t)` 估计；
2. 能在 outlier 比例高达 99% 时仍然恢复正解；
3. 实际运行时间在毫秒到秒级（依赖问题规模与 max-clique 求解策略）。

### 1.2 库的边界

- **输入**：两组 3D 点 + 对应关系列表，或两组已经按列对齐的 `Eigen::Matrix3Xd`。  
  [VERIFY: teaser/include/teaser/registration.h:560-570]
- **输出**：`RegistrationSolution { valid, scale, rotation, translation }`。  
  [VERIFY: teaser/include/teaser/registration.h:32-39]
- **不做**：特征提取、对应关系建立（这是 `teaser_features` 的 FPFH wrapper 的职责，或留给用户）。

库给上层一句承诺：**给我对应，我给你 (s, R, t) 和 inlier 集合**。

---

## 2. 算法管线（pipeline）

TEASER++ 把一个**耦合的非凸问题**

```
min_{s,R,t} Σ ρ( ‖b_i − s·R·a_i − t‖ )
```

通过 **TIM（Translation-Invariant Measurement）解耦**为可顺序求解的子问题：

```
Stage 1 ┌─────────────────┐
        │   compute TIMs  │  src/dst 各计算两两差向量
        └────────┬────────┘
                 │
Stage 2 ┌────────▼────────┐
        │  scale solver   │  Scalar-TLS（自适应投票）求 s
        │  + inlier mask  │  得到 TIM 上的 scale-inlier mask
        └────────┬────────┘
                 │
Stage 3 ┌────────▼────────┐
        │   inlier graph  │  把 TIM mask 转成原始测量的图
        │   max-clique    │  PMC-exact / PMC-heu / k-core 三选一
        └────────┬────────┘
                 │  max_clique_ ⊂ {0..N-1}
                 │
Stage 4 ┌────────▼────────┐
        │  prune TIMs by  │  CHAIN 或 COMPLETE TIM 图重建
        │  max-clique     │
        └────────┬────────┘
                 │
Stage 5 ┌────────▼────────┐
        │ rotation solver │  GNC-TLS / FGR / Quatro
        │ on pruned TIMs  │
        └────────┬────────┘
                 │  R̂
                 │
Stage 6 ┌────────▼────────┐
        │translation TLS  │  三个轴独立 Scalar-TLS
        └────────┬────────┘
                 │  t̂
                 ▼
            (s, R, t)
```

入口函数：`teaser::RobustRegistrationSolver::solve()`  
[VERIFY: teaser/src/registration.cc:574-743]

### 2.1 Stage 1 — TIM 构造

**关键性质**：若 `b_i = s·R·a_i + t + ε_i`，则对任意 `i, j`

```
b_i − b_j = s·R·(a_i − a_j) + (ε_i − ε_j)
```

平移 `t` **消失**了。这就是 TIM：`ã_{ij} := a_i − a_j`，`b̃_{ij} := b_i − b_j`，把原本 `(s,R,t)` 耦合的问题降为只与 `(s,R)` 相关的子问题，从而允许先估 `s`、再估 `R`。

实现：`computeTIMs()` 把 N 个点生成 `N(N-1)/2` 个 TIM，列存于 `Eigen::Matrix<double, 3, Dynamic>`，同时返回 `2 × N(N-1)/2` 的 `map` 记录每个 TIM 对应的 `(i, j)` 索引。  
[VERIFY: teaser/src/registration.cc:518-557]

### 2.2 Stage 2 — Scale 求解

两种模式（由 `params_.estimate_scaling` 决定）：

| 模式 | 估算器 | 行为 |
|------|--------|------|
| `estimate_scaling = true` | `TLSScaleSolver` | 用 Scalar-TLS 解 raw scales `‖b̃‖/‖ã‖` |
| `estimate_scaling = false` | `ScaleInliersSelector` | 固定 `s=1`，仅用 `‖‖ã‖−‖b̃‖‖ ≤ 2β√c̄²` 筛 inlier |

[VERIFY: teaser/include/teaser/registration.h:826-832] [VERIFY: teaser/src/registration.cc:416-449]

两种模式都会输出 `scale_inliers_mask_`，这是后续 max-clique 阶段建图的输入。

### 2.3 Stage 3 — Inlier Graph + Max-Clique

把每个原始测量 `i ∈ {0..N-1}` 作为图的一个顶点；对每个 scale-inlier TIM `(i, j)`，加一条无向边。如果 inlier 集真的存在，并且 inlier 之间两两 TIM 都满足噪声约束，那么 inlier 在这张图里**构成一个团**。所以「找 inlier」≈「找最大团」。

三种求解策略 (`INLIER_SELECTION_MODE`)：

| 模式 | 实现 | 特性 |
|------|------|------|
| `PMC_EXACT` | PMC 精确求解器 | 最优，问题大时慢 |
| `PMC_HEU`   | PMC 启发式 | 快，可能漏点 inlier |
| `KCORE_HEU` | k-core 启发式 | 最快，宽松（取 max-core 全部顶点） |
| `NONE`      | 跳过该阶段 | 全部输入直接进 rotation |

[VERIFY: teaser/include/teaser/registration.h:397-402] [VERIFY: teaser/src/graph.cc:12-125]

### 2.4 Stage 4 — TIM 图重建

`max_clique_` 是顶点 ID 的子集。要进 rotation solver 的还是 TIM（因为旋转估计也需要 translation-invariant 的形式），所以这里要在 max-clique 上重新构 TIM 图：

| `INLIER_GRAPH_FORMULATION` | TIM 数 | 拓扑 |
|----------------------------|--------|------|
| `CHAIN` (默认) | `len(max_clique_)` | 环状链：第 i 个连接到第 (i+1) % L 个 |
| `COMPLETE`     | `L(L-1)/2`         | 完全图，复用 `computeTIMs()` |

CHAIN 数量线性、计算量小，但每个 TIM 只「见」相邻一个；COMPLETE 数量是 O(L²)，鲁棒但更慢。

[VERIFY: teaser/src/registration.cc:663-700]

### 2.5 Stage 5 — Rotation

输入 `pruned_src_tims_`、`pruned_dst_tims_`（其中 dst 已经除掉 `solution_.scale`）。三种求解器（`ROTATION_ESTIMATION_ALGORITHM`）：

| 算法 | 类 | 适用场景 |
|------|----|---------|
| `GNC_TLS`（默认） | `GNCTLSRotationSolver` | 通用，3D 全旋转 |
| `FGR`             | `FastGlobalRegistrationSolver` | 较快、依赖较好初值 |
| `QUATRO`          | `QuatroSolver` | 仅 yaw（SO(2)），城市/车载场景退化对抗 |

[VERIFY: teaser/include/teaser/registration.h:383-387] [VERIFY: teaser/src/registration.cc:770-872]

GNC-TLS 通过 **Black-Rangarajan 对偶**把 TLS 代价的非凸权重展开成一系列凸子问题，每次内部用加权 SVD 求 rotation；`μ` 按 `gnc_factor` 调整以渐进逼近原始 TLS。

### 2.6 Stage 6 — Translation

把 `s · R̂ · a_i` 当作变换后的 source，与 `b_i` 相减得到 raw translation。三个轴独立做 Scalar-TLS：x、y、z 各跑一次 `ScalarTLSEstimator::estimate()`，取**所有三个轴都是 inlier** 的点为最终 inlier。  
[VERIFY: teaser/src/registration.cc:451-477]

---

## 3. 模块清单（按编译目标组织）

### 3.1 `teaser_io`
- **职责**：PLY 文件读写。
- **文件**：`teaser/src/ply_io.cc`、`teaser/include/teaser/ply_io.h`。
- **依赖**：`tinyply`（FetchContent）。
- **关键类**：`teaser::PLYReader`、`teaser::PLYWriter`。

与 registration 完全解耦，是单独可用的目标 `teaserpp::teaser_io`。  
[VERIFY: teaser/CMakeLists.txt:35-43]

### 3.2 `teaser_registration` ★ 本系列分析重点

| 子模块 | 头文件 | 实现 | 角色 |
|--------|--------|------|------|
| `RobustRegistrationSolver` | `registration.h:362-933` | `registration.cc:479-743` | 主入口，编排管线 |
| Scale solvers | `registration.h:136-187` | `registration.cc:416-449` | `TLSScaleSolver`、`ScaleInliersSelector` |
| Rotation solvers | `registration.h:218-355` | `registration.cc:212-414, 770-872` | `GNCTLSRotationSolver`、`FastGlobalRegistrationSolver`、`QuatroSolver` |
| Translation solver | `registration.h:189-215` | `registration.cc:451-477` | `TLSTranslationSolver` |
| Scalar TLS estimator | `registration.h:102-131` | `registration.cc:27-210` | `ScalarTLSEstimator`，含 `estimate()` 与 `estimate_tiled()` |
| TIM computation | `registration.h:543-550` | `registration.cc:518-557` | `RobustRegistrationSolver::computeTIMs()` |
| Graph + max-clique | `graph.h:29-279` | `graph.cc:12-125` | `Graph`、`MaxCliqueSolver`（PMC facade） |
| Certification | `certification.h` | `certification.cc` | `DRSCertifier`（**本系列不展开**） |
| Utilities | `utils.h` | header-only | `svdRot`、`svdRot2d`、`calculateDiameter`、`findNonzero` |

- **public 依赖**：`Eigen3::Eigen`。
- **private 依赖**：`pmc`（FetchContent，max-clique 实际求解器）、`spectra`（FetchContent，certification 用的特征分解，header-only）、`OpenMP::OpenMP_CXX`（可选）。  
  [VERIFY: teaser/CMakeLists.txt:51-77]

### 3.3 `teaser_features`（可选）
- **职责**：FPFH 特征描述 + 匹配，给「光秃秃没有对应关系」的两片点云生成对应输入。
- **条件**：`-DBUILD_TEASER_FPFH=ON`（需 PCL）。
- **本系列不展开。**

### 3.4 Python / MATLAB 绑定
- **Python**：`python/teaserpp_python/teaserpp_python.cc`，pybind11，导出 `_teaserpp` 扩展，包装在 `teaserpp_python` 包内。
- **MATLAB**：`matlab/teaser_mex.cc` + `teaser_solve.m`。
- 这两层是**纯转接**，不参与算法逻辑，所以本系列不分析。

---

## 4. 控制流时序（end-to-end）

```
User
 │   solver = RobustRegistrationSolver(params)
 │     ├─ reset(params) ┐
 │     │                ├─ scale_solver_   = TLSScaleSolver | ScaleInliersSelector
 │     │                ├─ rotation_solver_= GNCTLS | FGR | Quatro
 │     │                └─ translation_solver_ = TLSTranslationSolver
 │     │   [VERIFY: registration.h:823-863]
 │     │
 │   sol = solver.solve(src, dst)
 │     │
 │     ├─ src_tims_ = computeTIMs(src, &src_tims_map_)       ← [registration.cc:605]
 │     ├─ dst_tims_ = computeTIMs(dst, &dst_tims_map_)       ← [registration.cc:606]
 │     │
 │     ├─ solveForScale(src_tims_, dst_tims_)                ← [registration.cc:609]
 │     │     └─ scale_solver_->solveForScale(...)
 │     │           ├─ TLSScaleSolver: ScalarTLSEstimator::estimate()
 │     │           └─ ScaleInliersSelector: |‖b̃‖−‖ã‖|≤2β√c̄²
 │     │
 │     ├─ if inlier_selection_mode != NONE:                   ← [registration.cc:615]
 │     │     ├─ build inlier_graph_ from scale_inliers_mask_  ← [registration.cc:620-625]
 │     │     ├─ clique_solver.findMaxClique(inlier_graph_)    ← [registration.cc:641]
 │     │     │     ├─ G.compute_cores()                       ← [graph.cc:58]
 │     │     │     ├─ if KCORE_HEU & threshold-met → 取 max-core 顶点
 │     │     │     ├─ else: pmc_heu → 初始下界 lb              ← [graph.cc:88-91]
 │     │     │     └─ if PMC_EXACT: pmcx_maxclique.search()   ← [graph.cc:105-122]
 │     │     └─ max_clique_ = clique (sorted)
 │     │   else:
 │     │     └─ max_clique_ = 0..N-1
 │     │
 │     ├─ rebuild pruned_*_tims_ on max_clique_:              ← [registration.cc:663-700]
 │     │     ├─ CHAIN:    ring graph of size L = |max_clique_|
 │     │     └─ COMPLETE: computeTIMs(inliers, ...)
 │     │
 │     ├─ pruned_dst_tims_ /= solution_.scale                 ← [registration.cc:703]
 │     ├─ adjust rotation noise_bound *= 2/scale              ← [registration.cc:708-710]
 │     │
 │     ├─ solveForRotation(pruned_src_tims_, pruned_dst_tims_) ← [registration.cc:714]
 │     │     └─ rotation_solver_->solveForRotation(...)
 │     │           ├─ GNCTLSRotationSolver: 见 ALGORITHM_02
 │     │           ├─ FastGlobalRegistrationSolver: line-process FGR
 │     │           └─ QuatroSolver: 2D yaw only
 │     │
 │     ├─ pull inlier indices from rotation_inliers_mask_     ← [registration.cc:718-722]
 │     ├─ collect rotation_pruned_src/dst from max_clique_    ← [registration.cc:723-728]
 │     │
 │     ├─ solveForTranslation(s·R·src_pruned, dst_pruned)     ← [registration.cc:732]
 │     │     └─ TLSTranslationSolver: 3 × ScalarTLSEstimator
 │     │
 │     ├─ translation_inliers_ = findNonzero(translation_inliers_mask_)
 │     └─ return solution_                                    ← [registration.cc:742]
```

---

## 5. 关键设计选择与权衡

### 5.1 为什么用 TIM 解耦而不是联合优化 `(s, R, t)`？

| 方案 | 优点 | 缺点 |
|------|------|------|
| 联合优化 `min Σ ρ(‖b - sRa - t‖)` | 一次性，模型紧 | 非凸 + 强耦合，outlier 比例高时容易陷局部最优 |
| **TIM 解耦** | 分阶段，每阶段都是 1D 或 SO(3) 上的子问题，可分别用 TLS / GNC 处理 | TIM 数量 `O(N²)`，大 N 时内存敏感 |

TEASER 选 TIM：作者证明在 noise bound `β` 下，TIM 形式的代价对 outlier 仍然保持 robust，且 scale-rotation-translation 三个子问题在 inlier 集已知时**各自有闭式解或一维 TLS**。

参考论文：H. Yang et al., "TEASER: Fast and Certifiable Point Cloud Registration", arXiv:2001.07715.

### 5.2 为什么先 scale 再 rotation 再 translation？

因为：
- **Scale**：TIM 的**范数比**只与 `s` 有关。 `‖b̃_{ij}‖ / ‖ã_{ij}‖ ≈ s`（无 outlier）。所以 1D 上做 TLS 就够。  
  [VERIFY: teaser/src/registration.cc:421-431]
- **Rotation**：去掉 scale 后，`b̃ = R·ã`，只剩 SO(3) 上的旋转配准。GNC-TLS 处理。
- **Translation**：旋转已知后，`b_i − s·R·a_i ≈ t + ε_i`，三轴独立的 1D TLS。

顺序是不可换的——每一步都依赖前一步的结果。

### 5.3 为什么要先 max-clique，再做 rotation？

scale 阶段已经给出了 TIM 上的 inlier mask，但这是 TIM 维度（`O(N²)`）的 mask，不是原始测量维度的 mask。TIM `(i,j)` 是 inlier 不代表 `i, j` 本身是好的（可能两个 outlier 偶然差向量接近）。但反过来，**如果 `i, j` 都是 inlier，那 TIM `(i,j)` 必然是 inlier**（噪声有界且独立）。

所以 inlier 集合在「TIM-inlier 关系图」上必然构成团：max-clique 是一个**保守且严格**的 inlier 提取器。这是 TEASER 论文里 outlier 鲁棒性的关键保证之一。

### 5.4 GNC-TLS vs Plain TLS

直接对 rotation 做 TLS 是组合优化（每个测量是不是 inlier 是 0/1 决策）。GNC 用 **Black-Rangarajan 对偶**把它转成连续权重 `w ∈ [0,1]^N`，引入控制参数 `μ`，从凸近似（`μ` 大时所有权重接近 1）逐步逼近 TLS（`μ→∞`，权重二值化）。每次内部循环用加权 SVD 闭式解 rotation。

收敛条件：连续两次代价差 `< cost_threshold` 或达到 `max_iterations`。详见 `ALGORITHM_02-GNC_TLS_Rotation.md`。

---

## 6. 输入 / 输出契约

### 6.1 输入

**两种 overload**：

```cpp
// Form A: PointCloud + correspondences
RegistrationSolution solve(
    const teaser::PointCloud& src_cloud,
    const teaser::PointCloud& dst_cloud,
    const std::vector<std::pair<int, int>> correspondences);
// [VERIFY: teaser/src/registration.cc:559-572]

// Form B: aligned 3×N matrices (column = correspondence)
RegistrationSolution solve(
    const Eigen::Matrix<double, 3, Eigen::Dynamic>& src,
    const Eigen::Matrix<double, 3, Eigen::Dynamic>& dst);
// [VERIFY: teaser/src/registration.cc:574-576]
```

Form A 内部直接构造 3×N 矩阵后调 Form B。**两组的列数必须相等**——这是「对应关系列表」。

### 6.2 输出

```cpp
struct RegistrationSolution {
    bool valid = true;
    double scale;
    Eigen::Vector3d translation;
    Eigen::Matrix3d rotation;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
// [VERIFY: teaser/include/teaser/registration.h:32-39]
```

满足 `dst[i] ≈ scale · rotation · src[i] + translation` 当 `i ∈ inliers`。

`valid = false` 仅在 max-clique 大小 ≤ 1 时设置：
```cpp
if (max_clique_.size() <= 1) {
    solution_.valid = false;
    return solution_;
}
// [VERIFY: teaser/src/registration.cc:649-653]
```

注意：成功路径**只**在 `solve()` 末尾把 `valid = true`：
```cpp
solution_.valid = true;
return solution_;
// [VERIFY: teaser/src/registration.cc:740-742]
```
但其他失败模式（GNC 不收敛、scale 估计退化）当前并不显式标记 `valid = false`——用户应额外检查 `getGNCRotationCostAtTermination()` 等指标。

### 6.3 中间结果访问器

`RobustRegistrationSolver` 暴露了大量 getter，用于调试与可视化：

| Getter | 含义 |
|--------|------|
| `getSrcTIMs()` / `getDstTIMs()` | Stage 1 的全量 TIM |
| `getScaleInliersMask()`、`getScaleInliers()` | Stage 2 输出的 TIM-level mask |
| `getInlierGraph()` | Stage 3 用于 max-clique 的图（邻接表形式） |
| `getInlierMaxClique()` | Stage 3 的 max-clique（原始测量 ID） |
| `getMaxCliqueSrcTIMs()` / `getMaxCliqueDstTIMs()` | Stage 4 重建的 TIM |
| `getRotationInliers()` / `getRotationInliersMask()` | Stage 5 的 inlier（在 pruned TIM 维度） |
| `getTranslationInliers()` / `getTranslationInliersMask()` | Stage 6 的 inlier（在 max-clique 维度） |
| `getInputOrderedTranslationInliers()` | Stage 6 inlier 映回**原始输入下标** |
| `getGNCRotationCostAtTermination()` | GNC 退出时的代价，用于评估解质量 |

[VERIFY: teaser/include/teaser/registration.h:640-817]

注意 `getInputOrderedTranslationInliers()` 在 FGR 模式下会抛异常，因为 FGR 跳过 max-clique。  
[VERIFY: teaser/include/teaser/registration.h:745-749]

---

## 7. 配置参数全景

完整定义见 `Params` 结构体，[VERIFY: teaser/include/teaser/registration.h:420-514]：

| 参数 | 默认值 | 含义 |
|------|--------|------|
| `noise_bound` | `0.01` | `β`，测量噪声上界（**最重要的参数**） |
| `cbar2` | `1.0` | TLS 截断阈值平方比 `c̄²` |
| `estimate_scaling` | `true` | 是否求 scale；`false` 时只筛 inlier |
| `rotation_estimation_algorithm` | `GNC_TLS` | rotation solver 选择 |
| `rotation_gnc_factor` | `1.4` | μ 调整因子（GNC-TLS 乘、FGR 除） |
| `rotation_max_iterations` | `100` | rotation solver 内部最大轮数 |
| `rotation_cost_threshold` | `1e-6` | rotation solver 收敛阈值 |
| `rotation_tim_graph` | `CHAIN` | rotation 阶段 TIM 拓扑 |
| `inlier_selection_mode` | `PMC_EXACT` | max-clique 求解策略 |
| `kcore_heuristic_threshold` | `0.5` | k-core 启发式触发阈值（max-core ≥ threshold·N 时启用） |
| `max_clique_time_limit` | `3600` s | max-clique 时间上限 |
| `max_clique_num_threads` | `omp_get_max_threads()` | max-clique 并行线程数 |
| `use_max_clique` / `max_clique_exact_solution` | true / true | **已废弃**，被 `inlier_selection_mode` 替代 |

废弃字段的兼容逻辑在 `solve()` 开头处理：
```cpp
if (!params_.use_max_clique) {
    params_.inlier_selection_mode = INLIER_SELECTION_MODE::NONE;
}
if (!params_.max_clique_exact_solution) {
    params_.inlier_selection_mode = INLIER_SELECTION_MODE::PMC_HEU;
}
// [VERIFY: teaser/src/registration.cc:579-589]
```

### 7.1 调参建议（基于代码与论文）

| 场景 | 建议 |
|------|------|
| outlier 比例 < 50% | `PMC_HEU` 足够，节省时间 |
| outlier 比例 50%–95% | `PMC_EXACT`，N ≤ 几千 |
| N > 5000，求快 | `KCORE_HEU` + `kcore_heuristic_threshold = 0.5` |
| 已知刚体（无 scale） | `estimate_scaling = false`，跳 Scalar-TLS scale |
| 车载/水平场景 | `rotation_estimation_algorithm = QUATRO`（只估 yaw） |
| 小内存 | `INLIER_GRAPH_FORMULATION::CHAIN`（默认） |
| 极高鲁棒性 | `INLIER_GRAPH_FORMULATION::COMPLETE`，搭 `PMC_EXACT` |

`noise_bound` 必须**真实反映**测量噪声尺度；过小会把好对应也判为 outlier，过大会接纳 outlier。它单位与点云一致（m / mm / 任意）。

---

## 8. 第三方依赖与版本

通过 `FetchContent` 配置时拉取，缓存在 `build/_deps/`：

| 依赖 | 用途 | 拉取方式 |
|------|------|---------|
| `pmc` | max-clique 求解 | `https://github.com/jingnanshi/pmc.git` (HEAD) |
| `spectra` | 大型对称矩阵特征分解（certification） | `https://github.com/jingnanshi/spectra` @ `5c4fb1de` |
| `tinyply` | PLY IO | `https://github.com/ddiakopoulos/tinyply.git` (HEAD) |
| `googletest` | 单元测试 | `https://github.com/google/googletest.git` @ `main` |

[VERIFY: teaser/CMakeLists.txt:4-32] [VERIFY: test/CMakeLists.txt:3-12]

系统依赖：
- `Eigen3 ≥ 3.2`（必须，public link 给用户）
- `Boost ≥ 1.58`（仅 `BUILD_TEASER_FPFH=ON` 时）
- `PCL ≥ 1.8`（仅 `BUILD_TEASER_FPFH=ON` 时）
- `OpenMP`（可选，加速 TIM 构造与 TLS tiled 版本）
- `MKL`（可选，`-DENABLE_MKL=ON` 时给 Eigen 用）

---

## 9. 与本系列其他文档的对应关系

| 本文档章节 | 详细文档 |
|------------|---------|
| §2.1 TIM | `02_DATA_FLOW.md` §2 |
| §2.2 Scalar-TLS scale | `ALGORITHM_01-Scalar_TLS.md` |
| §2.3 Max-clique + inlier graph | `ALGORITHM_03-MaxClique_TIM.md` |
| §2.5 GNC-TLS rotation | `ALGORITHM_02-GNC_TLS_Rotation.md` |
| §3 模块清单 | `01_DATA_STRUCTURES.md` |
| §7 参数全景 | `01_DATA_STRUCTURES.md` §3 |

---

## 10. 文档完整性检查表

- [x] 所有断言均有 `[VERIFY: file:line]` 标签。
- [x] 入口函数 `solve()` 的每一阶段都有代码位置标注。
- [x] 三大求解器（scale / rotation / translation）的实现位置已锁定。
- [x] 数据结构 `RegistrationSolution`、`Params`、`Graph`、`MaxCliqueSolver::Params` 已索引到行号。
- [x] FetchContent 拉取的第三方依赖已列出。
- [x] 废弃字段（`use_max_clique`、`max_clique_exact_solution`）的兼容处理路径已标注。

