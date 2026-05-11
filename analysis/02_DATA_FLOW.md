# TEASER++ 数据流分析

> **范围**：跟踪一次 `RobustRegistrationSolver::solve()` 调用过程中，数据如何在各模块间流转——从用户输入的 N 对点云对应一路到最终 `(s, R, t)` 与 inlier 索引。
>
> 每一阶段都同时给出：输入数据形态 → 算子 → 输出数据形态 → 下游消费者。

---

## 1. 顶层数据流

```
                       ┌─────────────────────────────────────┐
   User input          │ src, dst : 3 × N (correspondences)  │
                       └──────────────┬──────────────────────┘
                                      │
                                      ▼
   Stage 1: TIM         ┌──────────────────────────────────────────────┐
   computeTIMs()        │ src_tims_, dst_tims_ : 3 × M  (M=N(N-1)/2)    │
                        │ src_tims_map_, dst_tims_map_ : 2 × M           │
                        └──────────────┬───────────────────────────────┘
                                       │
                                       ▼
   Stage 2: scale       ┌──────────────────────────────────────────────┐
   ScalarTLS / select   │ solution_.scale : double                      │
                        │ scale_inliers_mask_ : 1 × M (bool)             │
                        └──────────────┬───────────────────────────────┘
                                       │
                                       ▼
   Stage 3: inlier      ┌──────────────────────────────────────────────┐
   graph + max-clique   │ inlier_graph_ : Graph(N vertices)             │
                        │ max_clique_ : vector<int>, size L             │
                        └──────────────┬───────────────────────────────┘
                                       │
                                       ▼
   Stage 4: TIM rebuild ┌──────────────────────────────────────────────┐
   on max-clique        │ pruned_src_tims_, pruned_dst_tims_ : 3 × Lp   │
                        │   Lp = L         (CHAIN)                       │
                        │   Lp = L(L-1)/2  (COMPLETE)                    │
                        │ {src,dst}_tims_map_rotation_ : 2 × Lp          │
                        │ ★ pruned_dst_tims_ /= scale                    │
                        └──────────────┬───────────────────────────────┘
                                       │
                                       ▼
   Stage 5: rotation    ┌──────────────────────────────────────────────┐
   GNC-TLS / FGR /      │ solution_.rotation : 3 × 3                    │
   Quatro               │ rotation_inliers_mask_ : 1 × Lp                │
                        │ rotation_inliers_ : vector<int>                │
                        └──────────────┬───────────────────────────────┘
                                       │
                                       ▼
   Stage 6: translation ┌──────────────────────────────────────────────┐
   3 × ScalarTLS        │ solution_.translation : 3 × 1                 │
                        │ translation_inliers_mask_ : 1 × L              │
                        │ translation_inliers_ : vector<int>             │
                        └──────────────┬───────────────────────────────┘
                                       │
                                       ▼
                       ┌─────────────────────────────────────┐
   Output              │ RegistrationSolution { valid, s,    │
                       │                        R, t }       │
                       └─────────────────────────────────────┘
```

---

## 2. Stage 1 — TIM Construction

### 2.1 数据契约

**输入**：
```
v : Eigen::Matrix<double, 3, N>
```
每列一个 3D 点。在主入口 `solve()` 中，`src` 与 `dst` 分别走一次。  
[VERIFY: teaser/src/registration.cc:605-606]

**输出**：
```
vtilde : Eigen::Matrix<double, 3, N*(N-1)/2>
map    : Eigen::Matrix<int,    2, N*(N-1)/2>
```

### 2.2 实现细节

`computeTIMs()` 的完整实现 [VERIFY: teaser/src/registration.cc:518-557]：

```cpp
Eigen::Matrix<double, 3, Eigen::Dynamic>
RobustRegistrationSolver::computeTIMs(
    const Eigen::Matrix<double, 3, Eigen::Dynamic>& v,
    Eigen::Matrix<int, 2, Eigen::Dynamic>* map) {

    auto N = v.cols();
    Eigen::Matrix<double, 3, Eigen::Dynamic> vtilde(3, N * (N - 1) / 2);
    map->resize(2, N * (N - 1) / 2);

#pragma omp parallel for shared(N, v, vtilde, map)
    for (size_t i = 0; i < N - 1; i++) {
        // i=0: 加 N-1 个 TIM ; i=k: 加 N-1-k 个 TIM
        // 累加得起始下标公式：
        size_t segment_start_idx = i * N - i * (i + 1) / 2;
        size_t segment_cols      = N - 1 - i;

        Eigen::Matrix<double, 3, 1> m = v.col(i);
        Eigen::Matrix<double, 3, Eigen::Dynamic> temp
            = v - m * Eigen::MatrixXd::Ones(1, N);

        vtilde.middleCols(segment_start_idx, segment_cols)
            = temp.rightCols(segment_cols);

        Eigen::Matrix<int, 2, Eigen::Dynamic> map_addition(2, N);
        for (size_t j = 0; j < N; ++j) {
            map_addition(0, j) = i;
            map_addition(1, j) = j;
        }
        map->middleCols(segment_start_idx, segment_cols)
            = map_addition.rightCols(segment_cols);
    }
    return vtilde;
}
```

### 2.3 索引推导

每个 TIM 列 `c` 对应一对 `(i, j)`，其中 `i < j`，意义是 `vtilde.col(c) = v.col(j) - v.col(i)`。

按 `i = 0, 1, 2, ..., N-2` 分段：
- `i=0` 的段含 `j = 1, 2, ..., N-1`，共 `N-1` 列；
- `i=1` 的段含 `j = 2, 3, ..., N-1`，共 `N-2` 列；
- 一般地，`i=k` 的段含 `N-1-k` 列。

**起始下标公式**：
```
segment_start_idx(i) = (N-1) + (N-2) + ... + (N-i)
                     = i * N - i * (i+1) / 2
```

验算 `i=3, N=10`：
- `start = 3 × 10 - 3 × 4 / 2 = 30 - 6 = 24`
- 前 3 段大小 `9 + 8 + 7 = 24` ✓

代码 [VERIFY: teaser/src/registration.cc:537-538] 使用的就是这个公式。

### 2.4 `map` 字段语义

`map(0, c) = i`，`map(1, c) = j`，且总有 `i < j`。

这个约定决定了下游 Stage 3 建图代码可以无歧义地取两端点：

```cpp
inlier_graph_.addEdge(src_tims_map_(0, i), src_tims_map_(1, i));
// [VERIFY: teaser/src/registration.cc:623]
```

### 2.5 OpenMP 并行性

外层 `i` 循环被 `#pragma omp parallel for` 并行化。各 `i` 写入 `vtilde` 的**不相交列段**，所以无数据竞争。`map` 同理。

**性能影响**：N=1000 时 M ≈ 500k，矩阵分配和填充约占 TIM 阶段的 80%。多核加速通常 4–8 倍（取决于内存带宽）。

### 2.6 关键不变式

- `vtilde.cols() == N * (N - 1) / 2`，**不是** `N * N` 也**不是** `N * (N-1)`。代码用「上三角」存法。
- `vtilde.col(c) = v.col(map(1,c)) - v.col(map(0,c))`，约定方向：「第二个点」减去「第一个点」。
- `src_tims_map_` 与 `dst_tims_map_` **数值上相同**（因为 N 相同、构造逻辑相同），但作为成员各自存储。下游只用 `src_tims_map_`，但代码两者都构造。[VERIFY: teaser/src/registration.cc:605-606]

### 2.7 数据传递

输出 `src_tims_`、`dst_tims_` 直接进入 Stage 2 的 `solveForScale(src_tims_, dst_tims_)`。  
[VERIFY: teaser/src/registration.cc:609]

`src_tims_map_` 则在 Stage 3 的图构造里被消费。

---

## 3. Stage 2 — Scale Estimation

### 3.1 两条分支

`reset()` 时根据 `params_.estimate_scaling` 决定 scale solver 类型：

```cpp
if (params_.estimate_scaling) {
    setScaleEstimator(std::make_unique<TLSScaleSolver>(noise_bound, cbar2));
} else {
    setScaleEstimator(std::make_unique<ScaleInliersSelector>(noise_bound, cbar2));
}
// [VERIFY: teaser/include/teaser/registration.h:826-832]
```

调用入口统一为：

```cpp
double solveForScale(...) {
    scale_inliers_mask_.resize(1, v1.cols());
    scale_solver_->solveForScale(v1, v2, &(solution_.scale), &scale_inliers_mask_);
    return solution_.scale;
}
// [VERIFY: teaser/src/registration.cc:745-751]
```

`scale_inliers_mask_` 的列数 = TIM 数 M。

### 3.2 `TLSScaleSolver` 路径

```cpp
v1_dist = sqrt(sum(src.^2, dim=0))   // 1×M：每个 src TIM 的范数
v2_dist = sqrt(sum(dst.^2, dim=0))   // 1×M：每个 dst TIM 的范数
raw_scales = v2_dist ./ v1_dist      // 1×M：候选 scale
β = 2 * noise_bound_ * sqrt(cbar2_)  // 标量
alphas = β ./ v1_dist                // 1×M：每个候选的容许半径
tls_estimator_.estimate(raw_scales, alphas, &scale, inliers)
// [VERIFY: teaser/src/registration.cc:419-431]
```

**数学含义**：对 inlier TIM，`‖b̃‖ ≈ s · ‖ã‖`，所以 `‖b̃‖/‖ã‖ ≈ s`。噪声扰动下，`raw_scales` 是 `s` 附近的带噪测量；`alphas` 是每个测量的不确定半径，由 `β = 2β_noise · √c̄²` 和 `‖ã‖` 决定（大 ã 容许更小相对误差）。

Scalar-TLS 见 `ALGORITHM_01-Scalar_TLS.md`。

### 3.3 `ScaleInliersSelector` 路径

```cpp
*scale = 1;
v1_dist = ‖src.col(i)‖
v2_dist = ‖dst.col(i)‖
β = 2 * noise_bound_ * sqrt(cbar2_)
*inliers = |v1_dist - v2_dist| ≤ β  // element-wise
// [VERIFY: teaser/src/registration.cc:436-449]
```

这里不真的求 scale，只用「两端范数差是否在噪声范围内」判 inlier。适用于已知 `s = 1` 的刚体配准（绝大多数场景）。

### 3.4 数据传递

输出：
- `solution_.scale`（标量）→ Stage 4 用于 `pruned_dst_tims_ /= solution_.scale` 与 `rotation_solver_` noise bound 调整；
- `scale_inliers_mask_`（1×M）→ Stage 3 用于建图。

---

## 4. Stage 3 — Inlier Graph + Max-Clique

### 4.1 数据契约

**输入**：
- `scale_inliers_mask_` (1×M bool)
- `src_tims_map_` (2×M int)
- 原始测量数 `N = src.cols()`

**输出**：
- `inlier_graph_`（顶点数 N，边数 ≤ M）
- `max_clique_`（`vector<int>`，size L ≤ N）

### 4.2 图构造

```cpp
if (params_.inlier_selection_mode != INLIER_SELECTION_MODE::NONE) {
    inlier_graph_.populateVertices(src.cols());
    for (size_t i = 0; i < scale_inliers_mask_.cols(); ++i) {
        if (scale_inliers_mask_(0, i)) {
            inlier_graph_.addEdge(src_tims_map_(0, i), src_tims_map_(1, i));
        }
    }
    ...
}
// [VERIFY: teaser/src/registration.cc:615-625]
```

**算子语义**：每个 scale-inlier TIM `(i, j)` 在测量图上加一条无向边 `(i, j)`。

- 顶点是原始**测量索引**（`0` 到 `N-1`）；
- 边只在「该测量对的 TIM 通过 scale 检查」时存在。

### 4.3 max-clique 求解

```cpp
teaser::MaxCliqueSolver::Params clique_params;
// 根据 inlier_selection_mode 设置 solver_mode
// ...
clique_params.time_limit              = params_.max_clique_time_limit;
clique_params.kcore_heuristic_threshold = params_.kcore_heuristic_threshold;
clique_params.num_threads             = params_.max_clique_num_threads;

teaser::MaxCliqueSolver clique_solver(clique_params);
max_clique_ = clique_solver.findMaxClique(inlier_graph_);
std::sort(max_clique_.begin(), max_clique_.end());
// [VERIFY: teaser/src/registration.cc:627-642]
```

`findMaxClique()` 内部把 `inlier_graph_` 转 PMC CSR 表示，调 PMC 库求解。详见 `ALGORITHM_03-MaxClique_TIM.md`。

### 4.4 NONE 分支：跳过 max-clique

```cpp
} else {
    max_clique_.reserve(src.cols());
    for (size_t i = 0; i < src.cols(); ++i) {
        max_clique_.push_back(i);
    }
}
// [VERIFY: teaser/src/registration.cc:654-660]
```

把 `max_clique_` 填成 `0..N-1`，相当于「全部输入都是 inlier」，让 Stage 4/5 在 GNC-TLS 内部再做 outlier 剔除。这适用于已知输入质量很高的情况。

### 4.5 退化情况

```cpp
if (max_clique_.size() <= 1) {
    TEASER_DEBUG_INFO_MSG("Clique size too small. Abort.");
    solution_.valid = false;
    return solution_;
}
// [VERIFY: teaser/src/registration.cc:649-653]
```

max-clique 大小 ≤ 1 时无法估计 rotation（需要至少 3 个非共线点），直接返回 `valid = false`。

**注意**：max-clique 大小 = 2 时**不会**触发这个判断，但 2 点也只够 1 个 TIM（CHAIN 模式下），rotation 也是欠定的——库不会拒绝，但解必然不可靠。实务上 inlier 至少应有 6–10 个才稳定。

### 4.6 数据传递

`max_clique_` → Stage 4 的 TIM 重建与 Stage 6 的 inlier 索引映射。

---

## 5. Stage 4 — Pruned TIM Reconstruction

### 5.1 两种拓扑

由 `params_.rotation_tim_graph` 决定：

```cpp
if (params_.rotation_tim_graph == INLIER_GRAPH_FORMULATION::CHAIN) {
    // 环状链
    ...
} else {
    // 完全图
    ...
}
// [VERIFY: teaser/src/registration.cc:663-700]
```

### 5.2 CHAIN 实现

```cpp
pruned_src_tims_.resize(3, max_clique_.size());
pruned_dst_tims_.resize(3, max_clique_.size());
src_tims_map_rotation_.resize(2, max_clique_.size());
dst_tims_map_rotation_.resize(2, max_clique_.size());

for (size_t i = 0; i < max_clique_.size(); ++i) {
    const auto& root = max_clique_[i];
    int leaf;
    if (i != max_clique_.size() - 1) {
        leaf = max_clique_[i + 1];
    } else {
        leaf = max_clique_[0];  // ★ 环回到第 0 个
    }
    pruned_src_tims_.col(i) = src.col(leaf) - src.col(root);
    pruned_dst_tims_.col(i) = dst.col(leaf) - dst.col(root);

    dst_tims_map_rotation_(0, i) = leaf;
    dst_tims_map_rotation_(1, i) = root;
    src_tims_map_rotation_(0, i) = leaf;
    src_tims_map_rotation_(1, i) = root;
}
// [VERIFY: teaser/src/registration.cc:666-686]
```

**结构**：把 `max_clique_` 看作环，第 i 条 TIM 是相邻两个测量的差向量。最后一个测量与第 0 个连接，形成闭环。

| i | root | leaf | TIM |
|---|------|------|-----|
| 0 | max_clique_[0] | max_clique_[1] | leaf − root |
| 1 | max_clique_[1] | max_clique_[2] | leaf − root |
| ... | ... | ... | ... |
| L−1 | max_clique_[L−1] | max_clique_[0] | leaf − root |

总共 `L` 条 TIM，每条「跨度」是 1。

**为什么是环而不是开链？** 环保证每个顶点都参与两条 TIM，结构对称，无端点效应。

### 5.3 COMPLETE 实现

```cpp
Eigen::Matrix<double, 3, Eigen::Dynamic> src_inliers(3, max_clique_.size());
Eigen::Matrix<double, 3, Eigen::Dynamic> dst_inliers(3, max_clique_.size());
for (size_t i = 0; i < max_clique_.size(); ++i) {
    src_inliers.col(i) = src.col(max_clique_[i]);
    dst_inliers.col(i) = dst.col(max_clique_[i]);
}
pruned_dst_tims_ = computeTIMs(dst_inliers, &dst_tims_map_rotation_);
pruned_src_tims_ = computeTIMs(src_inliers, &src_tims_map_rotation_);
// [VERIFY: teaser/src/registration.cc:687-700]
```

复用 Stage 1 的 `computeTIMs()`，生成 `L(L-1)/2` 个 TIM。

### 5.4 Scale 去除

无论哪种拓扑，紧接着：

```cpp
pruned_dst_tims_ *= (1 / solution_.scale);
// [VERIFY: teaser/src/registration.cc:703]
```

把 dst TIM 缩放回与 src 同尺度，让 rotation solver 专注 SO(3)：

```
b̃ = s · R · ã + noise
b̃ / s = R · ã + noise / s
```

### 5.5 Rotation solver noise bound 调整

```cpp
auto params = rotation_solver_->getParams();
params.noise_bound *= (2 / solution_.scale);
rotation_solver_->setParams(params);
// [VERIFY: teaser/src/registration.cc:708-710]
```

**因子 `2`** 来源：原始测量噪声 `‖ε_i‖ ≤ β`，TIM 噪声为两点噪声之差，`‖ε_i − ε_j‖ ≤ 2β`（三角不等式）。

**除以 scale**：因为 dst 已经除掉了 scale，TIM 的实际噪声也变成 `2β / s`。

### 5.6 数据传递

- `pruned_src_tims_`、`pruned_dst_tims_` → Stage 5 rotation solver；
- `rotation_solver_` 的 `params_.noise_bound` 已就地更新。

---

## 6. Stage 5 — Rotation Estimation

### 6.1 入口

```cpp
Eigen::Matrix3d solveForRotation(const Eigen::Matrix<double, 3, Eigen::Dynamic>& v1,
                                 const Eigen::Matrix<double, 3, Eigen::Dynamic>& v2) {
    rotation_inliers_mask_.resize(1, v1.cols());
    rotation_solver_->solveForRotation(v1, v2, &(solution_.rotation), &rotation_inliers_mask_);
    return solution_.rotation;
}
// [VERIFY: teaser/src/registration.cc:762-768]
```

`rotation_inliers_mask_` 维度 = `pruned_src_tims_.cols()` = `Lp`（CHAIN 时为 L，COMPLETE 时为 L(L−1)/2）。

### 6.2 三种实现的数据形态

均以 `(src_3xLp, dst_3xLp)` 为输入，输出 `R̂ ∈ SO(3)` + mask。

#### GNCTLSRotationSolver
- 内部状态：`weights : 1 × Lp`、`residuals_sq : 1 × Lp`、`mu : double`、`cost_ : double`
- 每轮：加权 SVD 解 R → 计算残差 → 闭式更新 weights → μ ← μ · gnc_factor
- 输出：`*rotation`、`*inliers = (weights >= 0.5)`
- [VERIFY: teaser/src/registration.cc:806-871]

#### FastGlobalRegistrationSolver
- 内部状态：`l_pq : 1 × Lp`、`mu : double`、`cost_ : double`
- 每轮：line-process 权重 → SVD 解 R → μ ← μ / gnc_factor
- 输出：`*rotation`、`*inliers = l_pq.cast<bool>()`
- [VERIFY: teaser/src/registration.cc:212-284]

#### QuatroSolver
- 只用 src/dst 的前 2 行（XY 平面），求 SO(2) rotation
- 内部状态：`weights : 1 × Lp`、`rotation_2d : 2 × 2`
- 输出：`*rotation` 的左上 2×2 块设为 rotation_2d，其他保持 identity
- [VERIFY: teaser/src/registration.cc:286-414]

### 6.3 Inlier 提取

`solve()` 调完 rotation solver 后：

```cpp
for (size_t i = 0; i < rotation_inliers_mask_.cols(); ++i) {
    if (rotation_inliers_mask_[i]) {
        rotation_inliers_.emplace_back(i);
    }
}
// [VERIFY: teaser/src/registration.cc:718-722]
```

`rotation_inliers_` 是 **pruned TIM 维度的索引**，**不是原始测量索引**。要映回原始测量需要通过 `src_tims_map_rotation_`。

### 6.4 数据传递

- `solution_.rotation` → Stage 6 用于 `s · R · src`；
- `rotation_inliers_` 主要作调试/可视化（不参与 Stage 6 决策）。

---

## 7. Stage 6 — Translation Estimation

### 7.1 输入准备

```cpp
Eigen::Matrix<double, 3, Eigen::Dynamic> rotation_pruned_src(3, max_clique_.size());
Eigen::Matrix<double, 3, Eigen::Dynamic> rotation_pruned_dst(3, max_clique_.size());
for (size_t i = 0; i < max_clique_.size(); ++i) {
    rotation_pruned_src.col(i) = src.col(max_clique_[i]);
    rotation_pruned_dst.col(i) = dst.col(max_clique_[i]);
}
// [VERIFY: teaser/src/registration.cc:723-728]
```

**注意维度切换**：rotation 阶段用 TIM（`pruned_*_tims_`），但 translation 阶段重新切回**原始测量**（`max_clique_` 个），因为 translation 不能在 TIM 上估计（TIM 已经把 t 消掉了）。

### 7.2 调用

```cpp
solveForTranslation(solution_.scale * solution_.rotation * rotation_pruned_src,
                    rotation_pruned_dst);
// [VERIFY: teaser/src/registration.cc:732-733]
```

**变换**：把 `s · R · a_i` 当作 source，与 `b_i` 比对。如果估计正确，差值应近似 `t`：

```
b_i - s · R · a_i ≈ t + ε
```

### 7.3 TLSTranslationSolver 实现

```cpp
Eigen::Matrix<double, 3, Eigen::Dynamic> raw_translation = dst - src;
int N = src.cols();
double β = noise_bound_ * sqrt(cbar2_);
Eigen::Matrix<double, 1, Eigen::Dynamic> alphas = β * Eigen::MatrixXd::Ones(1, N);

*inliers = Ones(1, N);
Eigen::Matrix<bool, 1, Eigen::Dynamic> inliers_temp(1, N);
for (size_t i = 0; i < raw_translation.rows(); ++i) {
    tls_estimator_.estimate(raw_translation.row(i), alphas, &((*translation)(i)), &inliers_temp);
    *inliers = (*inliers).cwiseProduct(inliers_temp);
}
// [VERIFY: teaser/src/registration.cc:451-477]
```

**算子**：
1. `raw_translation = dst - src`（注意此时 `src` 已经是 `s · R · a_i`）；
2. 对每行（x、y、z 各一次）跑 `ScalarTLSEstimator::estimate()`；
3. 三轴 inlier mask 取 element-wise AND，得到最终 inlier。

**β 在 translation 中的取值**：`noise_bound · sqrt(cbar2)`——**没有因子 2**！与 scale 阶段的 `2β` 不同。

原因：scale 阶段处理 TIM（两点差），噪声放大 2 倍；translation 阶段处理单点（已经经过 `s · R` 变换），噪声只是单点噪声 `β`。

### 7.4 结果整合

```cpp
translation_inliers_ = utils::findNonzero<bool>(translation_inliers_mask_);
solution_.valid = true;
return solution_;
// [VERIFY: teaser/src/registration.cc:737-742]
```

`translation_inliers_` 是 **max-clique 维度的索引**（0 到 L−1）。

要拿到「**输入原始下标**」的 inlier，用 getter：

```cpp
std::vector<int> getInputOrderedTranslationInliers() {
    if (params_.rotation_estimation_algorithm == ROTATION_ESTIMATION_ALGORITHM::FGR) {
        throw std::runtime_error("Not supported when using FGR ...");
    }
    std::vector<int> translation_inliers;
    for (const auto& i : translation_inliers_) {
        translation_inliers.emplace_back(max_clique_[i]);
    }
    return translation_inliers;
}
// [VERIFY: teaser/include/teaser/registration.h:745-756]
```

---

## 8. 索引语义总结

不同阶段的 inlier mask 维度不同，索引语义也不同。下表总结：

| Mask / Index | 维度 | 索引指向 | 写入处 | 如何映回原始测量 |
|--------------|------|---------|--------|---------------|
| `scale_inliers_mask_` | 1×M (M=N(N-1)/2) | TIM 下标 | Stage 2 | TIM `c` ↔ 测量 `(src_tims_map_(0,c), src_tims_map_(1,c))` |
| `max_clique_` | L 个 int | 原始测量下标 | Stage 3 | 已经是原始下标 |
| `rotation_inliers_mask_` | 1×Lp | pruned TIM 下标 | Stage 5 | TIM 下标 `c` ↔ `(src_tims_map_rotation_(0,c), ..._rotation_(1,c))`（值已经是原始测量下标） |
| `rotation_inliers_` | 子集 of [0..Lp) | pruned TIM 下标 | Stage 5 | 同上 |
| `translation_inliers_mask_` | 1×L | max-clique 内下标 | Stage 6 | `i ↔ max_clique_[i]` |
| `translation_inliers_` | 子集 of [0..L) | max-clique 内下标 | Stage 6 | `i ↔ max_clique_[i]` |

### 8.1 工程坑点

```cpp
// ❌ 错误：直接用 rotation_inliers_ 索引原始 src
for (int idx : solver.getRotationInliers()) {
    auto p = src.col(idx);  // 错！idx 是 pruned TIM 索引，不是原始测量下标
}

// ✅ 正确：先经 dst_tims_map_rotation_ 反查
auto m = solver.getDstTIMsMapForRotation();
for (int idx : solver.getRotationInliers()) {
    int leaf_meas = m(0, idx);
    int root_meas = m(1, idx);
    // leaf_meas 和 root_meas 才是原始测量下标
}

// ✅ 正确：translation inlier 经 max_clique_ 映射
for (int input_idx : solver.getInputOrderedTranslationInliers()) {
    auto p = src.col(input_idx);  // 这是原始下标
}
```

`getInputOrderedTranslationInliers()` 是唯一一个**已经映回原始下标**的 getter。其他 inlier getter 都需要手动映射。

---

## 9. 多 solve 调用的状态影响

`solve()` 内部覆盖的字段（无 `reset()`）：

| 字段 | 每次 solve 是否覆盖 |
|------|-------------------|
| `solution_` | 是（每个 sub-solver 写入对应字段） |
| `src_tims_`、`dst_tims_` | 是（Stage 1 重赋值） |
| `*_tims_map_*` | 是 |
| `scale_inliers_mask_` 等 mask | 是（`resize` + 写入） |
| `pruned_*_tims_` | 是（resize + 赋值） |
| `max_clique_` | **NONE 模式下追加**！[VERIFY: teaser/src/registration.cc:656-659 使用 push_back]，其他模式是整体赋值 |
| `rotation_inliers_` | **追加**！[VERIFY: teaser/src/registration.cc:720 emplace_back]，每次 solve 累积下标 |
| `translation_inliers_` | 替换（`= utils::findNonzero<bool>(...)`）[VERIFY: teaser/src/registration.cc:737]，无累积 |
| `inlier_graph_` | **追加**！`populateVertices` 仅 resize，`addEdge` 不清空 [VERIFY: teaser/include/teaser/graph.h:67, 96-104] |

**重大隐患**：`max_clique_`、`rotation_inliers_`、`translation_inliers_`、`inlier_graph_` 在 `reset()` 之外**不会自动清空**。如果连续调用两次 `solve()` 而不 `reset()`，第二次的结果会**叠加**到第一次之上。

代码中的对策：构造函数路径 → `reset()` → 清空这些字段。但**调用者多次复用同一个 solver 实例时必须显式 `reset()`**。

最佳实践：每次 `solve()` 前先 `solver.reset(params)`，或每对点云用新的 solver 实例。

---

## 10. 数据流图（带维度标注）

```
Input (N=1000 correspondences)
 src : 3 × 1000          dst : 3 × 1000
   │                       │
   ├── computeTIMs() ──────┤
   ▼                       ▼
 src_tims_ : 3 × 499500   dst_tims_ : 3 × 499500
 src_tims_map_ : 2 × 499500
   │
   │      ┌────────── TLSScaleSolver / ScaleInliersSelector
   ▼      │
 raw_scales : 1 × 499500
 alphas    : 1 × 499500
   │
   │      ScalarTLSEstimator (sweep-line voting)
   ▼
 scale (scalar)
 scale_inliers_mask_ : 1 × 499500
   │
   │  ┌── populateVertices(1000) ────┐
   │  │                              ▼
   │  └── for each inlier TIM ──── addEdge(i, j)
   │                                  │
   ▼                                  ▼
 inlier_graph_ : 1000 vertices, ~M_inlier edges
   │
   │  MaxCliqueSolver::findMaxClique
   │  (PMC_EXACT / PMC_HEU / KCORE_HEU)
   ▼
 max_clique_ : sorted vector<int>, say L=200
   │
   ├── CHAIN  ────────► pruned_*_tims_ : 3 × 200
   │  (环状链)        src_tims_map_rotation_ : 2 × 200
   └── COMPLETE  ────► pruned_*_tims_ : 3 × 19900
                       (= L(L-1)/2 = 200·199/2)
   │
   │ ★ pruned_dst_tims_ /= scale
   │ ★ rotation_solver_.params.noise_bound *= 2/scale
   ▼
 GNCTLSRotationSolver / FGR / Quatro
   │
   ▼
 rotation : 3 × 3
 rotation_inliers_mask_ : 1 × Lp (Lp = 200 or 19900)
   │
   │
 rotation_pruned_src : 3 × 200    rotation_pruned_dst : 3 × 200
   │       (gather from src by max_clique_)
   ▼
 s · R · rotation_pruned_src : 3 × 200
   │
   │  TLSTranslationSolver (3 × ScalarTLSEstimator)
   ▼
 translation : 3 × 1
 translation_inliers_mask_ : 1 × 200
 translation_inliers_ : subset of [0..200)
   │
   ▼
 RegistrationSolution { valid=true, scale, rotation, translation }
```

---

## 11. 文档完整性检查表

- [x] 每个 Stage 的输入 / 算子 / 输出维度都有显式标注。
- [x] TIM 索引推导公式有数值验证。
- [x] scale 阶段 β 与 translation 阶段 β 的因子差异（2 vs 1）已解释。
- [x] 三种 inlier mask 的语义差异有专门的索引语义表。
- [x] 多次调用 solve 不 reset 的隐患（`max_clique_`、`inlier_graph_` 等追加）已标注。
- [x] `getInputOrderedTranslationInliers()` 是唯一映回原始下标的 getter 这一事实已点明。

