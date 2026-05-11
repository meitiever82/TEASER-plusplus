# ALGORITHM 02：GNC-TLS Rotation Estimation（图渐进式非凸 - 截断最小二乘旋转估计）

> **范围**：`teaser::GNCTLSRotationSolver`，TEASER++ 默认的 SO(3) 旋转估计器。本文档**不**展开 FGR 与 Quatro 实现（结构类似），只在对比表里点到为止。
>
> **核心代码**：
> - 算法实现 [VERIFY: teaser/src/registration.cc:770-872]
> - 类声明 [VERIFY: teaser/include/teaser/registration.h:257-278]
> - 父类 `GNCRotationSolver` [VERIFY: teaser/include/teaser/registration.h:218-247]
> - 加权 SVD 辅助函数 `svdRot` [VERIFY: teaser/include/teaser/utils.h:121-136]
>
> **理论来源**：
> - H. Yang, P. Antonante, V. Tzoumas, L. Carlone, "Graduated Non-Convexity for Robust Spatial Perception: From Non-Minimal Solvers to Global Outlier Rejection", RA-L 2020, arXiv:1909.08605
> - GNC（Graduated Non-Convexity）的 Black-Rangarajan 对偶展开。

---

## 1. 问题定义

### 1.1 输入

经过 max-clique 过滤的 TIM 对：
```
ã_i ∈ R³, b̃_i ∈ R³,  i = 1..N_p
```
其中 N_p = `pruned_src_tims_.cols()`，CHAIN 时 N_p = |max_clique_|，COMPLETE 时 N_p = |max_clique_|·(|max_clique_|−1)/2。

经过 stage 4 的尺度归一化 (`pruned_dst_tims_ /= scale`)，模型简化为：

```
b̃_i = R · ã_i + ε_i,    ‖ε_i‖ ≤ β'
```
其中 R ∈ SO(3)，β' = 2β/s 是 TIM 噪声上界（由调整后的 `params.noise_bound` 反映） [VERIFY: teaser/src/registration.cc:708-710]。

### 1.2 目标

求 R 与 inlier 集合 T 使 TLS 代价最小：

```
min_{R ∈ SO(3), T ⊆ {1..N_p}}
    Σ_{i ∈ T} ‖b̃_i − R·ã_i‖²  +  Σ_{i ∉ T} c̄² · β'²
```

`c̄² β'²` 是 outlier 单个罚分（实际代码里 `c̄²` 默认 1）。

### 1.3 为什么需要 GNC？

直接求解：
- T 的选择 → 离散变量，组合优化；
- R 在 SO(3) 上的优化 → 非凸（旋转的群结构）。

**Scalar TLS** 的 sweep-line 不适用：旋转不是 1D 量。
**RANSAC** 在 SO(3) 上效率低，且无最优性保证。

**GNC** 思路：把离散 T 替换为连续权重 `w_i ∈ [0,1]`，引入控制参数 μ，从 μ 很大（凸近似，权重接近 1）逐步推进到 μ → ∞（权重二值化，等价 TLS）。每个 μ 下问题对 (R, w) 都有简单结构。

---

## 2. Black-Rangarajan 对偶展开

### 2.1 TLS 代价的连续化

TLS 单点代价：
```
ρ(r²) = min(r², c̄²β'²)
      = { r²,        若 r² ≤ c̄²β'²
        { c̄²β'²,    若 r² > c̄²β'²
```

`ρ(r²)` 在 r² = c̄²β'² 处不光滑、不严格凸。直接最小化 `Σ ρ(r_i²)` 非常容易陷局部最优。

### 2.2 引入连续权重

Black-Rangarajan 定理：任何形如 `Σ ρ(r_i²)` 的非凸代价都可以重写为
```
Σ_i [ w_i · r_i² + Φ_ρ(w_i) ]
```
其中 `w_i ∈ [0, 1]`，`Φ_ρ` 是 ρ 的「对偶罚分」。**关键**：对固定 r²，对 w 求最优 `w*(r²)` 可以闭式算出来，且 `Σ ρ(r_i²) = min_w Σ [w_i r_i² + Φ_ρ(w_i)]`。

对 TLS：
```
Φ_TLS(w) = c̄²β'² · (w - 1) ·... [详细形式见论文]
```

直接最优 w*（TLS）：
```
w*_i(r²) = { 1,  若 r² ≤ c̄²β'²
           { 0,  否则
```
**二值化**——这就是 TLS 离散选择的等价连续形式。

### 2.3 GNC：渐进逼近

引入 control 参数 μ > 0，定义 ρ_μ（一族凸近似）：
```
ρ_μ(r²) → r²      as μ → 0
ρ_μ(r²) → ρ_TLS  as μ → ∞
```

对应的 w*_μ(r²) 也连续过渡：μ 小时接近全 1（凸），μ 大时接近 TLS 的 0/1。

TEASER 代码用的具体形式（**from registration.cc:834-849**）：

```cpp
double th1 = (mu + 1) / mu * noise_bound_sq;     // 上阈值
double th2 = mu / (mu + 1) * noise_bound_sq;     // 下阈值
for (size_t j = 0; j < match_size; ++j) {
    if (residuals_sq(j) >= th1) {
        weights(j) = 0;
    } else if (residuals_sq(j) <= th2) {
        weights(j) = 1;
    } else {
        weights(j) = sqrt(noise_bound_sq * mu * (mu + 1) / residuals_sq(j)) - mu;
    }
}
// [VERIFY: teaser/src/registration.cc:833-850]
```

写成数学形式：
```
                    ⎧ 0,                                          r² ≥ (μ+1)/μ · β²
w*_μ(r²) =          ⎨ √[β² μ(μ+1)/r²] − μ,                       β² μ/(μ+1) < r² < (μ+1)/μ · β²
                    ⎩ 1,                                          r² ≤ μ/(μ+1) · β²
```
其中 β² = `noise_bound_sq`。

### 2.4 三段函数的解读

- **下阈值 th2 = μ β²/(μ+1)**：r² 很小（残差小），完全 inlier；
- **上阈值 th1 = (μ+1)β²/μ**：r² 很大（残差大），完全 outlier；
- **中间过渡**：连续插值。

阈值随 μ 的演化：
| μ → 0 (极弱约束) | μ = 1 (平衡) | μ → ∞ (TLS) |
|----------------|-------------|-------------|
| th1 → ∞，th2 → 0：所有 w → 1 | th1 = 2β²，th2 = β²/2 | th1, th2 → β²：硬阈值 |

直观：μ 大时，过渡区间 [th2, th1] 收窄到 [β², β²]，权重函数变成阶跃；μ 小时过渡区间宽，权重函数平滑。

### 2.5 数学练习：验证三段连续

- 在 r² = th2 = μβ²/(μ+1) 处：

  ```
  √[β² μ(μ+1) / (μβ²/(μ+1))] − μ
  = √[μ(μ+1)·(μ+1)/μ] − μ
  = √[(μ+1)²] − μ
  = (μ+1) − μ = 1
  ```
  与 inlier 段一致。

- 在 r² = th1 = (μ+1)β²/μ 处：
  ```
  √[β² μ(μ+1) / ((μ+1)β²/μ)] − μ
  = √[μ·μ] − μ
  = μ − μ = 0
  ```
  与 outlier 段一致。

权重函数在端点连续。

### 2.6 中间段的来源

固定 R，对 w 求偏导（论文公式）：

```
∂/∂w_i [w_i r_i² + Φ_TLS(w_i; μ)] = 0
⇒ w*_μ(r²) = √[β² μ(μ+1) / r²] − μ   (中间区间内)
```

详细推导见论文 RA-L 2020 公式 (12)–(15)。

---

## 3. 算法主循环

### 3.1 整体伪代码

```
INPUT: src (3×N_p), dst (3×N_p), params
1. 初始化 R = I, weights = 1_N_p, μ = 1
2. for iter in 0..max_iterations:
     2.1 R = svdRot(src, dst, weights)        // 加权 SVD 闭式
     2.2 residuals_sq = colwise_sum((dst - R*src)²)
     2.3 if iter == 0:
            μ = 1 / (2 · max(residuals_sq)/β² − 1)
            if μ ≤ 0: BREAK
     2.4 prev_cost = cost
     2.5 cost = Σ weights[j] · residuals_sq[j]
     2.6 for j: weights[j] = w*_μ(residuals_sq[j])  // 三段函数
     2.7 μ = μ · gnc_factor
     2.8 if |cost − prev_cost| < cost_threshold: BREAK
3. 输出 R, inliers = (weights ≥ 0.5)
```

### 3.2 与代码逐行对照

完整代码 [VERIFY: teaser/src/registration.cc:770-872]：

#### 3.2.1 输入验证

```cpp
void teaser::GNCTLSRotationSolver::solveForRotation(
    const Eigen::Matrix<double, 3, Eigen::Dynamic>& src,
    const Eigen::Matrix<double, 3, Eigen::Dynamic>& dst, Eigen::Matrix3d* rotation,
    Eigen::Matrix<bool, 1, Eigen::Dynamic>* inliers) {
    assert(rotation);
    assert(src.cols() == dst.cols());
    assert(params_.gnc_factor > 1);
    assert(params_.noise_bound != 0);
    if (inliers) {
        assert(inliers->cols() == src.cols());
    }
// [VERIFY: teaser/src/registration.cc:770-780]
```

- `gnc_factor > 1`：μ 是递增的（与 FGR 相反），factor 必须真增大；
- `noise_bound != 0`：避免除零。

#### 3.2.2 初始化

```cpp
    size_t match_size = src.cols();
    double mu = 1; // arbitrary starting mu
    double prev_cost = std::numeric_limits<double>::infinity();
    cost_ = std::numeric_limits<double>::infinity();
    double noise_bound_sq = std::pow(params_.noise_bound, 2);
    if (noise_bound_sq < 1e-16) {
        noise_bound_sq = 1e-2;
    }
    Eigen::Matrix<double, 3, Eigen::Dynamic> diffs(3, match_size);
    Eigen::Matrix<double, 1, Eigen::Dynamic> weights(1, match_size);
    weights.setOnes(1, match_size);
    Eigen::Matrix<double, 1, Eigen::Dynamic> residuals_sq(1, match_size);
// [VERIFY: teaser/src/registration.cc:793-809]
```

- `mu = 1`：初始 μ。但**第一轮**会被覆盖（line 822-823）；
- `noise_bound_sq < 1e-16` 兜底为 `1e-2`：避免后续除零。这是个**重要的工程保护**——如果用户传 `noise_bound = 0`（虽 assert 拦截 `==0`，但接近 0 也危险），代码会自动用一个非零值。
- `weights.setOnes()`：初始权重全 1，对应「全 inlier」假设。

#### 3.2.3 第 1 步：固定 weights 解 R

```cpp
    for (size_t i = 0; i < params_.max_iterations; ++i) {
        *rotation = teaser::utils::svdRot(src, dst, weights);
// [VERIFY: teaser/src/registration.cc:812-815]
```

调 [VERIFY: teaser/include/teaser/utils.h:121-136]：

```cpp
inline Eigen::Matrix3d svdRot(const Eigen::Matrix<double, 3, Eigen::Dynamic>& X,
                              const Eigen::Matrix<double, 3, Eigen::Dynamic>& Y,
                              const Eigen::Matrix<double, 1, Eigen::Dynamic>& W) {
    Eigen::Matrix3d H = X * W.asDiagonal() * Y.transpose();    // H = Σ w_i X_i Y_i^T
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3d U = svd.matrixU();
    Eigen::Matrix3d V = svd.matrixV();
    if (U.determinant() * V.determinant() < 0) {
        V.col(2) *= -1;
    }
    return V * U.transpose();
}
```

这是经典的 Procrustes 问题加权解（Kabsch 算法的加权版本）：
```
H = X · diag(W) · Y^T = Σ_i w_i · x_i · y_i^T  (3×3)
SVD: H = U Σ V^T
R = V U^T,  with sign correction so det(R) = +1
```

数学根据：加权 Procrustes 求 R 最小化 `Σ w_i ‖y_i - R x_i‖²` = `Σ w_i (‖y_i‖² + ‖x_i‖² - 2 y_i^T R x_i)`。固定其他项，最大化 `Σ w_i y_i^T R x_i = trace(R H^T) = trace(V^T R U Σ)`。当 `R = V U^T` 时 trace = sum of singular values，最大。

**符号修正**：若 `det(U) det(V) < 0`，得到的 `V U^T` 是反射（det = -1），不在 SO(3) 中。把 V 最后一列翻号修正为旋转。这是 Kabsch 算法的标准做法。

#### 3.2.4 第 2 步：计算残差平方

```cpp
        diffs = (dst - (*rotation) * src).array().square();
        residuals_sq = diffs.colwise().sum();
// [VERIFY: teaser/src/registration.cc:817-819]
```

- `diffs(k, j) = (dst(k,j) - R*src(k,j))²`（element-wise）
- `residuals_sq(j) = ‖dst_j - R*src_j‖²`

#### 3.2.5 第 1 轮：初始化 μ

```cpp
        if (i == 0) {
            double max_residual = residuals_sq.maxCoeff();
            mu = 1 / (2 * max_residual / noise_bound_sq - 1);
            if (mu <= 0) {
                TEASER_DEBUG_INFO_MSG("GNC-TLS terminated because maximum residual at initialization is very small.");
                break;
            }
        }
// [VERIFY: teaser/src/registration.cc:820-831]
```

**μ 的初始化公式**：

```
μ₀ = 1 / (2 · max_r²/β² − 1)
```

为什么这个公式？

设 `r_max² = max_i r_i²`。希望初始 μ 让 `w*_μ(r_max²)` 处于「过渡边缘」——即恰好使 `r_max²` 落在 `th1 = (μ+1)β²/μ` 上。设此条件：
```
(μ+1)β²/μ = r_max²
⇒ μ + 1 = μ · r_max²/β²
⇒ μ = 1 / (r_max²/β² − 1)
⇒ μ = β² / (r_max² − β²)
```

但代码用的是 `μ₀ = 1 / (2·r_max²/β² − 1)`——多了个**因子 2**。

为什么是 2？这其实让初始 μ 更小（更「凸」）：
- 没有因子 2：`μ₀ = β²/(r_max² − β²)`，r_max² ≈ β² 时 μ₀ 很大，r_max² 远大于 β² 时 μ₀ 很小；
- 有因子 2：`μ₀ = β²/(2 r_max² − β²)`，更保守，永远比无因子 2 小（当 r_max² > β² 时）。

直观解释：用 r_max²/2 作为「典型」残差，μ 设得让典型残差刚好在过渡区中点。这个选择让初始权重整体接近 1（残差小的全是 inlier），同时给最大残差留一些「outlier 怀疑」。

**退出条件 μ ≤ 0**：当 `2·r_max²/β² ≤ 1`，即 `r_max² ≤ β²/2`，所有残差都小于 β²/2。此时已经全部是 inlier，TLS 退化为普通 LS——R 已经是闭式最优解，没必要继续 GNC。

#### 3.2.6 第 3 步：阈值计算与权重更新

```cpp
        double th1 = (mu + 1) / mu * noise_bound_sq;
        double th2 = mu / (mu + 1) * noise_bound_sq;
        cost_ = 0;
        for (size_t j = 0; j < match_size; ++j) {
            cost_ += weights(j) * residuals_sq(j);

            if (residuals_sq(j) >= th1) {
                weights(j) = 0;
            } else if (residuals_sq(j) <= th2) {
                weights(j) = 1;
            } else {
                weights(j) = sqrt(noise_bound_sq * mu * (mu + 1) / residuals_sq(j)) - mu;
                assert(weights(j) >= 0 && weights(j) <= 1);
            }
        }
// [VERIFY: teaser/src/registration.cc:834-850]
```

**注意 cost 计算的时序**：循环里 `cost_ += weights(j) * residuals_sq(j)` 用的是 **本轮更新前** 的 weights（上一轮留下来的）与**本轮**算出的 residuals_sq。所以 `cost_` 反映「上一轮的权重 + 本轮的 R」的代价。注释也明确说：
```cpp
// Note: the cost calculated is using the previously solved weights
// [VERIFY: teaser/src/registration.cc:838-839]
```

这是个**故意**的设计：用于和 `prev_cost` 做差，看「本轮 R 更新后的代价」是否相对「上一轮 R + 同样权重」改善了。

**`assert`**：中间段 `w*_μ(r²)` 数学上保证在 [0, 1] 内，但浮点误差可能溢出。这个 assert 在 Debug 下会捕获。

#### 3.2.7 第 4 步：μ 调整与收敛判定

```cpp
        double cost_diff = std::abs(cost_ - prev_cost);
        mu = mu * params_.gnc_factor;
        prev_cost = cost_;

        if (cost_diff < params_.cost_threshold) {
            TEASER_DEBUG_INFO_MSG("GNC-TLS solver terminated due to cost convergence.");
            break;
        }
    }
// [VERIFY: teaser/src/registration.cc:852-865]
```

**μ 的递增**：`mu *= gnc_factor`（默认 1.4，即每轮 μ × 1.4）。

**收敛条件**：连续两轮 cost 差 < cost_threshold（默认 `1e-6`）。

**注意顺序**：先算 `cost_diff`，再更新 μ，**再**判收敛。如果在 μ 已经更新到很大时（接近 TLS 极限）才检收敛，那 μ 已经走过了——但这无害，因为下一轮还会用新 μ 更新 weights 再算新 cost。

但有个**潜在问题**：如果在第 1 轮（`i==0`）就触发了 break（μ₀ ≤ 0），那 weights 还是全 1，R 是普通 LS 解。此时 `inliers` 是从 weights 推出来的，全部为 1（全 inlier）——合理，因为没有 outlier 怀疑。

#### 3.2.8 最终输出

```cpp
    if (inliers) {
        for (size_t i = 0; i < weights.cols(); ++i) {
            (*inliers)(0, i) = weights(0, i) >= 0.5;
        }
    }
}
// [VERIFY: teaser/src/registration.cc:867-871]
```

阈值 0.5：weights 是连续值，≥0.5 视为 inlier。这是论文形式：μ → ∞ 时 weights 会自然二值化到 0 或 1，0.5 是过渡区的中点。

**对比 FGR**：FGR 的 inlier 阈值是 `l_pq.cast<bool>()`——bool 转换是 `!= 0`，相当于阈值非零。差异：
- GNC-TLS：`w ≥ 0.5` 算 inlier，比较保守（要求权重明显偏 inlier 侧）；
- FGR：`l_pq != 0` 算 inlier，宽松（权重微小非零都算）。

---

## 4. svdRot 详细数学推导

### 4.1 加权 Procrustes 问题

求 R ∈ SO(3) 使

```
J(R) = Σ_i w_i · ‖y_i − R x_i‖²
```
最小。

### 4.2 展开

```
J(R) = Σ w_i [ ‖y_i‖² − 2 y_i^T R x_i + ‖x_i‖² ]
     = Σ w_i (‖y_i‖² + ‖x_i‖²)  − 2 · Σ w_i y_i^T R x_i
     = C − 2 · trace(R^T · H)         (常数 C 与 R 无关)
```
其中
```
H := Σ_i w_i · x_i · y_i^T
   = X · diag(w) · Y^T              (与代码一致)
```

### 4.3 最大化 trace(R^T H)

```
H = U Σ V^T   (SVD, σ_1 ≥ σ_2 ≥ σ_3 ≥ 0)
trace(R^T U Σ V^T) = trace(Σ V^T R^T U)
```
设 `M = V^T R^T U`，注意 R ∈ SO(3) 意味着 M 是正交阵（保 detM = ±1）：
- 若 det(U) det(V) > 0：`R = V U^T`，则 `M = V^T (V U^T)^T U = V^T U U^T U = V^T U U^T U`。

简化：直接由 `R = V U^T` 给出 `trace(R^T H) = trace(U V^T · U Σ V^T) = trace(Σ) = σ_1 + σ_2 + σ_3`，是最大可能值。

- 若 det(U) det(V) < 0：`V U^T` 行列式 −1，是反射不是旋转。需要把一个奇异向量翻号。代码翻 V 的第 3 列（对应最小奇异值 σ_3）：
  ```cpp
  if (U.determinant() * V.determinant() < 0) {
      V.col(2) *= -1;
  }
  ```
  这样 `R' = V' U^T` 满足 det(R') = +1，对应的 `trace(R'^T H) = σ_1 + σ_2 − σ_3`，是 SO(3) 内可达的最大值。

### 4.4 计算成本

`svdRot` 复杂度：H 计算 O(3·3·N_p) = O(N_p)，3×3 SVD 是 O(1)（固定大小）。所以 svdRot 本身 O(N_p)。

GNC-TLS 主循环每轮：svdRot O(N_p) + 残差计算 O(N_p) + 权重更新 O(N_p) = **O(N_p) per iter**。

总：O(max_iterations × N_p) ≈ O(100 N_p)。N_p = 200 时一次完整 GNC-TLS ~20000 浮点操作的 N_p 倍——毫秒级。

---

## 5. 关键设计选择

### 5.1 为什么 weights 是连续而不是 0/1？

直接 TLS 的 weights 应该是 0/1，但 0/1 的优化是组合的（NP-hard）。GNC 用连续松弛把组合优化转成「连续优化的极限」。

参考论文：Yang et al. RA-L 2020 证明了 GNC-TLS 在适当条件下**全局收敛到 TLS 全局最优**。这是 GNC-TLS 比 RANSAC 强的地方：
- RANSAC 是概率性的；
- GNC-TLS 是**deterministic** + **理论保证**。

### 5.2 为什么 μ 递增而不是递减？

代码 [VERIFY: registration.cc:856]：`mu = mu * params_.gnc_factor`，gnc_factor > 1，所以 μ 递增。

| Solver | μ 演化 | 原因 |
|--------|--------|------|
| **GNC-TLS** | μ ← μ · factor（递增） | 阈值函数 `1/(μ+1)μ`、`(μ+1)/μ` 形式：μ 大时阈值收敛到 β²（TLS 极限） |
| **FGR** | μ ← μ / factor（递减） | 不同的对偶形式：μ → 0 时阈值收紧到 TLS |

两者数学上等价（都从凸→TLS 过渡），只是参数化不同。

### 5.3 为什么 svdRot 的权重直接用 weights 而不是 weights · ranges？

在 1D Scalar-TLS 里，权重是 `1/α²`（结合 range）。但在 3D rotation 里，所有 TIM 共享同一个 noise_bound（因为它们的噪声同分布），所以权重直接用 GNC 给的连续 0–1 即可，不需要再乘 1/α²。

这是 1D 与 3D 在数据建模上的差异：1D 测量的 α 可以因测量而异（如 scale 中 α_i = β/‖ã_i‖），3D rotation 阶段假设所有 TIM 同质。

---

## 6. 复杂度与性能分析

### 6.1 时间复杂度

| 步骤 | 单轮成本 | 备注 |
|------|---------|------|
| svdRot | O(N_p) | 含 3×3 SVD 常数因子 |
| residuals_sq | O(N_p) | element-wise + colwise reduce |
| weights 更新 | O(N_p) | 三段函数 |
| cost 累加 | O(N_p) | 同上循环 |
| **每轮** | **O(N_p)** | |
| **总** | **O(max_iter · N_p)** | 默认 100 轮 |

### 6.2 收敛速度

实验观察（基于论文与代码注释）：
- 高 inlier 比例（>50%）：通常 10–30 轮收敛；
- 极高 outlier（<10% inlier）：50–100 轮，可能到 max_iterations 上限；
- 初始 μ 选得好时（公式给的）：前几轮 cost 下降快。

### 6.3 内存

主要矩阵：
- `weights`：1 × N_p
- `residuals_sq`：1 × N_p
- `diffs`：3 × N_p
- `R`：3 × 3
- `H`、`U`、`V`：3 × 3 各一

总：O(N_p)，可忽略。

---

## 7. ASCII 图示：μ 演化与权重函数

设 β² = 1.0：

```
μ = 0.1   (初始接近 0，凸近似)
  th2 = 0.1/1.1  = 0.091     th1 = 1.1/0.1 = 11.0
  权重函数：
   1 ┤████████████████░──────────
     │              ╲╲
  .5 ┤               ╲╲╲
     │                ╲╲╲╲
   0 ┤───────────────────╲╲░░░░░
     0          1          11        r²
  ★ 过渡区很宽（[0.091, 11.0]），大多数测量都得到中间权重

μ = 1.0   (中间，平衡)
  th2 = 0.5     th1 = 2.0
   1 ┤████░╲────
     │     ╲╲
  .5 ┤      ╲╲
     │       ╲╲
   0 ┤────────╲╲────
     0     1     2       r²

μ = 10    (大，接近 TLS)
  th2 = 0.909   th1 = 1.1
   1 ┤████████╲
     │         ╲
  .5 ┤          ╲
     │           ╲
   0 ┤───────────────
     0    1    1.1      r²
  ★ 过渡区很窄（[0.909, 1.1]），几乎二值化

μ = 100   (TLS 极限)
  th2 ≈ 0.99    th1 ≈ 1.01
  权重函数趋近阶跃：
   1 ┤████████│
            │
   0 ┤──────│──────
            β²            r²
```

---

## 8. 与 FGR / Quatro 的对比

### 8.1 三个 solver 的代码并排

| 维度 | GNC-TLS | FGR | Quatro |
|------|---------|-----|--------|
| 旋转维度 | SO(3) | SO(3) | SO(2)（仅 yaw） |
| 输入 src/dst 行数 | 3 | 3 | 用前 2 行（XY） |
| μ 演化 | 乘 gnc_factor（递增） | 除 gnc_factor（递减） | 乘 gnc_factor（同 GNC-TLS） |
| μ 初始化 | `1/(2 r_max²/β² − 1)` | `(diameter/β²)² · 1/β²` | `1/(2 r_max²/β² − 1)`（同 GNC-TLS） |
| 权重函数 | 三段：0 / 中间 / 1 | line process：`l = (μ̃/(μ̃+r²))²` | 三段（同 GNC-TLS） |
| 内层 R 解 | svdRot (3×3) | svdRot (3×3) | svdRot2d (2×2) |
| inlier 阈值 | `w ≥ 0.5` | `l > 0` | `w ≥ 0.4` |

代码位置：
- GNC-TLS：[VERIFY: teaser/src/registration.cc:770-872]
- FGR：[VERIFY: teaser/src/registration.cc:212-284]
- Quatro：[VERIFY: teaser/src/registration.cc:286-414]

### 8.2 行为差异

- **GNC-TLS** 最通用，收敛保证最强；
- **FGR** 在高 inlier 比例（>70%）下更快，但低 inlier 时容易陷局部；
- **Quatro** 只估 yaw，对城市/车载场景的「地平面已对齐」假设有用，能避免 SO(3) 完整估计在退化几何下的失败。

### 8.3 Quatro 的不寻常之处

Quatro 代码里有个 `static` 局部变量：

```cpp
static double rot_noise_bound = params_.noise_bound;
static double noise_bound_sq = std::pow(rot_noise_bound, 2);
// [VERIFY: teaser/src/registration.cc:335-336]
```

`static` 在 C++ 函数内意味着**只初始化一次**！第一次进函数时用 `params_.noise_bound`，之后即使 `setParams` 改了 noise_bound，这两个 static 变量**也不会更新**。

这是个**潜在 bug**——重复调用 Quatro 时若 noise_bound 改变了，行为会不一致。GNC-TLS 没有这个问题。

---

## 9. 工程坑点

### 9.1 noise_bound 已被调整

调用 `GNCTLSRotationSolver::solveForRotation` 前，主管线已经把 `params_.noise_bound *= 2/scale` [VERIFY: teaser/src/registration.cc:708-710]。所以 solver 内部看到的 noise_bound 不是用户传入的原始值。

调试时打印 `getParams().noise_bound` 看到的不是用户设置的——这是设计上的反映 TIM 噪声放大与 scale 归一化。

### 9.2 max_iterations 触顶不报错

代码在循环 100 轮没收敛后正常退出，不会抛异常或设置 `valid = false`。用户应额外检查 `getGNCRotationCostAtTermination()`：

```cpp
double final_cost = solver.getGNCRotationCostAtTermination();
if (final_cost > 1e-3 * N_p) {
    // 警告：rotation 可能没收敛
}
```

阈值 `1e-3 * N_p` 是经验值，没有理论保证。

### 9.3 第 1 轮 mu ≤ 0 break 时 weights 全 1

```cpp
if (mu <= 0) {
    break;
}
// [VERIFY: teaser/src/registration.cc:826-830]
```

此时 weights 还没经过 GNC 更新，全是初始的 1。所以 `inliers` 输出全 true——这是**故意**的，因为所有残差都很小（max_residual < β²/2），合理判全 inlier。

### 9.4 cost 累加时序

`cost_ += weights(j) * residuals_sq(j)` 用「上一轮的 weights」与「本轮的 residuals_sq」。这导致：
- 第 1 轮：weights 全 1，cost = Σ residuals_sq；
- 第 2 轮：weights 是第 1 轮 GNC 更新后的，cost ≈ Σ inlier_residuals_sq；
- ...

收敛判据 `|cost_ − prev_cost| < threshold` 实际比较的是「相邻两轮的 (上轮 weights · 本轮 R) 代价」。这不完全等价于「TLS 代价收敛」，但实际效果接近。

---

## 10. 关键不变量

| 不变量 | 验证 |
|--------|------|
| `R ∈ SO(3)`，即 det(R) = +1 | `svdRot` 的符号修正保证 [VERIFY: utils.h:131-133] |
| `weights ∈ [0, 1]^N_p` | 中间段 assert [VERIFY: registration.cc:848]，0/1 段显式赋值 |
| μ 单调递增 | `mu *= gnc_factor`，gnc_factor > 1 |
| `cost_` 与 `weights` 同步更新（同一循环） | 见 §3.2.6 |
| 最终 inliers = (weights ≥ 0.5) | [VERIFY: registration.cc:868-870] |

---

## 11. 与 Scalar-TLS 的算法学对照

| 维度 | Scalar-TLS | GNC-TLS Rotation |
|------|------------|-------------------|
| 决策变量 | x ∈ R | R ∈ SO(3) |
| 闭式解（固定 T） | 加权均值 | 加权 Procrustes (SVD) |
| 优化方法 | sweep-line（枚举有限 T） | GNC 渐进逼近 |
| 复杂度 | O(N log N) | O(max_iter · N) |
| 最优性 | 全局（论文形式下） | 全局（论文条件下） |
| 应用 | scale + translation 三轴 | rotation |

两者都是 TEASER 的「robust kernel」，只是不同维度的问题需要不同的解法。

---

## 12. 文档完整性检查表

- [x] Black-Rangarajan 对偶展开有具体公式（包括 TLS 的三段权重函数）。
- [x] 三段权重函数的端点连续性数学验证。
- [x] μ 初始化公式的推导（注释代码中因子 2 的几何意义）。
- [x] svdRot 的加权 Procrustes 完整推导（含符号修正）。
- [x] 三个 rotation solver 的对比表，并标注 Quatro 的 static 变量隐患。
- [x] 主循环每行有 [VERIFY:] 标签。
- [x] noise_bound 在主管线中被调整的事实（2/scale）已强调。
- [x] cost 累加时序的「上轮 weights × 本轮 R」语义点明。
- [x] ASCII 图示 μ 演化对权重函数形状的影响。

