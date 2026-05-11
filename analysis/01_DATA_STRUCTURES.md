# TEASER++ 关键数据结构

> **范围**：覆盖 registration pipeline 用到的全部 C++ 类型，包括公开 API、内部状态、参数枚举、第三方互操作面（PMC graph）。
>
> **阅读顺序建议**：先看 §1 几何基元（point cloud），再看 §2 主求解器（`RobustRegistrationSolver`）的字段，再看 §3 参数结构，最后看 §4 求解器接口族与 §5 图结构。

---

## 1. 几何基元

### 1.1 `teaser::PointXYZ`

最小的点表示：

```cpp
struct PointXYZ {
    float x;
    float y;
    float z;
    // operator==, operator!=
};
// [VERIFY: teaser/include/teaser/geometry.h:15-24]
```

注意 **`float` 而非 `double`**。这是公开输入接口的精度——一旦数据进 `solve()`，会被复制到 `Eigen::Matrix<double, 3, Dynamic>` 后全程双精度运算。  
[VERIFY: teaser/src/registration.cc:563-570]

`float` 输入的设计动机：与常见点云格式（PLY、PCD、传感器原始数据）的存储精度对齐，避免在 IO 层做不必要的 promotion。

### 1.2 `teaser::PointCloud`

包装 `std::vector<PointXYZ>` 的简单 STL-style 容器：

```cpp
class PointCloud {
public:
    using value_type      = PointXYZ;
    using reference       = PointXYZ&;
    using const_reference = const PointXYZ&;
    using iterator        = std::vector<PointXYZ>::iterator;
    using const_iterator  = std::vector<PointXYZ>::const_iterator;

    // capacity
    size_t size() const;
    void reserve(size_t n);
    bool empty();

    // element access
    PointXYZ& operator[](size_t i);
    PointXYZ& at(size_t n);
    PointXYZ& front();  PointXYZ& back();

    void push_back(const PointXYZ& pt);
    void clear();

private:
    std::vector<PointXYZ> points_;
};
// [VERIFY: teaser/include/teaser/geometry.h:26-70]
```

满足 C++ 命名容器要求（`value_type`、`iterator`、`size()`），所以可以直接被 range-for 遍历，也能配合 STL 算法。但没有 `random_access_iterator` 概念，因为它只暴露 `iterator` typedef 而没有 trait——实质上是顺序容器。

**用途**：只在 `solve(PointCloud&, PointCloud&, correspondences)` 的 overload 中作为输入参数。内部立刻被解包到 Eigen 矩阵，不再使用。所以**不要在性能敏感路径上反复构造 PointCloud**。  
[VERIFY: teaser/src/registration.cc:559-571]

---

## 2. 主求解器：`RobustRegistrationSolver`

### 2.1 公开接口面

类定义跨越 [VERIFY: teaser/include/teaser/registration.h:362-933]，分为 4 类成员：

1. **嵌套枚举与参数类型**（§3 详述）
2. **构造 / 重置**
3. **求解函数**（管线分段 + 总入口）
4. **结果访问器**（大量 inline getter）

#### 构造路径

```cpp
RobustRegistrationSolver();  // 默认构造，必须 reset() 后才能用
RobustRegistrationSolver(double noise_bound, double cbar2, bool estimate_scaling,
                         ROTATION_ESTIMATION_ALGORITHM rotation_estimation_algorithm,
                         double rotation_gnc_factor, size_t rotation_max_iterations,
                         double rotation_cost_threshold,
                         INLIER_GRAPH_FORMULATION rotation_tim_graph,
                         INLIER_SELECTION_MODE inlier_selection_mode,
                         double kcore_heuristic_threshold,
                         bool use_max_clique, bool max_clique_exact_solution,
                         double max_clique_time_limit, int max_clique_num_threads = 0);
// [VERIFY: teaser/include/teaser/registration.h:523-532]

RobustRegistrationSolver(const Params& params);  // 推荐
// [VERIFY: teaser/include/teaser/registration.h:541]
```

两个构造函数都最终调用 `reset(params)`：

```cpp
void reset(const Params& params) {
    params_ = params;
    if (params_.estimate_scaling) {
        setScaleEstimator(std::make_unique<TLSScaleSolver>(params_.noise_bound, params_.cbar2));
    } else {
        setScaleEstimator(std::make_unique<ScaleInliersSelector>(params_.noise_bound, params_.cbar2));
    }
    teaser::GNCRotationSolver::Params rotation_params{
        params_.rotation_max_iterations, params_.rotation_cost_threshold,
        params_.rotation_gnc_factor, params_.noise_bound};
    switch (params_.rotation_estimation_algorithm) {
        case ROTATION_ESTIMATION_ALGORITHM::GNC_TLS:
            setRotationEstimator(std::make_unique<GNCTLSRotationSolver>(rotation_params)); break;
        case ROTATION_ESTIMATION_ALGORITHM::FGR:
            setRotationEstimator(std::make_unique<FastGlobalRegistrationSolver>(rotation_params)); break;
        case ROTATION_ESTIMATION_ALGORITHM::QUATRO:
            setRotationEstimator(std::make_unique<QuatroSolver>(rotation_params)); break;
    }
    setTranslationEstimator(std::make_unique<TLSTranslationSolver>(params_.noise_bound, params_.cbar2));
    max_clique_.clear();
    rotation_inliers_.clear();
    translation_inliers_.clear();
    inlier_graph_.clear();
}
// [VERIFY: teaser/include/teaser/registration.h:823-863]
```

**关键观察**：每次 `reset()` 都会**重新构造**三个 solver 的 `unique_ptr`，这意味着前一次 `solve()` 留下的 solver 内部状态不会跨调用泄露——但 `solution_`、TIM 矩阵、`scale_inliers_mask_` 等并不在 `reset()` 里清除（除了显式 `clear()` 的 `max_clique_` 等）。后续 `solve()` 会就地覆盖这些字段，所以单线程下安全；**多线程下不能共享一个 `RobustRegistrationSolver` 实例**。

### 2.2 私有字段（内部状态）

完整定义见 [VERIFY: teaser/include/teaser/registration.h:892-933]：

```cpp
private:
    Params params_;
    RegistrationSolution solution_;

    // Stage 2 / 5 / 6 输出的 inlier mask
    Eigen::Matrix<bool, 1, Eigen::Dynamic> scale_inliers_mask_;
    Eigen::Matrix<bool, 1, Eigen::Dynamic> rotation_inliers_mask_;
    Eigen::Matrix<bool, 1, Eigen::Dynamic> translation_inliers_mask_;

    // Stage 1 全量 TIM（3 × N(N-1)/2）
    Eigen::Matrix<double, 3, Eigen::Dynamic> src_tims_;
    Eigen::Matrix<double, 3, Eigen::Dynamic> dst_tims_;

    // Stage 4 重建的 TIM（rotation solver 输入）
    Eigen::Matrix<double, 3, Eigen::Dynamic> pruned_src_tims_;
    Eigen::Matrix<double, 3, Eigen::Dynamic> pruned_dst_tims_;

    // TIM 索引映射：2 × M，每列 [i, j] 表示 ã_col = a_i - a_j
    Eigen::Matrix<int, 2, Eigen::Dynamic> src_tims_map_;          // Stage 1 全量
    Eigen::Matrix<int, 2, Eigen::Dynamic> dst_tims_map_;          // 等价 src_tims_map_，构造时分开存
    Eigen::Matrix<int, 2, Eigen::Dynamic> src_tims_map_rotation_; // Stage 4 CHAIN/COMPLETE
    Eigen::Matrix<int, 2, Eigen::Dynamic> dst_tims_map_rotation_;

    // Stage 3 输出（原始测量索引）
    std::vector<int> max_clique_;

    // Stage 5 输出（pruned TIM 维度的 inlier 索引）
    std::vector<int> rotation_inliers_;

    // Stage 6 输出（max-clique 维度的 inlier 索引，最终结果）
    std::vector<int> translation_inliers_;

    // Stage 3 用的 inlier graph
    teaser::Graph inlier_graph_;

    // Stage 2 / 5 / 6 注入式 solver（pimpl-ish）
    std::unique_ptr<AbstractScaleSolver> scale_solver_;
    std::unique_ptr<GNCRotationSolver> rotation_solver_;
    std::unique_ptr<AbstractTranslationSolver> translation_solver_;
```

#### 字段语义索引表

| 字段 | 维度 | 写入处 | 读取处 |
|------|------|--------|--------|
| `solution_` | scalar+9 doubles | `solveForScale/Rotation/Translation`、`solve()` 末尾 | `getSolution()` |
| `scale_inliers_mask_` | `1 × M` (M=N(N-1)/2) | `solveForScale()` | `getScaleInliersMask()`、Stage 3 建图 |
| `rotation_inliers_mask_` | `1 × L_pruned` | `solveForRotation()` | `getRotationInliersMask()`、`rotation_inliers_` 提取 |
| `translation_inliers_mask_` | `1 × L_pruned` | `solveForTranslation()` | `getTranslationInliersMask()` |
| `src_tims_` / `dst_tims_` | `3 × M` | `solve()` 中 `computeTIMs()` | Stage 2 输入 |
| `pruned_src_tims_` / `pruned_dst_tims_` | `3 × L_pruned` | Stage 4 | Stage 5 输入 |
| `src_tims_map_` / `dst_tims_map_` | `2 × M` | `solve()` 中 `computeTIMs()` | Stage 3 建图（取 `src_tims_map_(0,i)` 与 `src_tims_map_(1,i)`） |
| `src_tims_map_rotation_` / `dst_tims_map_rotation_` | `2 × L_pruned` | Stage 4 | 仅 getter 暴露 |
| `max_clique_` | `vector<int>`，大小 `L = |C|` | Stage 3 | Stage 4、Stage 6 |
| `rotation_inliers_` | `vector<int>` | Stage 5 末尾 | `getRotationInliers()` |
| `translation_inliers_` | `vector<int>` | Stage 6 末尾 | `getTranslationInliers()` |
| `inlier_graph_` | `Graph`，N 顶点 | Stage 3 入口 | `findMaxClique` 内部 |

注意三个 mask 的**维度不同**：
- `scale_inliers_mask_` 是 **TIM 维度** `M`；
- `rotation_inliers_mask_` 是 **pruned TIM 维度** `L_pruned`（CHAIN 时为 `L`，COMPLETE 时为 `L(L-1)/2`）；
- `translation_inliers_mask_` 是 **max-clique 维度** `L`（每个 max-clique 顶点对应一个原始测量）。

这是工程上很容易踩坑的地方——三种 mask 不能直接做 logical AND，必须先经 `dst_tims_map_` 等索引映射回原始测量。

#### 维度推导举例

设输入 `N = 1000` 个对应：

```
Stage 1: M = N(N-1)/2 = 499,500  TIM 列数
         src_tims_, dst_tims_     : 3 × 499500
         src_tims_map_, dst_tims_map_ : 2 × 499500
Stage 2: scale_inliers_mask_       : 1 × 499500
Stage 3: max_clique_               : 假设 L = 200
Stage 4 (CHAIN):
         pruned_src_tims_          : 3 × 200
         pruned_dst_tims_          : 3 × 200
         src_tims_map_rotation_    : 2 × 200
Stage 5: rotation_inliers_mask_    : 1 × 200
Stage 6: translation_inliers_mask_ : 1 × 200
         final inlier 数 ≤ 200
```

Stage 1 的内存占用是 `O(N²)`，N=10000 时单个矩阵就是 ~3.6 GB，所以**实务上 N 不宜 > 10k**——这是 TEASER 算法本身的瓶颈，与实现无关。

### 2.3 `RegistrationSolution`

最终输出：

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

**`valid` 字段的语义**（基于代码而非文档）：
- 构造时默认 `true`；
- max-clique 大小 ≤ 1 时被显式置为 `false` 并提前返回 [VERIFY: teaser/src/registration.cc:649-653]；
- 成功路径末尾再次显式赋 `true` [VERIFY: teaser/src/registration.cc:740-742]；
- **其他失败模式**（GNC 不收敛、scale 退化）**不会**自动置 `false`。

**`EIGEN_MAKE_ALIGNED_OPERATOR_NEW`** 是因为含 `Eigen::Matrix3d`，按 16 字节对齐分配 placement new 以利 SIMD。如果用户在堆上 `new RegistrationSolution`，必须使用 Eigen 提供的 aligned new。栈上声明或作为成员都没问题。

---

## 3. 参数结构体族

### 3.1 `RobustRegistrationSolver::Params`

主参数结构，14 个字段：

```cpp
struct Params {
    double noise_bound = 0.01;
    double cbar2 = 1;
    bool   estimate_scaling = true;
    ROTATION_ESTIMATION_ALGORITHM rotation_estimation_algorithm = GNC_TLS;
    double rotation_gnc_factor = 1.4;
    size_t rotation_max_iterations = 100;
    double rotation_cost_threshold = 1e-6;
    INLIER_GRAPH_FORMULATION rotation_tim_graph = CHAIN;
    INLIER_SELECTION_MODE inlier_selection_mode = PMC_EXACT;
    double kcore_heuristic_threshold = 0.5;
    bool   use_max_clique = true;            // deprecated
    bool   max_clique_exact_solution = true; // deprecated
    double max_clique_time_limit = 3600;
    int    max_clique_num_threads = teaser_default_max_threads();
};
// [VERIFY: teaser/include/teaser/registration.h:420-514]
```

#### 字段语义与使用位置

| 字段 | 使用位置 | 备注 |
|------|----------|------|
| `noise_bound` | TLS 阈值 `β`，所有 solver 都用 | 单位与点云一致，**最重要的参数** |
| `cbar2` | TLS 阈值平方比 `c̄²` | 通常 1 |
| `estimate_scaling` | `reset()` 决定 scale solver 类型 | false 时跳 TLS scale，只筛 inlier |
| `rotation_estimation_algorithm` | `reset()` 决定 rotation solver 类型 | 3 选 1 |
| `rotation_gnc_factor` | GNC-TLS：μ ← μ·factor；FGR：μ ← μ/factor | 必须 > 1 |
| `rotation_max_iterations` | rotation solver 外层循环上限 | |
| `rotation_cost_threshold` | rotation solver 收敛阈值 | GNC: cost 差；FGR: cost 绝对值 |
| `rotation_tim_graph` | Stage 4 决定 CHAIN 还是 COMPLETE | CHAIN 默认 |
| `inlier_selection_mode` | Stage 3 决定 max-clique 算法 | PMC_EXACT / PMC_HEU / KCORE_HEU / NONE |
| `kcore_heuristic_threshold` | Stage 3 KCORE_HEU 触发阈值 | 仅在 inlier_selection_mode = KCORE_HEU 时生效 |
| `max_clique_time_limit` | PMC 求解器超时 (秒) | 默认 1 小时 |
| `max_clique_num_threads` | PMC 并行线程数 | 默认 OpenMP 最大线程 |
| `use_max_clique` | `solve()` 开头，false 时强制 inlier_selection_mode = NONE | **deprecated** |
| `max_clique_exact_solution` | `solve()` 开头，false 时强制 inlier_selection_mode = PMC_HEU | **deprecated** |

废弃字段的兼容逻辑：
```cpp
if (!params_.use_max_clique) {
    params_.inlier_selection_mode = INLIER_SELECTION_MODE::NONE;
}
if (!params_.max_clique_exact_solution) {
    params_.inlier_selection_mode = INLIER_SELECTION_MODE::PMC_HEU;
}
// [VERIFY: teaser/src/registration.cc:579-589]
```

**注意优先级**：先检查 `use_max_clique`，再检查 `max_clique_exact_solution`。所以：

| `use_max_clique` | `max_clique_exact_solution` | 实际 `inlier_selection_mode` |
|------------------|-----------------------------|--------------------------------|
| true | true | 用户显式设置的值 |
| true | false | PMC_HEU |
| false | * | NONE（覆盖第二个分支） |

### 3.2 三个嵌套枚举

#### `ROTATION_ESTIMATION_ALGORITHM`

```cpp
enum class ROTATION_ESTIMATION_ALGORITHM {
    GNC_TLS = 0,
    FGR     = 1,
    QUATRO  = 2,
};
// [VERIFY: teaser/include/teaser/registration.h:383-387]
```

只在 `reset()` 的 switch 里被消费，决定 `rotation_solver_` 的具体类型。

#### `INLIER_SELECTION_MODE`

```cpp
enum class INLIER_SELECTION_MODE {
    PMC_EXACT = 0,
    PMC_HEU   = 1,
    KCORE_HEU = 2,
    NONE      = 3,
};
// [VERIFY: teaser/include/teaser/registration.h:397-402]
```

| 值 | 行为 | 选择时机 |
|----|------|---------|
| `PMC_EXACT` | PMC 精确算法，含 k-core 剪枝 + 邻居核排序 + 动态着色 | N ≤ 几千且 outlier 比例高，鲁棒性优先 |
| `PMC_HEU`   | PMC 启发式（`pmc_heu::search`） | 速度优先，inlier 比例较高 |
| `KCORE_HEU` | 直接取 max-core 全部顶点 | 大规模，可接受宽松 inlier |
| `NONE`      | 跳过 max-clique | 已知 inlier 几乎全是 inlier 或想绕过此阶段做对比 |

PMC_EXACT / PMC_HEU 之间还有个奇怪的中间状态：当 `solver_mode = KCORE_HEU` **且** `max_core > threshold·N` 时才真的走 k-core 启发式路径；否则会继续走 PMC_HEU 的下界搜索，并在 `lb < ub` 时**不**调用 exact——结果就是返回启发式找到的 clique。  
[VERIFY: teaser/src/graph.cc:66-122]

具体策略选择详见 `ALGORITHM_03-MaxClique_TIM.md`。

#### `INLIER_GRAPH_FORMULATION`

```cpp
enum class INLIER_GRAPH_FORMULATION {
    CHAIN    = 0,
    COMPLETE = 1,
};
// [VERIFY: teaser/include/teaser/registration.h:410-413]
```

Stage 4 用：决定从 max-clique 重新生成 TIM 的拓扑。

| 值 | TIM 数 | 实现方式 |
|----|--------|----------|
| `CHAIN` | `L = |max_clique_|` | 环：第 i 个连接到 (i+1) % L，最后一个连回第 0 个 |
| `COMPLETE` | `L(L-1)/2` | 重新跑 `computeTIMs()` |

代码：[VERIFY: teaser/src/registration.cc:663-700]。CHAIN 实现明确显示了环结构：

```cpp
for (size_t i = 0; i < max_clique_.size(); ++i) {
    const auto& root = max_clique_[i];
    int leaf;
    if (i != max_clique_.size() - 1) {
        leaf = max_clique_[i + 1];
    } else {
        leaf = max_clique_[0];  // 环回到第 0 个
    }
    pruned_src_tims_.col(i) = src.col(leaf) - src.col(root);
    pruned_dst_tims_.col(i) = dst.col(leaf) - dst.col(root);
    ...
}
// [VERIFY: teaser/src/registration.cc:670-686]
```

### 3.3 `GNCRotationSolver::Params`

子类共享的 GNC 参数：

```cpp
struct Params {
    size_t max_iterations;
    double cost_threshold;
    double gnc_factor;
    double noise_bound;
};
// [VERIFY: teaser/include/teaser/registration.h:223-228]
```

由 `RobustRegistrationSolver::reset()` 用 main `Params` 填充：

```cpp
teaser::GNCRotationSolver::Params rotation_params{
    params_.rotation_max_iterations,
    params_.rotation_cost_threshold,
    params_.rotation_gnc_factor,
    params_.noise_bound
};
// [VERIFY: teaser/include/teaser/registration.h:835-837]
```

注意 `noise_bound` 在 `solve()` 进 rotation 之前会被**临时改写**：

```cpp
auto params = rotation_solver_->getParams();
params.noise_bound *= (2 / solution_.scale);
rotation_solver_->setParams(params);
// [VERIFY: teaser/src/registration.cc:708-710]
```

含义：`pruned_dst_tims_ /= scale` 之后（[registration.cc:703]），dst 的尺度变小（scale 倍），但 src 已经是「未缩放」的。TIM 的差向量噪声本身是两点噪声差，所以 noise bound 变成 `2 * β`（因子 2），再除以 scale 把它映射到「dst 已缩放后的」坐标系。这是 GNC-TLS / FGR / Quatro 收到的真实 noise bound。

### 3.4 `MaxCliqueSolver::Params`

```cpp
struct Params {
    CLIQUE_SOLVER_MODE solver_mode = PMC_EXACT;
    bool   solve_exactly = true;            // deprecated
    double kcore_heuristic_threshold = 1;
    double time_limit = 3600;
    int    num_threads = 1;
};
// [VERIFY: teaser/include/teaser/graph.h:233-262]
```

由 `solve()` 从 main `Params` 派发：

```cpp
teaser::MaxCliqueSolver::Params clique_params;
if (params_.inlier_selection_mode == INLIER_SELECTION_MODE::PMC_EXACT) {
    clique_params.solver_mode = CLIQUE_SOLVER_MODE::PMC_EXACT;
} else if (params_.inlier_selection_mode == INLIER_SELECTION_MODE::PMC_HEU) {
    clique_params.solver_mode = CLIQUE_SOLVER_MODE::PMC_HEU;
} else {
    clique_params.solver_mode = CLIQUE_SOLVER_MODE::KCORE_HEU;
}
clique_params.time_limit = params_.max_clique_time_limit;
clique_params.kcore_heuristic_threshold = params_.kcore_heuristic_threshold;
clique_params.num_threads = params_.max_clique_num_threads;
// [VERIFY: teaser/src/registration.cc:627-638]
```

**注意默认值差异**：
- `MaxCliqueSolver::Params::kcore_heuristic_threshold` 默认 `1`（意为「只有 max-core 真等于全部顶点时才走 k-core」，几乎不触发）；
- `RobustRegistrationSolver::Params::kcore_heuristic_threshold` 默认 `0.5`（中等触发）。

这俩字段同名但默认不同。主求解器**会显式把自己的值复制给** `MaxCliqueSolver::Params`（registration.cc:637），所以最终生效的是主求解器的值。

`solver_mode` 与已废弃的 `solve_exactly`：
```cpp
if (!params_.solve_exactly) {
    params_.solver_mode = CLIQUE_SOLVER_MODE::PMC_HEU;
}
// [VERIFY: teaser/src/graph.cc:14-17]
```

### 3.5 `MaxCliqueSolver::CLIQUE_SOLVER_MODE`

```cpp
enum class CLIQUE_SOLVER_MODE {
    PMC_EXACT = 0,
    PMC_HEU   = 1,
    KCORE_HEU = 2,
};
// [VERIFY: teaser/include/teaser/graph.h:224-228]
```

**与 `RobustRegistrationSolver::INLIER_SELECTION_MODE` 区别**：前者只有 3 个值（没有 `NONE`），因为 `NONE` 是在更高层（`solve()`）被处理的——`inlier_selection_mode == NONE` 时**根本不会构造** `MaxCliqueSolver`，直接把 `max_clique_` 填成 `0..N-1`。  
[VERIFY: teaser/src/registration.cc:654-660]

---

## 4. 求解器接口族（abstract base classes）

设计模式：**Strategy 模式**。主求解器持有三个 unique_ptr 指向抽象基类，具体实现可替换。每个抽象基类只有 1 个纯虚函数。

### 4.1 `AbstractScaleSolver`

```cpp
class AbstractScaleSolver {
public:
    virtual ~AbstractScaleSolver() {}
    virtual void solveForScale(
        const Eigen::Matrix<double, 3, Eigen::Dynamic>& src,
        const Eigen::Matrix<double, 3, Eigen::Dynamic>& dst,
        double* scale,
        Eigen::Matrix<bool, 1, Eigen::Dynamic>* inliers) = 0;
};
// [VERIFY: teaser/include/teaser/registration.h:44-58]
```

#### 实现 1: `TLSScaleSolver`

```cpp
class TLSScaleSolver : public AbstractScaleSolver {
    explicit TLSScaleSolver(double noise_bound, double cbar2);
    void solveForScale(...) override;
private:
    double noise_bound_;
    double cbar2_;
    ScalarTLSEstimator tls_estimator_;
};
// [VERIFY: teaser/include/teaser/registration.h:136-160]
```

实现 [VERIFY: teaser/src/registration.cc:416-431]：把 TIM 的范数比作为 1D 测量，调 `ScalarTLSEstimator::estimate()`：

```cpp
v1_dist = ‖src.col(i)‖     // = ‖ã_i‖
v2_dist = ‖dst.col(i)‖     // = ‖b̃_i‖
raw_scales = v2_dist ./ v1_dist
β = 2 * noise_bound_ * sqrt(cbar2_)
alphas = β ./ v1_dist
tls_estimator_.estimate(raw_scales, alphas, &scale, inliers)
```

#### 实现 2: `ScaleInliersSelector`

```cpp
class ScaleInliersSelector : public AbstractScaleSolver { ... };
// [VERIFY: teaser/include/teaser/registration.h:167-187]
```

实现 [VERIFY: teaser/src/registration.cc:433-449]：固定 `*scale = 1`，仅做 inlier 筛选：

```cpp
β = 2 * noise_bound_ * sqrt(cbar2_)
*inliers = |‖src‖ - ‖dst‖| ≤ β
```

公式：当 `s = 1` 时，inlier TIM 满足 `‖b̃ − ã‖ ≤ 2β`，再由三角不等式 `|‖b̃‖ − ‖ã‖| ≤ ‖b̃ − ã‖`，可得保守阈值 `2β`。

### 4.2 `AbstractRotationSolver`

```cpp
class AbstractRotationSolver {
public:
    virtual ~AbstractRotationSolver() {}
    virtual void solveForRotation(
        const Eigen::Matrix<double, 3, Eigen::Dynamic>& src,
        const Eigen::Matrix<double, 3, Eigen::Dynamic>& dst,
        Eigen::Matrix3d* rotation,
        Eigen::Matrix<bool, 1, Eigen::Dynamic>* inliers) = 0;
};
// [VERIFY: teaser/include/teaser/registration.h:64-79]
```

#### `GNCRotationSolver`（中间层）

```cpp
class GNCRotationSolver : public AbstractRotationSolver {
public:
    struct Params { size_t max_iterations; double cost_threshold; double gnc_factor; double noise_bound; };
    GNCRotationSolver(Params params) : params_(params) {}
    Params getParams() { return params_; }
    void   setParams(Params params) { params_ = params; }
    double getCostAtTermination() { return cost_; }
protected:
    Params params_;
    double cost_;
};
// [VERIFY: teaser/include/teaser/registration.h:218-247]
```

三个具体 rotation 实现都继承 `GNCRotationSolver`：

| 实现 | 头位置 | 行为 |
|------|--------|------|
| `GNCTLSRotationSolver` | registration.h:257-278 | 完整 SO(3) GNC-TLS |
| `FastGlobalRegistrationSolver` | registration.h:290-322 | line-process FGR（μ 递减） |
| `QuatroSolver` | registration.h:329-355 | 仅 SO(2) yaw |

### 4.3 `AbstractTranslationSolver`

```cpp
class AbstractTranslationSolver {
public:
    virtual ~AbstractTranslationSolver() {}
    virtual void solveForTranslation(
        const Eigen::Matrix<double, 3, Eigen::Dynamic>& src,
        const Eigen::Matrix<double, 3, Eigen::Dynamic>& dst,
        Eigen::Vector3d* translation,
        Eigen::Matrix<bool, 1, Eigen::Dynamic>* inliers) = 0;
};
// [VERIFY: teaser/include/teaser/registration.h:85-100]
```

只有一个实现 `TLSTranslationSolver`，三个坐标轴各跑一次 `ScalarTLSEstimator::estimate()`：

```cpp
*inliers = Ones(1, N);
for (i = 0..2) {
    tls_estimator_.estimate(raw_translation.row(i), alphas, &((*translation)(i)), &inliers_temp);
    *inliers = (*inliers).cwiseProduct(inliers_temp);  // logical AND
}
// [VERIFY: teaser/src/registration.cc:469-476]
```

`cwiseProduct` 在 `bool` 上等价于 element-wise AND——保证 x、y、z 三轴**都**是 inlier 的点才算 final inlier。

### 4.4 `ScalarTLSEstimator`

不继承任何 abstract base，是个独立的工具类，被 `TLSScaleSolver` 和 `TLSTranslationSolver` 拥有为成员。

```cpp
class ScalarTLSEstimator {
public:
    ScalarTLSEstimator() = default;
    void estimate(const Eigen::RowVectorXd& X, const Eigen::RowVectorXd& ranges,
                  double* estimate, Eigen::Matrix<bool, 1, Eigen::Dynamic>* inliers);
    void estimate_tiled(const Eigen::RowVectorXd& X, const Eigen::RowVectorXd& ranges,
                        const int& s, double* estimate, Eigen::Matrix<bool, 1, Eigen::Dynamic>* inliers);
};
// [VERIFY: teaser/include/teaser/registration.h:105-131]
```

`estimate()` 是 sweep-line adaptive voting 实现（详见 `ALGORITHM_01-Scalar_TLS.md`），`estimate_tiled()` 是 loop-tiling + OpenMP 并行版本，理论结果等价、量级大时更快。

调用面：
- `TLSScaleSolver::solveForScale()` → `estimate()` [VERIFY: teaser/src/registration.cc:430]
- `TLSTranslationSolver::solveForTranslation()` → `estimate()` 三次 [VERIFY: teaser/src/registration.cc:472]

代码中**没有调用 `estimate_tiled()` 的地方**（除非用户绕过封装直接调用）——它是个备用 API，可能用于 benchmark 或第三方扩展。

---

## 5. 图结构

### 5.1 `teaser::Graph`

简单邻接表无向图：

```cpp
class Graph {
public:
    Graph() : num_edges_(0) {}
    explicit Graph(const std::map<int, std::vector<int>>& adj_list);

    void addVertex(int id);
    void populateVertices(int num_vertices);
    bool hasEdge(int v1, int v2);
    bool hasVertex(int v);
    void addEdge(int v1, int v2);
    void removeEdge(int v1, int v2);
    int  numVertices() const { return adj_list_.size(); }
    int  numEdges() const { return num_edges_; }
    const std::vector<int>& getEdges(int id) const;
    std::vector<int> getVertices() const;
    Eigen::MatrixXi getAdjMatrix() const;
    std::vector<std::vector<int>> getAdjList() const;
    void reserve(int num_vertices);
    void clear();
    void reserveForCompleteGraph(int num_vertices);

private:
    std::vector<std::vector<int>> adj_list_;
    size_t num_edges_;
};
// [VERIFY: teaser/include/teaser/graph.h:29-207]
```

**约束**：顶点 ID 必须从 0 连续到 `numVertices()-1`。`addVertex(id)` 会把 `adj_list_` resize 到 `id+1`，所以中间不能有空缺。  
[VERIFY: teaser/include/teaser/graph.h:55-61]

**`addEdge` 维护对称**：会同时把 v1 加到 v2 的邻接列表和 v2 加到 v1 的邻接列表 [VERIFY: teaser/include/teaser/graph.h:96-104]。所以 `num_edges_` 表示无向边数，而 `adj_list_` 里每条无向边出现两次。

#### 在管线中的角色

Stage 3 构造 inlier graph：

```cpp
inlier_graph_.populateVertices(src.cols());                    // N 个顶点
for (i = 0..scale_inliers_mask_.cols()-1) {
    if (scale_inliers_mask_(0, i)) {
        inlier_graph_.addEdge(src_tims_map_(0, i), src_tims_map_(1, i));
    }
}
// [VERIFY: teaser/src/registration.cc:620-625]
```

注意 `src_tims_map_(k, i)` 取第 i 个 TIM 的第 k 个端点（k ∈ {0, 1}）。`src_tims_map_` 和 `dst_tims_map_` 在 `computeTIMs()` 内部其实是用相同方式构造的，所以等价。

### 5.2 `teaser::MaxCliqueSolver`

```cpp
class MaxCliqueSolver {
public:
    enum class CLIQUE_SOLVER_MODE { PMC_EXACT, PMC_HEU, KCORE_HEU };
    struct Params { ... };  // 见 §3.4

    MaxCliqueSolver() = default;
    MaxCliqueSolver(Params params);

    std::vector<int> findMaxClique(Graph graph);

private:
    Graph graph_;
    Params params_;
};
// [VERIFY: teaser/include/teaser/graph.h:219-279]
```

**注意 `findMaxClique` 按值接收 `Graph`**——意味着调用时图会被拷贝一次。对 inlier graph 不是问题（邻接表小），但调用者不应假设 `graph_` 被存进了 solver 实例。

返回类型是 `vector<int>`：max-clique 的**顶点 ID 列表**，未排序。`solve()` 收到后会显式排序：

```cpp
max_clique_ = clique_solver.findMaxClique(inlier_graph_);
std::sort(max_clique_.begin(), max_clique_.end());
// [VERIFY: teaser/src/registration.cc:641-642]
```

排序的目的是让 CHAIN 阶段的「环连接」按数值升序遍历，避免顺序歧义。

#### PMC interop

`findMaxClique` 内部把 `teaser::Graph` 转换成 `pmc::pmc_graph` 的 CSR-like 表示：

```cpp
std::vector<int> edges;            // 邻接列表展平
std::vector<long long> vertices;   // 每个顶点在 edges 中的起始偏移
vertices.push_back(0);
for (i : all_vertices) {
    edges.insert(edges.end(), c_edges.begin(), c_edges.end());
    vertices.push_back(edges.size());
}
pmc::pmc_graph G(vertices, edges);
// [VERIFY: teaser/src/graph.cc:19-32]
```

这是经典 CSR（Compressed Sparse Row）：`vertices[i]` 到 `vertices[i+1]` 是顶点 i 的邻居在 `edges` 中的范围。`vertices` 用 `long long` 是因为 PMC 期望 64-bit 偏移（支持大图）。

`pmc::input` 配置（[VERIFY: teaser/src/graph.cc:35-52]）：

| 字段 | 值 | 含义 |
|------|----|------|
| `algorithm` | 0 | PMC 主算法 |
| `threads` | `params_.num_threads` | 并行度 |
| `lb` | 0 | clique size 下界（搜索后会被填） |
| `ub` | 0 | 上界（被 `max_core + 1` 覆盖） |
| `adj_limit` | 20000 | 顶点数小于该值时用 dense adjacency 形式 |
| `time_limit` | `params_.time_limit` | 求解超时 |
| `heu_strat` | `"kcore"` | 启发式策略：k-core 排序 |
| `vertex_search_order` | `"deg"` | 搜索顺序按度 |

---

## 6. PointCloud 与 Eigen Matrix 的转换

由于公开 API 接受两种输入，需要一致的转换语义：

```cpp
// teaser::PointCloud → Eigen::Matrix3Xd
Eigen::Matrix<double, 3, Eigen::Dynamic> src(3, correspondences.size());
for (size_t i = 0; i < correspondences.size(); ++i) {
    auto src_idx = std::get<0>(correspondences[i]);
    auto dst_idx = std::get<1>(correspondences[i]);
    src.col(i) << src_cloud[src_idx].x, src_cloud[src_idx].y, src_cloud[src_idx].z;
    dst.col(i) << dst_cloud[dst_idx].x, dst_cloud[dst_idx].y, dst_cloud[dst_idx].z;
}
// [VERIFY: teaser/src/registration.cc:563-570]
```

`PointXYZ` 的 `float` 字段被隐式 promoted 到 `double` 进 Eigen 矩阵。

---

## 7. 内存布局与生命周期

### 7.1 主要矩阵的内存占用

设 `N` 为输入对应数，`L = |max_clique_|`：

| 字段 | 维度 | 元素类型 | 字节数 (N=1000, L=200) |
|------|------|----------|------------------------|
| `src_tims_`, `dst_tims_` | 3 × N(N-1)/2 | double | 3 × 499500 × 8 = 11.4 MB 各 |
| `src_tims_map_`, `dst_tims_map_` | 2 × N(N-1)/2 | int | 2 × 499500 × 4 = 3.8 MB 各 |
| `scale_inliers_mask_` | 1 × N(N-1)/2 | bool | 499500 字节 ≈ 0.5 MB |
| `pruned_src_tims_`, `pruned_dst_tims_` | 3 × L | double | 3 × 200 × 8 = 4.7 KB 各 (CHAIN) |
| `max_clique_` | L 个 int | int | 800 字节 |
| `inlier_graph_` | 邻接表 | 多个 `vector<int>` | 取决于 edge 数 |

**N = 10000 时**：
- TIM 矩阵单个 = 3 × ~5×10⁷ × 8 = ~1.2 GB
- src + dst + map 加起来 ~5 GB

这是 TEASER 算法本身的 O(N²) 内存复杂度，与实现无关。N 超过 5000 应优先考虑稀疏化输入或换算法。

### 7.2 生命周期与重入

- `RobustRegistrationSolver` 实例每次 `solve()` 后，内部矩阵保留（getter 仍可访问），下次 `solve()` 会就地覆盖；
- `reset()` 会**重建** solver 子组件并清空若干 vector，但**不**清空矩阵和 mask；
- **不可重入**：多线程不能共享同一实例。如需并行处理多对点云，每线程一个 `RobustRegistrationSolver`。

### 7.3 RAII 与异常安全

- 三个 solver `unique_ptr` 保证析构时释放。
- 没有自定义 `operator new`，全部 STL 容器自动管理。
- 唯一手动对齐 new 是 `RegistrationSolution` 的 `EIGEN_MAKE_ALIGNED_OPERATOR_NEW`——但只在用户 `new RegistrationSolution` 时才有用。库内部全部是栈或成员，不需要 placement new。

---

## 8. 数据结构关系图

```
                  ┌─────────────────────────────────────┐
                  │      RobustRegistrationSolver       │
                  │                                     │
                  │  Params params_                     │
                  │  RegistrationSolution solution_     │
                  │                                     │
                  │  ┌─────────────┐ ┌────────────────┐ │
                  │  │ src_tims_   │ │ src_tims_map_  │ │
                  │  │ dst_tims_   │ │ dst_tims_map_  │ │
                  │  │ 3×N(N-1)/2  │ │ 2×N(N-1)/2     │ │
                  │  └─────────────┘ └────────────────┘ │
                  │  ┌─────────────┐ ┌────────────────┐ │
                  │  │ pruned_*    │ │ *_map_rotation │ │
                  │  │ 3×L_pruned  │ │ 2×L_pruned     │ │
                  │  └─────────────┘ └────────────────┘ │
                  │  ┌─────────────────────────────┐    │
                  │  │ scale_inliers_mask_  (1×M)  │    │
                  │  │ rotation_inliers_mask_(1×Lp)│    │
                  │  │ translation_inliers_mask_   │    │
                  │  └─────────────────────────────┘    │
                  │  ┌─────────────────────────────┐    │
                  │  │ max_clique_     vector<int> │    │
                  │  │ rotation_inliers_           │    │
                  │  │ translation_inliers_        │    │
                  │  └─────────────────────────────┘    │
                  │  ┌─────────────────────────────┐    │
                  │  │ inlier_graph_  : Graph      │────┼──┐
                  │  └─────────────────────────────┘    │  │
                  │  ┌─────────────────────────────┐    │  │
                  │  │ scale_solver_      unique_ptr│   │  │
                  │  │   AbstractScaleSolver        │   │  │
                  │  │     ├─ TLSScaleSolver        │   │  │
                  │  │     └─ ScaleInliersSelector  │   │  │
                  │  │                              │   │  │
                  │  │ rotation_solver_  unique_ptr │   │  │
                  │  │   GNCRotationSolver          │   │  │
                  │  │     ├─ GNCTLSRotationSolver  │   │  │
                  │  │     ├─ FastGlobalRegistration│   │  │
                  │  │     └─ QuatroSolver          │   │  │
                  │  │                              │   │  │
                  │  │ translation_solver_          │   │  │
                  │  │   AbstractTranslationSolver  │   │  │
                  │  │     └─ TLSTranslationSolver  │   │  │
                  │  └──────────────────────────────┘   │  │
                  └─────────────────────────────────────┘  │
                                                           │
        ┌──────────────────────────────────────────────────┘
        │
        ▼
  ┌────────────────────┐         ┌───────────────────────┐
  │ Graph              │  ────▶  │ MaxCliqueSolver       │
  │  adj_list_         │         │   params_             │
  │  num_edges_        │         │   findMaxClique(g)    │
  └────────────────────┘         │     └─ pmc::pmc_graph │
                                 └───────────────────────┘
```

---

## 9. 文档完整性检查表

- [x] 每个公开类型都有 `[VERIFY: file:line]` 标签。
- [x] 字段语义与维度都有显式说明。
- [x] 三种 inlier mask 的维度差异已强调。
- [x] 三个废弃字段（`use_max_clique`、`max_clique_exact_solution`、`solve_exactly`）的兼容路径已分别标注。
- [x] 算法上界与内存量级有具体数字示例。
- [x] PMC interop 的 CSR 转换已说明。
- [x] 类层次（abstract base + 具体实现）与 Strategy 模式关系已图示。

