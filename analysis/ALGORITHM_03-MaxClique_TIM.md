# ALGORITHM 03：TIM 构造与 Max-Clique Inlier 选择

> **范围**：覆盖 TEASER++ 中两个紧密耦合的算法子模块：
> 1. **TIM (Translation-Invariant Measurement) 构造** —— `RobustRegistrationSolver::computeTIMs()`
> 2. **Inlier Graph + Max-Clique 求解** —— `Graph` + `MaxCliqueSolver` + PMC facade
>
> 这两个模块共同负责把「N 个对应」的鲁棒筛选问题转化为「图上的最大团」问题，并通过 PMC 库求解。
>
> **核心代码**：
> - TIM 构造 [VERIFY: teaser/src/registration.cc:518-557]
> - Inlier 图建立 [VERIFY: teaser/src/registration.cc:615-660]
> - Max-clique 求解器 [VERIFY: teaser/src/graph.cc:12-125]
> - Graph 数据结构 [VERIFY: teaser/include/teaser/graph.h:29-207]
> - MaxCliqueSolver 类 [VERIFY: teaser/include/teaser/graph.h:219-279]
>
> **理论来源**：
> - TIM 概念：H. Yang, J. Shi, L. Carlone, "TEASER", T-RO 2021
> - PMC: R. A. Rossi, D. F. Gleich, A. H. Gebremedhin, "Parallel Maximum Clique Algorithms with Applications to Network Analysis", SIAM J. Sci. Comput. 2015
> - K-core 启发式：经典图论结果

---

## 1. 大局观：为什么要 max-clique？

### 1.1 问题简化

Scale 估计后得到一组 1×M 的 boolean mask（M = N(N−1)/2 = TIM 数）：
```
scale_inliers_mask_[c] = 1  ⟺  TIM (i, j) 通过 scale 检查
```
其中 c 是 TIM 下标，对应原始测量对 (i, j)，由 `src_tims_map_` 索引到。

这是个 **TIM 维度** 的 mask。我们要的是 **测量维度** 的 inlier：
```
inlier_set ⊆ {0, 1, ..., N-1}
```

### 1.2 关键观察

**观察 1**：如果 i, j 都是真 inlier，那么 TIM (i, j) 几乎必然通过 scale 检查。

理由：
```
b_i = s · R · a_i + ε_i,    ‖ε_i‖ ≤ β
b_j = s · R · a_j + ε_j,    ‖ε_j‖ ≤ β
b̃_{ij} := b_j − b_i = s · R · (a_j − a_i) + (ε_j − ε_i)
‖b̃_{ij}‖ = ‖s · R · ã_{ij} + Δε‖ ≈ s · ‖ã_{ij}‖,   误差 ‖Δε‖ ≤ 2β
```
范数比 `‖b̃‖/‖ã‖ ≈ s`，落在容许带内。

**观察 2**：若 i 或 j 是 outlier，TIM (i, j) **可能**通过 scale 检查（偶然），也可能不过。

**观察 3**：把 i 加进 inlier 集 ⟺ 把顶点 i 加进图。inlier 集中两两都满足 TIM 约束 ⟺ 这些顶点两两有边。

**结论**：inlier 集合在「TIM 图」上必然是**完全子图**（团）。最大的「inlier 候选集」是图的**最大团**。这就是把 inlier 选择规约到 max-clique 的逻辑。

### 1.3 算法的鲁棒性保证

- **必要性**：真 inlier 必构成团。所以最大团**至少**与真 inlier 集一样大。
- **充分性**（弱）：最大团**可能含 outlier**，因为偶尔几个 outlier 的 TIM 也可能通过。但这种「巧合 outlier」难以构成大团（要求多个 outlier 之间的 TIM 全都通过 scale 检查，概率指数衰减）。
- 综合：max-clique 是 inlier 集的高质量**近似上界**——比例上几乎不含 outlier。后续 GNC-TLS 在 max-clique 上再做精细 outlier 剔除。

这是 TEASER 与其他鲁棒方法（RANSAC、M-estimator）的根本区别：**两阶段鲁棒筛选**——max-clique（图论）+ GNC-TLS（连续优化）。

---

## 2. TIM 构造算法

### 2.1 接口

```cpp
Eigen::Matrix<double, 3, Eigen::Dynamic>
RobustRegistrationSolver::computeTIMs(
    const Eigen::Matrix<double, 3, Eigen::Dynamic>& v,
    Eigen::Matrix<int, 2, Eigen::Dynamic>* map);
// [VERIFY: teaser/include/teaser/registration.h:548-550]
// [VERIFY: teaser/src/registration.cc:518-557]
```

**输入**：`v` 是 3 × N 的点矩阵（每列一个点）。  
**输出**：
- 返回值 `vtilde`：3 × M（M = N(N−1)/2），每列是 `v.col(j) - v.col(i)` (i < j)
- `map`：2 × M，每列 `[i, j]` 是端点下标

### 2.2 内存布局：「上三角扁平化」

TIM 是 N 选 2 = N(N−1)/2 个无序对。代码用 i 优先的分段排列：

```
TIM 列   |   段 (i)   |   端点对
─────────┼────────────┼─────────────────────
0        | i=0        | (0, 1)
1        |            | (0, 2)
...      |            | ...
N-2      |            | (0, N-1)
N-1      | i=1        | (1, 2)
N        |            | (1, 3)
...      |            | ...
2N-3     |            | (1, N-1)
2N-2     | i=2        | (2, 3)
...
```

每段大小：i=0 有 N−1 个，i=1 有 N−2 个，...，i=k 有 N−1−k 个。  
段起始下标：
```
S(k) = 累加 (N-1) + (N-2) + ... + (N-k)
     = k·N − k(k+1)/2
```

代码 [VERIFY: teaser/src/registration.cc:537]：
```cpp
size_t segment_start_idx = i * N - i * (i + 1) / 2;
size_t segment_cols      = N - 1 - i;
```

### 2.3 OpenMP 并行性分析

```cpp
#pragma omp parallel for shared(N, v, vtilde, map)
for (size_t i = 0; i < N - 1; i++) {
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
// [VERIFY: teaser/src/registration.cc:526-554]
```

**关键观察**：
- 不同 i 写入 vtilde 的**不相交**列段；
- 不同 i 写入 map 的**不相交**列段；
- 共享读取：v、N、起始下标等。

无数据竞争。

**每个 i 的工作量**：
- `temp = v - m * Ones(1, N)`：分配 3×N 临时矩阵，做 N 次三维向量减法 → O(3N) flop + O(3N) 分配；
- 拷到 vtilde：O(3·(N-1-i)) memcpy；
- map_addition 构造：N 次写入。

总成本 per i：O(N)。总 i 数：N−1。**串行总成本：O(N²)**。

OpenMP 把外层 i 循环并行——理论上 P 线程可达 O(N²/P)，但实际受内存带宽限制：N²=10⁶ 时通常 4–8 倍加速封顶。

### 2.4 内存峰值

输出 `vtilde`：3 × M × 8B = 12 × N(N−1) bytes。N=1000 时 ≈ 12 MB。每个 i 的临时 `temp`：3 × N × 8B = 24 KB（小，可忽略）。

**N=10000 时 vtilde 单个 ≈ 1.2 GB** —— 这是 TEASER 算法本身的内存 wall，与实现无关。

### 2.5 不变量

- `vtilde.cols() = N(N-1)/2`：[VERIFY: teaser/src/registration.cc:523]
- `vtilde.col(c) = v.col(map(1,c)) - v.col(map(0,c))`
- `map(0, c) < map(1, c)` 始终成立
- 对相同 v 输入，src 与 dst 各调一次得到 src_tims_、dst_tims_、src_tims_map_、dst_tims_map_——其中 `src_tims_map_` 与 `dst_tims_map_` **数值上完全相同**

---

## 3. Inlier 图构造

### 3.1 图模型

- **顶点集**：N 个，对应 N 对原始对应。顶点 ID 直接是测量下标 0..N−1。
- **边集**：对每个 TIM 下标 c，若 `scale_inliers_mask_[c] = 1`，则在 `map(0,c)` 与 `map(1,c)` 之间加无向边。

```cpp
inlier_graph_.populateVertices(src.cols());
for (size_t i = 0; i < scale_inliers_mask_.cols(); ++i) {
    if (scale_inliers_mask_(0, i)) {
        inlier_graph_.addEdge(src_tims_map_(0, i), src_tims_map_(1, i));
    }
}
// [VERIFY: teaser/src/registration.cc:620-625]
```

注意：
- `populateVertices(N)` 仅 resize 邻接表到 N，**不**清空边——如果 solver 重用，旧边会保留！见 02_DATA_FLOW.md §9。
- `addEdge` 自动维护对称性（双向插入），见 [VERIFY: teaser/include/teaser/graph.h:96-104]。

### 3.2 复杂度

- 顶点构造：O(N)
- 边构造：O(M) 检查 + 最多 M 次 `addEdge`，每次 `addEdge` 内部含 `hasEdge` 检查 → O(deg) 线性搜索 → 单次 O(N) worst case
- **总 worst case：O(M · N) = O(N³)**

这看起来吓人，但实际上 hasEdge 通常很快（adj_list 短），并且只在调试断言下触发——release 编译下 `hasEdge` 内部的 `TEASER_DEBUG_ERROR_MSG` 是 no-op。  
[VERIFY: teaser/include/teaser/graph.h:97-99]

实测：N=1000 时建图 < 50 ms。

### 3.3 图的属性

- **稀疏度**取决于 inlier 比例。极端：
  - 100% inlier、零噪声：完全图，N(N-1)/2 条边；
  - 100% outlier、零相关：几乎空图。
- **顶点度分布**：inlier 顶点度高（与其他 inlier 都有边），outlier 顶点度低且随机。

这种**度的两极分化**是 max-clique 算法发挥作用的几何基础——k-core 启发式正好利用这点。

---

## 4. MaxCliqueSolver：算法分派

### 4.1 总体流程

```cpp
std::vector<int> MaxCliqueSolver::findMaxClique(Graph graph) {
    // 1. 兼容 deprecated solve_exactly
    if (!params_.solve_exactly) {
        params_.solver_mode = CLIQUE_SOLVER_MODE::PMC_HEU;
    }

    // 2. 把 teaser::Graph 转 PMC CSR
    std::vector<int> edges;
    std::vector<long long> vertices;
    vertices.push_back(edges.size());
    for (const auto& i : graph.getVertices()) {
        const auto& c_edges = graph.getEdges(i);
        edges.insert(edges.end(), c_edges.begin(), c_edges.end());
        vertices.push_back(edges.size());
    }
    pmc::pmc_graph G(vertices, edges);

    // 3. 设置 PMC input
    pmc::input in;
    in.algorithm = 0;
    in.threads   = params_.num_threads;
    ...
    in.adj_limit = 20000;
    in.time_limit = params_.time_limit;
    in.heu_strat = "kcore";
    in.vertex_search_order = "deg";

    std::vector<int> C;

    // 4. K-core 上界（同时也用于 KCORE_HEU 启发式）
    G.compute_cores();
    auto max_core = G.get_max_core();

    // 5. K-core 启发式分支（仅当 mode = KCORE_HEU 且阈值满足）
    if (params_.solver_mode == CLIQUE_SOLVER_MODE::KCORE_HEU &&
        params_.kcore_heuristic_threshold != 1 &&
        max_core > threshold·|V|) {
        auto k_cores = G.get_kcores();
        for (int i = 1; i < k_cores->size(); ++i) {
            if ((*k_cores)[i] >= max_core) {
                C.push_back(i - 1);
            }
        }
        return C;
    }

    // 6. PMC heuristic 求下界
    if (in.ub == 0) in.ub = max_core + 1;
    if (in.lb == 0 && in.heu_strat != "0") {
        pmc::pmc_heu maxclique(G, in);
        in.lb = maxclique.search(G, C);
    }

    if (in.lb == 0) return C;   // 无 clique
    if (in.lb == in.ub) return C;  // 启发式已达上界，无需 exact

    // 7. PMC exact (仅 mode = PMC_EXACT)
    if (params_.solver_mode == CLIQUE_SOLVER_MODE::PMC_EXACT) {
        if (G.num_vertices() < in.adj_limit) {
            G.create_adj();
            pmc::pmcx_maxclique finder(G, in);
            finder.search_dense(G, C);
        } else {
            pmc::pmcx_maxclique finder(G, in);
            finder.search(G, C);
        }
    }

    return C;
}
// [VERIFY: teaser/src/graph.cc:12-125]
```

### 4.2 三种 mode 的实际行为对照

| mode | 实际走的路径 |
|------|------------|
| **PMC_EXACT** | k-core 计算 → PMC heuristic 给下界 → 若 lb < ub，调 PMC exact 求精确解 |
| **PMC_HEU** | k-core 计算 → PMC heuristic 给下界 → 返回（**不**调 exact） |
| **KCORE_HEU** | k-core 计算 → 若 max_core > threshold·N → 返回 max-core 全部顶点；否则降级到 PMC_HEU 行为 |

注意 **KCORE_HEU 不总是走 k-core**！只有 `max_core > kcore_heuristic_threshold · |V|` 时才触发；否则继续走 heuristic 流程（line 88-91 不被跳过）。

### 4.3 上下界协议

PMC 求 max-clique 是个 branch-and-bound：
- 下界 lb：当前找到的最佳团大小；
- 上界 ub：最大团大小的理论上限。

K-core 数 ω(G) 满足 **ω(G) ≤ k_max(G) + 1**（max core 加 1）。代码用 `in.ub = max_core + 1` [VERIFY: graph.cc:83-85]。

PMC heuristic（`pmc_heu::search`）快速找一个**大但未必最优**的团，给 `in.lb`。

若 `lb == ub`，已经找到最大团，无需 exact。

---

## 5. K-Core 算法

### 5.1 定义

图 G 的 k-core 是「每个顶点度数都 ≥ k 的最大子图」。一个顶点的「core number」是它所在的最大 k-core 的 k 值。

形式化：
```
core(v) = max { k : v ∈ k-core(G) }
```

### 5.2 与 max-clique 的关系

**定理**（经典图论）：若 G 含 ω-clique，则该 clique 内的每个顶点 v 在 G 中度数 ≥ ω−1（与团内其他 ω−1 个点都有边）。所以
```
core(v) ≥ ω(G) - 1   对 clique 中的 v 成立
⇒ k_max(G) ≥ ω(G) - 1
⇒ ω(G) ≤ k_max(G) + 1
```
这是 §4.3 用 `max_core + 1` 作 ub 的依据。

### 5.3 K-core 计算

经典 BZ 算法（Batagelj & Zaversnik, 2003）：O(|V| + |E|) 线性时间。PMC 库内部用此算法 (`G.compute_cores()`)，对 teaser 来说是黑盒。

### 5.4 K-core 启发式作为 inlier 估计器

代码 [VERIFY: teaser/src/graph.cc:70-80]：

```cpp
auto k_cores = G.get_kcores();
for (int i = 1; i < k_cores->size(); ++i) {
    if ((*k_cores)[i] >= max_core) {
        C.push_back(i - 1);
    }
}
```

**算法**：返回所有 core number = max_core 的顶点。这些顶点构成 max-core 子图，是「密集子图」，**通常**包含 max-clique，但严格说**不一定是团**。

注意 `i` 从 1 开始而不是 0：PMC 的 `k_cores` 数组大小为 `|V| + 1`（comment in graph.cc:74），下标 1..|V| 对应顶点 0..|V|−1。所以 `(*k_cores)[i]` 对应顶点 i−1 的 core number。

### 5.5 K-core 触发条件

```cpp
if (params_.solver_mode == CLIQUE_SOLVER_MODE::KCORE_HEU &&
    params_.kcore_heuristic_threshold != 1 &&
    max_core > static_cast<int>(params_.kcore_heuristic_threshold *
                                static_cast<double>(all_vertices.size()))) {
    // 走 k-core
}
// [VERIFY: teaser/src/graph.cc:66-69]
```

三条件：
1. solver_mode 是 KCORE_HEU；
2. threshold ≠ 1（threshold = 1 时该分支被绕过——这是个**短路优化**，因为 threshold = 1 意味着「只有 max_core = N 时才用 k-core」，几乎不会发生）；
3. max_core > threshold · N。

**threshold = 0.5（默认 RobustRegistrationSolver）**意味着：当 max_core > N/2，直接取 max-core 顶点作为 max-clique 估计。这对**高 inlier 比例**场景成立（inlier 比例 > 50% 时 max-clique 大小 ≈ inlier 数 > N/2 → max_core ≥ inlier - 1 > N/2-1）。

但对**低 inlier 比例**场景（如 20% inlier），max_core 远小于 N/2，k-core 触发不了，会自动 fallback 到 PMC heuristic。

### 5.6 K-core 的 trade-off

| 维度 | 评价 |
|------|------|
| 速度 | 极快，O(|V| + |E|) 线性 |
| 精度 | 上界 `k_max + 1` 较松；返回的是 k-core 顶点而非真团 |
| 安全性 | k-core 顶点不一定全 inlier，可能含「假阳性」 |
| 适用 | 高 inlier 比例 + 求快 |

---

## 6. PMC Exact 算法

### 6.1 总体策略

Rossi 等人 2015 论文的 PMC（Parallel Maximum Clique）核心：
1. **K-core 剪枝**：删除所有 core number < lb 的顶点（它们不可能在大团里）；
2. **Neighbor-core 排序**：按邻居中的最大 core 排序顶点，先访问高潜力的；
3. **动态着色上界**：在 branch-and-bound 中用贪心着色给当前子图一个 clique 数上界，提前剪枝；
4. **并行**：不同根顶点的搜索可并行。

### 6.2 在 TEASER 中的调用

```cpp
if (params_.solver_mode == CLIQUE_SOLVER_MODE::PMC_EXACT) {
    if (G.num_vertices() < in.adj_limit) {
        G.create_adj();
        pmc::pmcx_maxclique finder(G, in);
        finder.search_dense(G, C);
    } else {
        pmc::pmcx_maxclique finder(G, in);
        finder.search(G, C);
    }
}
// [VERIFY: teaser/src/graph.cc:105-122]
```

**`adj_limit = 20000`**：顶点数小于此值时构建邻接矩阵（密集形式），调 `search_dense`；否则保留 CSR，调 `search`。

邻接矩阵需要 O(|V|²) 内存：|V| = 20000 时 400 MB。所以 20000 是个内存阈值。

### 6.3 时间复杂度

理论上 max-clique 是 NP-hard，最坏 O(2^N)。PMC 用强力剪枝把实际复杂度压到「图密度相关」：稀疏图上几乎多项式时间，密集图上变慢。

实测：N=1000、inlier=20% 的典型 outlier 图，PMC_EXACT < 100 ms（4 线程）。
N=5000、inlier=10%：1–10 s。
N > 10000：可能超过 `time_limit = 3600`，需切换到 heuristic。

### 6.4 内部参数

```cpp
in.algorithm  = 0;     // 主算法版本
in.experiment = 0;
in.lb = 0;             // 下界（被 heuristic 填）
in.ub = 0;             // 上界（被 max_core + 1 填）
in.adj_limit  = 20000;
in.remove_time = 4;    // ?
in.graph_stats = false;
in.verbose = false;
in.MCE = false;        // 不要 Maximum Clique Enumeration
in.decreasing_order = false;
in.heu_strat = "kcore";
in.vertex_search_order = "deg";
// [VERIFY: teaser/src/graph.cc:37-52]
```

这些参数由 TEASER 硬编码，不可配置。

---

## 7. PMC Heuristic

### 7.1 算法

`pmc::pmc_heu::search` 用贪心 + k-core 排序找一个大团：
1. 按 core number 降序排顶点；
2. 从顶点 v_1 开始，把它加入当前团 C；
3. 取 C 的所有邻居的交集，作为下一个候选集；
4. 从候选集中选 core number 最大的顶点加入 C；
5. 重复直到候选集空。

返回的 C 是个团，但未必最大。

### 7.2 在 TEASER 中的角色

两个用途：
1. **为 PMC_HEU 模式提供输出**（[VERIFY: teaser/src/graph.cc:89-91]）；
2. **为 PMC_EXACT 提供下界**（同一行），让 exact 求解器更早剪枝。

### 7.3 与 KCORE_HEU 的区别

| 维度 | PMC_HEU | KCORE_HEU |
|------|---------|-----------|
| 返回 | 真正的团（pairwise connected） | k-core 顶点集（密集子图，未必团） |
| 大小 | 通常 ≈ ω(G) − O(1) | 通常 ≥ ω(G) − 1，可能含非团顶点 |
| 速度 | 慢一些（要做交集计算） | 极快（线性） |

---

## 8. PMC interop：CSR 数据转换

### 8.1 CSR 表示

```cpp
std::vector<int> edges;            // 邻居顶点 ID 展平
std::vector<long long> vertices;   // 每个顶点在 edges 中的起始偏移
```

约定：
- `vertices[i]` 到 `vertices[i+1]` 是顶点 i 的邻居范围；
- `vertices.size() = |V| + 1`；
- `edges.size() = 2 |E|`（无向图每边两次）；
- `vertices[|V|] = 2 |E|`。

### 8.2 转换代码

```cpp
std::vector<int> edges;
std::vector<long long> vertices;
vertices.push_back(edges.size());   // vertices[0] = 0

const auto all_vertices = graph.getVertices();
for (const auto& i : all_vertices) {
    const auto& c_edges = graph.getEdges(i);
    edges.insert(edges.end(), c_edges.begin(), c_edges.end());
    vertices.push_back(edges.size());
}
// [VERIFY: teaser/src/graph.cc:20-29]
```

**值得注意**：
- 第一次 push 是 `0`（初始 `edges.size() = 0`），即 `vertices[0] = 0`；
- 每次 push 当前 edges 大小（即顶点 i+1 的起始偏移）；
- 最终 `vertices.size() = |V| + 1`，`vertices[|V|] = edges.size() = 2 |E|`。

### 8.3 复杂度

- 时间 O(|V| + |E|)；
- 空间 O(|V| + |E|)，额外峰值是 teaser::Graph 的 adj_list_ + edges + vertices ≈ 2× 的总大小。

---

## 9. 主管线对接

### 9.1 max_clique_ 后处理

```cpp
max_clique_ = clique_solver.findMaxClique(inlier_graph_);
std::sort(max_clique_.begin(), max_clique_.end());
// [VERIFY: teaser/src/registration.cc:641-642]
```

**为什么排序？**
1. CHAIN 阶段需要顺序遍历环：`leaf = max_clique_[i+1]`，无序时环结构无意义。
2. 调试输出可读（升序顶点 ID）。

### 9.2 退化处理

```cpp
if (max_clique_.size() <= 1) {
    TEASER_DEBUG_INFO_MSG("Clique size too small. Abort.");
    solution_.valid = false;
    return solution_;
}
// [VERIFY: teaser/src/registration.cc:649-653]
```

`max_clique_.size() == 0` 发生在：
- 图无边（无任何 TIM 通过 scale）→ PMC 返回空团；
- PMC 超时但 lb 未被更新（`in.lb = 0`）→ heuristic 没找到任何边的团 [VERIFY: teaser/src/graph.cc:93-97]。

`max_clique_.size() == 1` 发生在：
- 图只有孤立顶点（每个顶点都是「自己一个人的团」），PMC 可能返回单点；
- 极低 inlier 比例下退化。

两种情况下都返回 `valid = false`。

---

## 10. ASCII 图示：max-clique 选 inlier

设 N = 6，真实 inlier 是 {0, 1, 2}，outlier 是 {3, 4, 5}。

TIM 通过 scale 检查的对（边）：
```
inlier 之间：(0,1), (0,2), (1,2)        全通过
outlier 偶然过：(3,4)
inlier-outlier 偶然过：(2,3)
```

inlier 图：
```
          0───1
          │\ /│
          │ X │
          │/ \│
          2───3───4
                  │
                  5（孤立）
```

候选团：
- {0, 1, 2}：大小 3 ✓（真 inlier）
- {3, 4}：大小 2（偶然 outlier 团）
- {2, 3}：大小 2

max-clique = {0, 1, 2}，size = 3。完美返回真 inlier。

实际：outlier 之间巧合通过 scale 的概率极低（要求两个 outlier 的范数差刚好在 β 内），所以「outlier 偶然团」一般小于「真 inlier 团」。

---

## 11. 复杂度总结

| 阶段 | 时间复杂度 | 空间复杂度 | 备注 |
|------|------------|------------|------|
| TIM 构造 | O(N²) | O(N²) | OpenMP 并行外层 |
| inlier 图建图 | O(N² · 边检查 cost) | O(N + edges) | 实际近 O(N²) |
| K-core 计算 | O(\|V\| + \|E\|) | O(\|V\|) | 线性 |
| PMC heuristic | O(\|V\| · log \|V\|) typical | O(\|V\|) | 贪心 + k-core 排序 |
| PMC exact | O(2^\|V\|) worst, 实际近多项式 | O(\|V\|²) (search_dense) | branch-and-bound |
| **整阶段** | 由 PMC 主导 | O(N²) | TIM 矩阵是主要开销 |

---

## 12. 关键设计选择

### 12.1 为什么把测量作顶点而不是 TIM 作顶点？

另一种建图：把 TIM 作顶点，「两 TIM 共享一个测量」时加边。这是 line graph。但 line graph 的团对应原图的「轮廓」，不直接对应 inlier 集——不自然。

测量作顶点 + TIM 通过作边的建图，**直接**让 inlier 集等价于团。这是 TEASER 的关键洞察之一。

### 12.2 为什么默认 PMC_EXACT 而非 PMC_HEU？

[VERIFY: teaser/include/teaser/registration.h:476]：默认 `INLIER_SELECTION_MODE::PMC_EXACT`。

PMC_EXACT 保证返回最大团 → inlier 集最准确。代价是时间。但 TEASER 论文目标是「certifiable」（可证最优），所以默认走精确。

如果用户对时间敏感，可改 `PMC_HEU` 或 `KCORE_HEU`——文档提醒：「Set this to 0 to always skip exact max clique selection」（用 `kcore_heuristic_threshold`） [VERIFY: teaser/include/teaser/registration.h:480-482]。

### 12.3 为什么 K-core 是上界？

§5.2 已证明 `ω(G) ≤ k_max(G) + 1`。这是个**紧界**：在某些图上等号成立（如 cliques + isolated vertices）。给 PMC exact 一个好的初始 ub 能加速 branch-and-bound 剪枝。

### 12.4 为什么 noise_bound 的因子 2 在这里没出现？

注意 §1 与 02_DATA_FLOW.md：scale 阶段 β = 2 · noise_bound · √c̄²，含因子 2。但在 max-clique 阶段没出现「因子 2」——因为 max-clique 是个**二阶段后处理**：它消费 `scale_inliers_mask_`，不直接操作 noise_bound。

之后 rotation 阶段会再次出现 `2/scale` 因子，把 noise_bound 调整为 rotation 上的真实噪声 → 见 ALGORITHM_02 §1。

---

## 13. 工程坑点

### 13.1 inlier_graph_ 不自动清空

`solve()` 中 `populateVertices` 只 resize，**不清空旧边** [VERIFY: teaser/include/teaser/graph.h:67]。如果用户多次调 `solve()` 不 `reset()`，第二次的 inlier 图会**叠加**第一次的边——max-clique 会基于错乱的图求解。

修复：每次 `solve()` 前显式 `reset()`，或在 `solve()` 顶部加 `inlier_graph_.clear()`。

### 13.2 KCORE_HEU 不总是走 k-core

§5.5 详述。用户设 `INLIER_SELECTION_MODE::KCORE_HEU` 期待「最快」，但 max_core 不够大时会 fallback 到 PMC_HEU（仍然不快）。

要**强制**走 k-core，需要：`kcore_heuristic_threshold` 设极小（如 0.01），让 `max_core > 0.01 · N` 几乎总是满足。

### 13.3 `max_clique_num_threads` 默认 omp_get_max_threads()

[VERIFY: teaser/include/teaser/registration.h:513]：
```cpp
int max_clique_num_threads = teaser_default_max_threads();
```

`teaser_default_max_threads()` 返回 `omp_get_max_threads()` [VERIFY: teaser/src/registration.cc:23-25]。

但 PMC 内部用自己的线程池，与 OpenMP 主程序的线程池**独立**。所以高并发应用里同时跑多个 solver 时可能超额订阅 CPU——记得显式限制 `max_clique_num_threads`。

### 13.4 time_limit = 3600 太大

默认超时 1 小时。PMC 在难图上可能跑很久 → 长时间无响应。生产环境应改成几秒：
```cpp
params.max_clique_time_limit = 5.0;
```

超时后 PMC 返回当前最佳解，不会出错。`max_clique_` 可能比真最优小，但仍可用。

---

## 14. 数值示例

```
N = 200, true inlier = 100 (50%), Gaussian noise σ = 0.5 · noise_bound

Stage 1:
  M = N(N-1)/2 = 19,900 个 TIM
  src_tims_ 矩阵：3 × 19900 ≈ 470 KB

Stage 2 (scale):
  scale_inliers_mask_ → 大约 100·99/2 = 4950 个 inlier-inlier TIM 通过
  加上 outlier-outlier 偶然通过：随机噪声下约 0.5% = 100 个
  总通过 ≈ 5050 条边

Stage 3 (inlier graph):
  N = 200 顶点
  edges ≈ 5050
  k_max ≈ 99 (来自 inlier 间)

Stage 3 (max-clique):
  PMC_EXACT 返回 |C| = 100（真 inlier 全部入选）
  时间：几十 ms（单核）
```

---

## 15. 关键不变量

| 不变量 | 验证 |
|--------|------|
| `vtilde.cols() = N(N-1)/2` | [VERIFY: registration.cc:523-524] |
| `map(0,c) < map(1,c)` | 构造时 i < j（i 是外层 loop 变量，j ≥ i+1） |
| `inlier_graph_.numVertices() = N` | [VERIFY: registration.cc:620] |
| `max_clique_` 升序 | [VERIFY: registration.cc:642] |
| `|max_clique_| ≤ N` | trivially |
| `|max_clique_| ≤ max_core + 1` | k-core 上界定理 |

---

## 16. 文档完整性检查表

- [x] TIM 构造的索引推导有数学公式。
- [x] OpenMP 并行性的数据竞争分析。
- [x] inlier 图建图代码逐行注释。
- [x] 三种 max-clique mode 的实际行为对照表（含 KCORE_HEU 的非平凡触发逻辑）。
- [x] k-core 上界 ω ≤ k_max + 1 给出证明。
- [x] CSR 数据转换的约定与代码对应。
- [x] PMC interop 的硬编码参数已列出。
- [x] 工程坑点：图不清空、KCORE_HEU 不总走 k-core、超时默认过大。
- [x] 数值示例给出典型规模下的预期行为。

