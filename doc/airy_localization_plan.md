# Airy 雷达在已知地图中的定位技术方案

> **本版本状态**（截至 2026-05-11）：
> - ✅ Phase 1（需求 2 / 局部精化）：已实现 `examples/local_refinement/`
> - ✅ Phase 2（需求 1 / 全局定位，SC 路径）：已实现 `examples/global_localization_sc/`，实测 19/22 = 86%
> - ✅ Phase 3（STD 召回率评估）：完成，**STD 在该场景 0/22 = 0%，不进 Phase 4**。详见 `analysis/experiments/std_airy_22submaps.md`
> - ❌ Phase 4（STD 集成）：**取消**，路线确定为 SC + TEASER + submap 多帧累积
> - ✅ Phase 5（BBS 备选路径评估）：已实现 `examples/global_localization_bbs/`，实测 20/22 = 91%（中等阈值，rot<5° AND trans<1m）、pose 返回率 22/22 = 100%。详见 `analysis/experiments/bbs_airy_22submaps.md`。**决策：作为 SC 主线的兜底，进 hybrid pipeline**（理由：rot 精度 1.87° vs SC 的 0.5° 偏弱，不适合做默认主线）
> - ✅ Phase 6（ISC 评估）：已实施，新 session 23-submap 上 ISC vs SC 中等阈值 17/23 = 17/23 完全打平。详见 `analysis/experiments/isc_airy_23submaps.md`。**决策：归档**（理由：聚合零净收益；ensemble 延迟翻倍换 +2 case 不划算；保留 ISC 代码 + mapping intensity 通路作未来候选基础设施）
>
> 修订历史：在 v1 基础上根据 TEASER++ 实际代码分析（`analysis/` 目录下 6 篇文档）做了若干技术更正，最显著的是 Quatro 的 4DOF 误解、N² 内存墙、solver 生命周期等。Phase 5 在 v2 基础上新增 BBS 备选评估。

---

## 1. 需求定义

### 需求 1：全局定位（Global Localization）

- **输入**：完整地图（PLY）+ 单帧 Airy 扫描
- **约束**：完全无先验位姿
- **目标**：6DOF 位姿，目标精度 ±0.1–0.2 m

### 需求 2：局部优化（Local Refinement）

- **输入**：完整地图 + 单帧扫描 + 粗位姿
- **约束**：粗位姿误差 3-10 m / 几度
- **目标**：精化到 ±0.05–0.1 m / ±0.5°

---

## 2. 实施现状

### 2.1 需求 2 → `examples/local_refinement/` ✅

**架构**：FPFH 特征匹配 + TEASER（Quatro 或 GNC-TLS）鲁棒配准。

**实测**（demo 模式 / bunny 数据）：92% inlier、3.7 mm RMSE、与真值差 ~1° / ~2 mm、耗时 130 ms。

**关键工程要点（已在代码里落实）**：

1. **每次 solve 新建 solver 实例**——避开 `RobustRegistrationSolver` 内部 `inlier_graph_` 累加和 Quatro 的 `static double noise_bound_sq` 缓存隐患。
2. **`max_correspondences` 硬截断（默认 3000）**——TEASER TIM 矩阵是 `3 × N(N-1)/2 × 8` byte 的 double，N=3000 时单矩阵 ~100 MB，N=5000 时 ~300 MB。超过 3000 一律按 FPFH 距离 top-N 截断。
3. **`max_clique_time_limit = 10 s`**——TEASER 默认 3600 秒，难场景下会卡死，必须改小。
4. **`PointCloud::empty()` 不是 const**——上游 bug；代码里改用 `size() == 0` 绕过。
5. **`PLYReader::read` 返回 `int`**（0=success），非枚举。
6. **质量指标三件套**：correspondences 数、inlier 数+比例、inlier RMSE 同时输出，建议阈值见 README。

### 2.2 需求 1 → `examples/global_localization_sc/` ✅

**架构**：Scan Context 离线建库 + 在线检索 top-K + 对每个候选跑一次 FPFH+TEASER 验证 + 按 inlier 数选胜。

**模块组织**：
- `third_party/scancontext/`：vendor 自 `gisbi-kim/scancontext_tro`，仅去掉了 .h 中未使用的 ROS / OpenCV / cv_bridge / pcl_conversions includes（算法主体未改）。
- 主程序三个子命令：`--mode build | query | localize`。
- DB 文件是自定义二进制（magic `SCDB` + 描述子 + anchor 位置），不依赖外部序列化库。

**SC 关键设计选择**：
- **不**走 `SCManager::detectLoopClosureID()` —— 那里硬编码了 `NUM_EXCLUDE_RECENT = 50` 与 `SC_DIST_THRES = 0.13`，是 loop-closure 场景设计，不适合 cold-start 全局定位。
- 直接用 SC 的三个 public 原语 `makeScancontext` / `makeRingkeyFromScancontext` / `distanceBtnScanContext`，自己实现 top-K 检索（暴力对比所有 anchor，~1000 anchor 内毫秒级足够）。
- anchor 处取局部云时**平移到原点**——SC 假设 scan 在 scanner 局部坐标系。

**SC 参数边界**：`SCManager::PC_MAX_RADIUS = 80m` 是为车载 LiDAR 调的硬编码常量。煤矿巷道（巷道断面 4-6m、单帧有效半径 30-80m、巷道沿向 100+ m）在 SC 舒适区；室内小场景（< 10m）不适用，要改 .h 文件常量。

---

## 3. 关键技术更正（相对 v1）

### 3.1 ❌ Quatro 不估 4DOF

v1 写「Quatro 估 4DOF (x, y, z, yaw)」。**错**。Quatro 内部（`registration.cc:286-414`）只对 src/dst 的 `topRows(2)` 跑 SO(2) GNC-TLS，输出 R 的左上 2×2 块为 yaw 旋转，其余位置保持 identity。

- 平移 (x, y, z) 是后续 `TLSTranslationSolver` 给的，与 Quatro 无关；
- Quatro **不接受** IMU 的 roll/pitch 作输入。

**正确用法**（融合 IMU 时）：
```
1. 用 IMU 算 R_rp = Ry(pitch) · Rx(roll)
2. 构造 rough_pose: T_rough = [R_yaw_guess · R_rp | t_guess]
3. 程序内部 scan_init = T_rough · scan（roll/pitch 提前对齐到地图坐标系）
4. Quatro 只估 yaw 微调 + TLS 估 t 微调
5. T_refined = T_correction · T_rough 自动把 R_rp 复合上
```

### 3.2 ❌ "煤矿场景 roll/pitch ≈ 0" 过强

煤矿巷道常有坡度（运输大巷、盘区下山）。**默认走 GNC_TLS（完整 SO(3)）**，仅在 IMU 测得 |roll|, |pitch| < 2° 时降级到 Quatro。

### 3.3 ❌ 未考虑 N² 内存墙

TEASER TIM 矩阵 `3 × N(N-1)/2 × 8` byte：

| N | 单 TIM 矩阵 | 整体 |
|---|------|------|
| 1000 | 12 MB | ~60 MB ✅ |
| 3000 | 108 MB | ~300 MB ⚠️ |
| 5000 | 300 MB | ~1 GB ❌ |
| 10000 | 1.2 GB | ~4 GB ❌ |

**对策**：FPFH 匹配后**强制截断到 N ≤ 3000**（按描述子距离 top-N），且 `max_clique_time_limit` 改小到 5-10 s。

### 3.4 ❌ `RegistrationSolution.valid` 不完整

`valid` **只**在 max-clique ≤ 1 时被显式置 false（`registration.cc:649-653`）。GNC 不收敛、scale 退化等情况都不会标记。

**对策**：上层加额外质量检查（inlier 数、RMSE、`getGNCRotationCostAtTermination()`）。

### 3.5 ❌ 多次 `solve()` 不 `reset()` 会泄漏状态

- `inlier_graph_`：`populateVertices` 只 resize 不清空 edges，第二次 `solve()` 边会**累加**；
- `rotation_inliers_`：每次 `emplace_back` 不清空；
- Quatro 的 `static double noise_bound_sq`：函数级 static，只初始化一次。

**对策**：每次 `solve` 新建 solver 实例（Quatro 模式必须，其他模式建议）。

### 3.6 ❌ inlier 比例阈值不够

v1 用 "inlier 比例 > 30%" 作判据。问题：稀疏 Airy + cap=3000 时可能只有 50-200 个对应，30% 也才几十个绝对 inlier，未必够稳定估计 SO(3)。

**对策**：同时检查 **inlier 比例 ≥ 30% 且 inlier 绝对数 ≥ 10**。

---

## 4. Phase 3-4：STD 路径（已结案，**Phase 4 取消**）

### 4.1 实施情况

- **源码 vendor + ROS 剥离已完成**：`examples/global_localization_std/third_party/std/`
- **构建依赖**：Ceres 2.1+（系统 apt 装的 2.0 不行，已装 2.1+ 到 `/usr/local`）
- **评估工具**：`examples/global_localization_std/std_eval_22submaps.cc`，做 22-submap leave-one-out

### 4.2 Phase 3 结果（2026-05-11）—— 详见 `analysis/experiments/std_airy_22submaps.md`

| 方法 | 准确率 | 平移误差中位数 |
|------|--------|---------------|
| SC + TEASER（Phase 2） | **19/22 = 86.4%** | 0.06 m |
| STD top-1（Phase 3） | **0/22 = 0%** | 23 m |

STD 在我们的 Airy submap 数据上**完全失败**——平移误差 4-44m，多个 case rotation error ≈ 180°（场景对称导致 yaw flip）。

**失败机理**（机理分析见实验报告）：
1. 稀疏点云的 plane 检测退化（默认 `voxel_init_num=10` 在 2-8k 点 submap 上得 0 描述子，要降到 5）；
2. 走廊几何重复性 → 三角形哈希查表 false positive 极多；
3. 几何验证投票在对称场景下被 180°-flipped 解迷惑。

STD 论文的舒适区（KITTI 64-line LiDAR、室外、异构几何）与我们的场景（Airy 稀疏、室内走廊、重复几何）正好相反。

### 4.3 最终决定：Phase 4 取消

v2 plan 的判据：「STD 反而差 → 暂停，重审场景适配性」。判据满足。

**最终路线**：SC + TEASER + submap 多帧累积。剩下的 14% 失败（3 个稀疏 submap）通过两条兜底链：
1. **BBS + TEASER 兜底**（Phase 5 已实现，详见 `analysis/experiments/bbs_airy_22submaps.md`）：触发条件「SC 返回 NO_SOLUTION」→ 调 `examples/global_localization_bbs/bbs_localize`；后置验证 `rot_err < 5° AND inliers ≥ 10`，否则继续往下兜。实测能稳救 1/3 稀疏 case（submap 7）。
2. **multi-submap consensus**：余下仍救不下的 case（submap 10、17 这类极稀疏 < 3500 点）合并相邻 ±1 邻居后重跑 SC 或 BBS。

### 4.4 旧版规划（仅作历史归档）

<details>
<summary>Phase 3 / 4 原计划展开</summary>

原方案是：
1. catkin_make STD，跑 KITTI demo
2. 对比 SC 召回率
3. 若 STD 召回率高 +10pp → 剥离 ROS + 集成

实际：Phase 3 直接剥离 + 在我们数据上实测，0% 准确率终结路线。

</details>

---

## 5. 模块全景图

```
~/Documents/GitHub/tools/
├── TEASER-plusplus/
│   ├── analysis/                          # 6 份代码精读 + 索引
│   │   └── experiments/                   # 实验报告
│   │       ├── sc_airy_22submaps.md       # Phase 2 报告 (19/22 = 86%)
│   │       ├── std_airy_22submaps.md      # Phase 3 报告 (0/22, STD 不适用)
│   │       ├── bbs_airy_22submaps.md      # Phase 5 报告 (20/22 中等阈值, 兜底使用)
│   │       └── isc_airy_23submaps.md      # Phase 6 报告 (ISC ≈ SC 完全打平, 归档)
│   ├── doc/
│   │   ├── airy_localization_plan.md      # 本文档
│   │   └── phase5_bbs_plan.md             # Phase 5 方案细则
│   └── examples/
│       ├── local_refinement/              # ✅ Phase 1 输出
│       │   └── local_refinement.cc
│       ├── global_localization_sc/        # ✅ Phase 2 主线 + Phase 6 ISC 扩展
│       │   ├── global_localization_sc.cc  # 支持 --descriptor sc|isc (DB v3 格式)
│       │   ├── third_party/scancontext/   # vendored, ROS-stripped, + makeIntensityScancontext
│       │   ├── tools/submap_to_pcd.py     # glim_ros submap → PCD + 重力对齐 + intensity
│       │   └── run_sc_isc_eval.py         # 通用 leave-one-out 评测脚本
│       ├── global_localization_std/       # Phase 3 输出（评估失败，归档）
│       │   ├── std_eval_22submaps.cc
│       │   └── third_party/std/           # vendored, ROS-stripped
│       └── global_localization_bbs/       # ✅ Phase 5 输出（兜底）
│           ├── bbs_smoke_test.cc          # bunny self-test (BBS API + 链接)
│           ├── bbs_localize.cc            # BBS coarse + TEASER refine 主程序
│           ├── refine_pipeline.hpp        # FPFH+TEASER 精化（从 Phase 1 抽出）
│           ├── run_22submaps.sh           # leave-one-out 批跑脚本
│           └── third_party/               # vendored KOKIAOKI/3d_bbs（add_subdirectory in-tree 编译）
├── scancontext_tro/                       # 上游源（仅用于参考与 vendor）
└── STD/                                   # 上游源（已 vendor，主线不使用）
```

---

## 6. 性能预期（基于实测 + 理论）

### 需求 2（local_refinement）

| 数据规模 | 耗时 | 精度 |
|---------|------|------|
| Bunny (~1k pts) | ~130 ms | < 1°, < 2 mm（demo 验证） |
| 真实 80m 局部地图（降采样后 ~50k） | 估计 0.5-1.5 s | 目标 < 0.5°, < 0.1 m |

成功率：粗位姿误差 ≤ 10 m 时 > 95%（基于 TEASER 论文与 FPFH 文献）。

### 需求 1（global_localization_sc）

| 阶段 | 估计耗时（200m×200m 地图，~400 anchors）|
|------|-------|
| 离线 build DB | 5-30 s |
| Stage A SC retrieval | 50-200 ms |
| Stage B per-candidate TEASER | 200-500 ms |
| **K=5 端到端** | **~1-3 s** |

成功率（基于 SC 论文）：
- 开阔区域、特征丰富 → > 90%；
- 长走廊 / 对称场景 → 70-80%（**Phase 3 STD 评估的主要动机**）。

---

## 7. 工程使用建议

### 7.1 数据预处理（在调用任一 example 之前）

- **多帧累积**：Airy 单帧太稀，强烈建议用 IMU 去畸变 + 累积 1-2 秒（10-20 帧）后再做 FPFH/SC。
- **动态目标剔除**：人、车、设备在重定位中是干扰，建议用语义/几何先滤掉。
- **地面/天花板抑制**（可选）：极对称的几何（如长走廊里只剩平直地面+天花板）会让 FPFH 描述子退化，按高度筛点能改善。

### 7.2 参数调优起点

| 参数 | 煤矿巷道场景 |
|------|-------------|
| `noise_bound` | 0.05-0.1 m（多帧累积后） |
| `fpfh_normal_radius` | 0.5-1.0 m |
| `fpfh_radius` | 1.0-2.0 m |
| `voxel_size_map` | 0.3 m |
| `voxel_size_scan` | 0.2 m |
| `anchor_step` | 5-10 m |
| `sc_radius` | 80 m（默认） |

### 7.3 失败模式诊断

| 现象 | 可能原因 | 检查项 |
|------|---------|--------|
| inlier_ratio < 5% | FPFH 匹配几乎全错 | 多帧累积是否够？运动畸变是否补偿？ |
| inlier_ratio > 30% 但 RMSE > 0.5m | FPFH 系统偏差 | `--fpfh-r` 是否过大（吃进对侧帮）？ |
| TEASER 报 max-clique <= 1 | 几乎无 inlier-inlier TIM | noise_bound 是否过小？ |
| Quatro 残差 ~1° | scene 含非平凡 roll/pitch | 改 GNC_TLS 或检查 IMU 预对齐 |
| SC 检索都是 dist=0 | scene 尺度远小于 80m | 改 PC_MAX_RADIUS 常量 |

---

## 8. 文档与代码索引

- `analysis/` —— TEASER++ 代码精读 6 篇 + 索引（共 4154 行，154 个 [VERIFY:] 引用）。读 example 代码前建议先看 `analysis/02_DATA_FLOW.md` 与 `analysis/ALGORITHM_03-MaxClique_TIM.md`。
- `examples/local_refinement/README.md` —— Phase 1 用法、参数表、IMU 预对齐说明。
- `examples/global_localization_sc/README.md` —— Phase 2 用法、三个子命令详解、SC 参数边界。
- `examples/teaser_cpp_fpfh/quatro_cpp_fpfh.cc` —— TEASER 官方 FPFH+Quatro 模板（local_refinement 与 global_localization_sc 都从这里裂变出来）。

## 9. 下一步触发器

- ~~**Phase 3 触发条件**~~：已完成（2026-05-11）。
- ~~**Phase 4 触发条件**~~：已取消，路线最终为 SC + TEASER。
- ~~**Phase 5 触发条件**~~：已完成（2026-05-11），作为兜底接入。
- **Hybrid pipeline 集成**：在 `global_localization_sc.cc` 加 `--fallback-bbs` 选项；触发条件 `--fallback-bbs` 打开 + SC 返回 NO_SOLUTION；调 `bbs_localize` 并做后置 `rot<5° AND inliers≥10` 验证。约半天工作量。
- **Multi-submap consensus**：救 submap 10/17 这类极稀疏 case。合并 ±1 邻居后重跑 SC（或 BBS）。
- **回头看 Phase 1/2**：用真实 Airy 数据跑一遍 demo，校准默认参数。
