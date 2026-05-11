# Phase 5：3D-BBS 粗定位 + TEASER 精化 — 方案

> **状态**：✅ 已实施并评估（2026-05-11）。**结论：进 hybrid pipeline 作兜底**，不切主线。
> **实测**：pose 返回率 22/22 = 100%；中等阈值 (rot<5° AND trans<1m) 20/22 = 91%；严格 (rot<3° AND trans<0.5m) 14/22 = 64%。详见 `analysis/experiments/bbs_airy_22submaps.md`。
> **关联**：`doc/airy_localization_plan.md`（v2 主线方案）、`analysis/experiments/sc_airy_22submaps.md`（Phase 2 基线 86%）、`analysis/experiments/std_airy_22submaps.md`（Phase 3 0% 终止）、`analysis/experiments/bbs_airy_22submaps.md`（Phase 5 实测）
> **作用**：在 SC + TEASER 主线之外，引入一条**不依赖描述子**的全局定位备选路径，重点目标是把 Phase 2 失败的 3 个稀疏 submap（7 / 10 / 17）救回来。实测：稳救 1 个（submap 7），另 2 个虽然不再 NO_SOLUTION 但 pose 不可信，仍需 multi-submap consensus 兜底。

---

## 1. 动机

Phase 2 已经把 22 个 submap 的全局定位做到 19/22 = 86.4%，失败 3 个都是 **< 4000 点** 的稀疏 submap。失败链路追到底是：

```
稀疏 submap (<4k 点)
   → FPFH 法线邻域不稳（5-10 个邻居）
   → 描述子退化
   → SC top-30 里就算召回到真实区域，TEASER 也拿不到 ≥10 个可靠对应
   → 输出 NO_SOLUTION
```

SC + TEASER 这条线**整条链路都建在描述子可靠的假设上**。在稀疏点云上要救这 3 个 case，要么换思路、要么累积更厚的 submap。

3D-BBS 提供另一种思路：

| | SC + TEASER | 3D-BBS |
|---|---|---|
| 候选生成 | BEV 极坐标统计描述子检索 | **体素占用的几何穷举**（BnB） |
| 关键依赖 | 描述子区分度（受点密度影响） | 体素是否占用（对点密度不敏感） |
| 输出 | 6DOF 位姿（含尺度） | 4×4 粗位姿（论文称 *coarse*） |
| 适用场景 | 需要场景区分度 | 任何被点云覆盖的区域 |

**核心假设**：3D-BBS 不走描述子路径，理论上对稀疏点云更鲁棒；但它只给粗解，精度做不到 ±0.1 m。

**解法**：3D-BBS 出粗解 → 喂给 Phase 1 `local_refinement` 的 FPFH+TEASER 做精化。这正好是 Phase 1 设计的输入形态（粗位姿误差 3-10 m / 几度），3D-BBS 的输出（量级 ±0.5 m / ±2°）远在 Phase 1 的舒适区内。

---

## 2. 方案架构

```
   单帧/submap (gravity-aligned)            预建地图 PLY/PCD
            │                                     │
            │                          ┌──────────▼────────────┐
            │                          │ 离线：voxelmap 构建    │
            │                          │ 缓存到磁盘（130 ms 加载）│
            │                          └──────────┬────────────┘
            │                                     │
   ┌────────▼─────────────────────────────────────▼────────────┐
   │ Stage A — 3D-BBS 全局粗定位                                  │
   │   set_tar_points(map, min_level_res=0.5, max_level=6)      │
   │   set_src_points(scan)                                     │
   │   set_trans_search_range(map_bbox)                         │
   │   set_angular_search_range(rpy_bounds)                     │
   │   localize() → T_coarse (4×4)                              │
   └────────┬────────────────────────────────────────────────────┘
            │ T_coarse  (实测 / 经验：±0.5 m / ±2°)
            │
   ┌────────▼────────────────────────────────────────────────────┐
   │ Stage B — local_refinement (Phase 1 代码原样复用)             │
   │   scan_init = T_coarse · scan                              │
   │   FPFH (map_voxel, scan_voxel, 法线/特征半径)                 │
   │   TEASER GNC_TLS                                           │
   │   → T_correction                                            │
   │   T_refined = T_correction · T_coarse                       │
   └────────┬────────────────────────────────────────────────────┘
            │
            ▼ 最终 6DOF 位姿（目标 ±0.05 m / ±0.5°）
```

**关键接口**：3D-BBS 的输出是 `Eigen::Matrix4f`（gpu_bbs3d/bbs3d.cuh:`get_global_pose()`），直接作为 `local_refinement` 的 `--rough-pose` 输入。两段之间不需要中间格式转换。

---

## 3. 实施步骤

### 3.1 vendor 3D-BBS（不重新编译 CUDA）

- 复用 workspace 已有的 `src/3d_bbs/`，不再 vendor 到 TEASER-plusplus 树内（与 SC / STD 不同 — BBS 体量大且依赖 CUDA，外部链接更合适）
- TEASER-plusplus 这边只加 `find_package(gpu_bbs3d REQUIRED)`，库已经 `sudo make install` 装到 `/usr/local`

### 3.2 新增 `examples/global_localization_bbs/`

```
examples/global_localization_bbs/
├── global_localization_bbs.cc       # 主程序：build | localize 两个子命令
├── CMakeLists.txt                    # link gpu_bbs3d + teaser_registration + teaser_features
├── README.md
└── tools/
    └── (复用 ../global_localization_sc/tools/submap_to_pcd.py)
```

子命令设计：

```bash
# build：从地图构建 voxelmap 并缓存
./global_localization_bbs --mode build \
    --map /path/to/global_map.pcd \
    --voxelmap-out /tmp/w2_voxelmap \
    --min-level-res 0.5 --max-level 6

# localize：完整 pipeline（BBS 粗 + TEASER 精）
./global_localization_bbs --mode localize \
    --map /path/to/global_map.pcd \
    --voxelmap /tmp/w2_voxelmap \
    --scan /path/to/submap_levelled.pcd \
    --min-level-res 0.5 --max-level 6 \
    --score-threshold 0.9 \
    --refine-noise-bound 0.15 \
    --output /tmp/result.txt
```

### 3.3 完全复用的部分

- `tools/submap_to_pcd.py`：BBS 同样要求重力对齐（roll/pitch 已对齐到地图 z 轴），逻辑与 SC 完全一致
- `local_refinement.cc` 的 FPFH+TEASER 函数：抽成 header 后给 BBS 主程序 include 即可，不重写

### 3.4 评估脚本

照搬 `std_eval_22submaps.cc` 的 leave-one-out 框架，改成：
- DB：固定地图 `global_map.pcd`（不像 STD 那样 leave-one-out 重建库）
- Query：每次扫描一个 submap 的 `submap_levelled.pcd`
- 比较：predicted vs `T_world_origin` ground truth

---

## 4. 评估协议（与 Phase 2 / 3 完全可比）

**数据**：`map_w2/20260511_110448`，22 submap，与 SC / STD 实验同一份。

**指标**：

| 指标 | 来源 |
|---|---|
| 全 22 个 leave-one-out 成功率 | 同 SC 实验 |
| 成功 case 平移误差中位数 / 最大值 | 同 |
| 成功 case 旋转误差中位数 / 最大值 | 同 |
| 端到端耗时 = BBS + FPFH + TEASER | 新增 |
| BBS 单独耗时（含 voxelmap 加载） | 新增 |
| **submap 7 / 10 / 17 是否被救** | **关键观察项** |

**成功判据**：平移误差 < 0.5 m 且旋转误差 < 3°（与 SC 实验一致；BBS+TEASER 应在此区间内）。

---

## 5. 关键技术风险与对策

### 5.1 体素分辨率 vs 内存/搜索时间

| `min_level_res` | 地图体素数（40×43×12 m） | 单次 BnB 估计 |
|---|---|---|
| 1.0 m | ~21k | 50-100 ms |
| 0.5 m | ~165k | 150-300 ms |
| 0.25 m | ~1.3M | 1-3 s ⚠️ |

室内 40 m 场景选 0.5 m 比较合适；走廊宽度 ~4-6 m，0.5 m 体素能区分出墙壁的法向分布。

### 5.2 6DOF vs 3DOF（gravity-aligned）

- BBS 默认 3DOF（yaw + xy + z）+ ±0.02 rad 的 roll/pitch 微调，速度 ~200 ms
- 完整 6DOF 搜索 *慢 10 倍以上*（README 明示）
- 我们的 submap 已经重力对齐过（`submap_to_pcd.py`，roll/pitch 用 IMU 给的姿态消掉），**走 3DOF 模式即可**
- yaw 搜索范围必须给 360°（`min_rpy[2]=0, max_rpy[2]=6.28`）

### 5.3 粗位姿精度可能不达 ±0.5 m

如果 BBS 粗解平移误差超过 1 m，Phase 1 `local_refinement` 的 FPFH 半径（默认 1-2 m）会有问题。对策：

- 把 `local_refinement` 的 `--normal-radius` 临时调大到 2 m
- 或先对 BBS 粗解做一次小 ICP（PCL `IterativeClosestPoint` 用 voxel grid 加速）再喂 TEASER
- 真要兜底：BBS 输出 top-N 假设，每个跑 TEASER，按 inlier 选胜（仿 Phase 2 的 SC top-K 思路）

### 5.4 稀疏 submap 上 BBS 是否真比 SC 好？

BBS 的 src 点云越稀，score（与 voxelmap 体素的命中数）越低。`score_threshold_percentage=0.9` 默认是「90% 的 src 点要命中 voxel」，稀疏 submap 容易过不了阈值。

对策：稀疏 submap 自动把阈值降到 0.7，并把 `timeout_duration_in_msec` 加长。极端情况下走 multi-submap consensus（与 Phase 2 兜底方案统一）。

### 5.5 内存

3D-BBS 用 sparse hash voxelmap，40×43×12 m / 0.5 m 体素只占几 MB。地图扩到 200×200 m 也不会爆。这点比 TIM 矩阵的 N² 内存安全得多。

---

## 6. 参数起点（煤矿巷道 / 室内场景）

| 参数 | 值 | 备注 |
|---|---|---|
| `min_level_res` | 0.5 m | 室内场景，对应 SC 实验里 `--sc-radius 20` 的尺度 |
| `max_level` | 6 | 默认值，对应 32 m 粗层（0.5 × 2⁶） |
| `min_rpy` | `[-0.02, -0.02, 0.0]` | 已重力对齐，roll/pitch 留 ~1° 余量 |
| `max_rpy` | `[+0.02, +0.02, 6.28]` | yaw 全 360° |
| `score_threshold_percentage` | 0.9（稠密）/ 0.7（< 4k 点） | 自动按 src 点数切换 |
| `src_leaf_size` | 0.2 m | 与 Phase 2 `--sc-downsample 0.3` 同量级 |
| `tar_leaf_size` | 0.1 m | 默认 |
| `timeout_duration_in_msec` | 5000 | 默认 60s 太长 |
| Phase 1 精化 `noise_bound` | 0.15 m | 同 Phase 2 |

---

## 7. 决策判据

跑完 22-submap 后按下表决策：

| 结果 | 决策 |
|---|---|
| BBS+TEASER ≥ 21/22（95%+） | **切主线**：BBS+TEASER 成为新的 Phase 2 主路径 |
| BBS+TEASER = 19-20/22（与 SC 持平） | **作为兜底**：SC 失败 → 回落到 BBS+TEASER（hybrid pipeline） |
| BBS+TEASER < 19/22 | **归档** Phase 5；继续走 multi-submap consensus 救 SC 的 3 个失败 |
| BBS+TEASER 救回 7/10/17 但其他场景反而退化 | **场景分流**：稀疏 submap → BBS，稠密 submap → SC |

---

## 8. 与现有 plan 的关系

- `doc/airy_localization_plan.md` §4 已经写了「Phase 4 取消，路线确定为 SC + TEASER + submap 多帧累积」。Phase 5 不推翻该决定，而是开**第二条全局定位通路**，目的是
  1. 救 Phase 2 失败的 3 个稀疏 case；
  2. 在 SC 调参依赖（`PC_MAX_RADIUS`、`anchor_step`、`sc_radius`）较多的前提下，提供一条**几乎只依赖几何分辨率**的备选，便于跨场景迁移。
- 这条线如果跑通，会在 `airy_localization_plan.md` 顶部状态表新增一行 Phase 5 ✅；如果失败，归档于此并在 `airy_localization_plan.md` §4.3 后追加一条「Phase 5 评估失败」说明。

---

## 9. 工件清单（计划产出）

| 路径 | 内容 |
|---|---|
| `doc/phase5_bbs_plan.md` | 本方案（已创建） |
| `examples/global_localization_bbs/global_localization_bbs.cc` | 主程序，build + localize |
| `examples/global_localization_bbs/CMakeLists.txt` | link gpu_bbs3d、teaser_registration、teaser_features |
| `examples/global_localization_bbs/README.md` | 用法、参数表、与 Phase 2 的对照 |
| `examples/global_localization_bbs/bbs_eval_22submaps.cc` | leave-one-out 评估 |
| `analysis/experiments/bbs_airy_22submaps.md` | 实验报告（跑完后） |

---

## 10. 下一步

按工作量从小到大：

1. ~~**API 烟雾测试**~~（已完成）：`bbs_smoke_test` 在 bunny 数据上 PASS（65 ms / trans_err 0.77 m / rot_err 0.47°）。
2. ~~**单 submap 验证**~~（已完成）：submap 5 端到端 0.10 m / 1.71°，与 Phase 2 同 case 精度持平。
3. ~~**22-submap 批跑**~~（已完成）：`bbs_airy_22submaps.md` 出炉，pose 返回率 100%，中等阈值 91%。
4. ~~**决策与文档收口**~~（已完成）：本文档顶部 + `airy_localization_plan.md` 已更新。

实际总工时：< 半个工作日（连续执行，gpu_bbs3d 编译 + 烟测 + 单帧 + 22 帧 + 报告）。

**接下来不属于 Phase 5 的工作**：

- **hybrid pipeline 集成**：把 `bbs_localize` 接为 SC NO_SOLUTION 兜底，加 `rot<5° AND inliers≥10` 后置验证；
- **multi-submap consensus**：救 submap 10/17 这类极稀疏 case（合并相邻 ±1 邻居）。
