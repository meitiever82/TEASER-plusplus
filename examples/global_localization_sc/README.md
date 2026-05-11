# global_localization_sc —— Scan Context 全局检索 + TEASER 精配

对应 `airy_localization_plan.md` 的**需求 1**：完全没有先验位姿，只有一张全局地图和一帧扫描，要求恢复扫描在地图中的 6DOF 位姿。

两阶段方案：
1. **Stage A** —— Scan Context 检索：离线把地图按网格采样若干 anchor，每个 anchor 算一个 SC 描述子存库。在线时把扫描算成 SC，与库内 anchor 比对，取 top-K。
2. **Stage B** —— TEASER 精配：对每个候选用 `local_refinement` 同款 FPFH+TEASER 流程跑一次，按 inlier 数与 RMSE 选胜者。

## 编译

需要 `BUILD_TEASER_FPFH=ON`：

```bash
cmake -S . -B build -DBUILD_TEASER_FPFH=ON
cmake --build build -j
```

可执行文件：`build/examples/global_localization_sc/global_localization_sc`。

## 三个子命令

### `--mode build` 构建 SC 库（离线一次）

```bash
./global_localization_sc --mode build \
    --map /path/to/map.ply \
    --db /path/to/sc_database.bin \
    --anchor-step 10.0 \
    --sc-radius 80.0 \
    --sc-downsample 0.5
```

把地图按 `--anchor-step` 米的 XY 网格采样 anchor，每个 anchor 抽取半径 `--sc-radius` 米的局部点云，平移到原点后送 SC 算描述子。结果写入二进制 DB 文件。

**关键参数**：

| 参数 | 默认 | 含义 |
|------|------|------|
| `--anchor-step` | 10 | 网格步长（米）。煤矿巷道宽 4-6m 时建议 5；露天大场地可 15-20 |
| `--sc-radius` | 80 | **必须与 `SCManager::PC_MAX_RADIUS` 一致**。是 SC 内部硬编码常量，改了得动 `Scancontext.h` |
| `--sc-downsample` | 0.5 | SC 输入降采样体素（米）。SC 是形状描述子，密度变化对它是噪声 |

### `--mode query` 仅做 SC 检索（调试用）

```bash
./global_localization_sc --mode query \
    --db sc_database.bin \
    --scan scan.ply \
    --top-k 5 \
    --sc-downsample 0.5
```

只跑 Stage A，打印 top-K 候选 anchor 位置 + SC 距离 + yaw 估计。**不跑 TEASER**，用来调 SC 参数或排查检索召回率。

### `--mode localize` 端到端定位

```bash
./global_localization_sc --mode localize \
    --map map.ply \
    --db sc_database.bin \
    --scan scan.ply \
    --top-k 5 \
    --algo Quatro \
    --noise-bound 0.1 \
    --output refined_pose.txt
```

完整 Stage A + Stage B。对每个 top-K 候选独立跑一次 FPFH+TEASER，最后输出 inlier 数最大、且通过 `--min-inliers` 阈值的解。

#### 全部可调参数

| 参数 | 默认 | 含义 |
|------|------|------|
| `--top-k` | 5 | 进入 Stage B 的候选数。增大召回率但耗时线性增长 |
| `--min-inliers` | 10 | inlier 绝对数下限，低于此值的候选视为失败 |
| `--algo` | Quatro | rotation solver。Quatro 仅 yaw，TEASER 完整 SO(3) |
| `--local-radius` | 60 | TEASER 阶段局部地图半径（米） |
| `--voxel-map` | 0.3 | TEASER 地图体素 |
| `--voxel-scan` | 0.2 | TEASER 扫描体素 |
| `--fpfh-normal-r` | 1.0 | FPFH 法线半径 |
| `--fpfh-r` | 2.0 | FPFH 描述子半径 |
| `--noise-bound` | 0.1 | TEASER 噪声上界 |
| `--max-corres` | 3000 | 对应数硬上限（防 N² 爆内存） |

### 输入输出格式

- **map / scan**：标准 PLY，含 `x y z`。
- **DB**：二进制 `sc_database.bin`，由本程序自己读写，格式见 `global_localization_sc.cc` 中 `SCDatabase::save/load`。
- **refined_pose.txt** 输出：
  ```
  # winning anchor: 47
  # refined pose: x y z roll pitch yaw (radians)
  12.34 56.78 1.23 0.01 -0.02 1.57
  # 4x4 matrix:
  [R t]
  [0 1]
  # quality: inlier_ratio rmse_meters elapsed_seconds
  0.42 0.075 1.456
  ```

## 流程图

```
                          ┌─ build 模式 ──────────────────────────┐
                          │                                       │
                          │ map.ply                               │
                          │   │                                   │
                          │   ├─ grid sample (anchor_step)        │
                          │   │   for each (x, y):                │
                          │   │     local = crop(map, x,y, R)     │
                          │   │     local -= (x, y, 0)            │
                          │   │     SC = SCManager.makeSC(local)  │
                          │   │     ringkey = makeRingkey(SC)     │
                          │   │     push to DB                    │
                          │   ▼                                   │
                          │ sc_database.bin                       │
                          └───────────────────────────────────────┘

                          ┌─ localize 模式 ───────────────────────┐
                          │                                       │
                          │ scan.ply  + map.ply  + DB             │
                          │   │                                   │
                          │   ├─ Stage A: SC 检索                  │
                          │   │   SC_q = makeSC(scan)             │
                          │   │   for each anchor:                │
                          │   │     d = distanceBtnScanContext()  │
                          │   │   top-K by sc_distance asc        │
                          │   │                                   │
                          │   ├─ Stage B: per-candidate TEASER    │
                          │   │   for each (anchor_pos, yaw):     │
                          │   │     T_rough = (anchor_pos, yaw)   │
                          │   │     scan_init = T_rough * scan    │
                          │   │     local = crop(map, anchor_pos) │
                          │   │     FPFH + Matcher + TEASER       │
                          │   │     -> (T_correction, quality)    │
                          │   │     T_final = T_correction*T_rough│
                          │   │                                   │
                          │   ▼  argmax_quality                   │
                          │ refined_pose.txt                      │
                          └───────────────────────────────────────┘
```

## 性能预期

基于 `SCManager` 默认参数（`PC_NUM_RING=20, PC_NUM_SECTOR=60, PC_MAX_RADIUS=80m`）：

| 阶段 | 200m × 200m 地图、anchor_step=10m（~400 anchors）|
|------|----|
| build DB（离线一次） | ~5-30 s |
| SC 检索 1 query（暴力对比所有 anchor） | ~50-200 ms |
| 单候选 TEASER 验证 | ~200-500 ms |
| **端到端定位（K=5 候选）** | **~1-3 s** |

大幅大于 400 anchors 时可以加 KD 树用 ringkey 粗筛——目前用暴力，原因是 anchor 数 ≤ 1000 时暴力反而省事且更稳。

## SC 默认参数 vs 室内数据

`SCManager` 的 `PC_MAX_RADIUS = 80m` 是为**车载 LiDAR**调的。如果你的场景明显小于这个量级（如 3DMatch 室内 ~3m 房间），**SC 没区分度**——所有 anchor 看到的内容几乎一样，距离都是 0。

煤矿巷道场景的尺度大致是：
- 巷道断面 4×4–6×6 米；
- 单帧 Airy 有效半径 ~30–80 米；
- 长大巷道沿向 100+ 米。

这正是 SC 的舒适区。如果你的场景偏小，要改 `third_party/scancontext/Scancontext.h` 的 `PC_MAX_RADIUS` 常量（属于硬编码），并保持 `--sc-radius` 一致。

## 与 `local_refinement` 的关系

Stage B 的代码（`verifyOneCandidate`）和 `local_refinement` 的 `refineLocalPose` 几乎相同——都是「rough pose → local crop → FPFH → TEASER → compose」。两个 example 各自独立成文件是保留 example 间无依赖的库内惯例，不是技术约束。

如果你要在生产代码中合并，建议把 `extractLocalMap` / `applyTransform` / `voxelDownsample` / `verifyOneCandidate` 抽到一个共享 header。

## 已知限制 & 后续

1. **SC 在长走廊场景召回率会下降**——这就是为什么 `airy_localization_plan.md` 的 Phase 3-4 评估 STD（基于关键点三角形，对走廊更鲁棒）。
2. **单帧 Airy 太稀**——建议先用 IMU 累积 1–2 秒再调用本程序，提高 SC 与 FPFH 的稳定性。
3. **没有 multi-frame consensus**——一次定位失败时没有再尝试机制，需要上层逻辑（如多帧投票）。
4. **DB 没有时间戳/版本管理**——重建地图后需要手动重新跑 `--mode build`。
5. **`SC_DIST_THRES = 0.13`** 在 `SCManager` 中是硬编码常量，但本 example 没有使用 `detectLoopClosureID`，所以这个阈值不生效。我们用「top-K + TEASER 验证」代替了「单阈值判定」，更鲁棒。

## 上游致谢

`third_party/scancontext/` 整套 vendor 自 [gisbi-kim/scancontext_tro](https://github.com/gisbi-kim/scancontext_tro)：
- Giseop Kim et al., "Scan Context: Egocentric Spatial Descriptor for Place Recognition within 3D Point Cloud Map", IROS 2018
- Giseop Kim et al., "Scan Context++: Structural Place Recognition Robust to Rotation and Lateral Variations in Urban Environments", IEEE T-RO 2022

vendor 时仅去掉了 .h 中**未使用**的 ROS / OpenCV / cv_bridge / pcl_conversions includes，算法主体未改。
