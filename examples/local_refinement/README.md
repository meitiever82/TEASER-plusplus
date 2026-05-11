# local_refinement —— 基于粗位姿的精确配准

对应 `airy_localization_plan.md` 的**需求 2**：给定地图 + 单帧扫描 + 粗位姿（3-10m / 几度误差），用 FPFH + TEASER（默认 Quatro）精化到 ±0.1m / ±0.5°。

## 编译

需要 `BUILD_TEASER_FPFH=ON`（因为依赖 `teaser_features` 的 FPFH + Matcher）。在仓库根目录：

```bash
cmake -S . -B build -DBUILD_TEASER_FPFH=ON
cmake --build build -j
```

可执行文件输出在 `build/examples/local_refinement/local_refinement`。

## 运行

### 演示模式（无参数）

加载 `example_data/bun_zipper_res3.ply` 作为地图与扫描，应用合成的真值变换，加随机噪声形成粗位姿，跑全流程，最后打印与真值的误差。

```bash
cd build/examples/local_refinement
./local_refinement
```

预期输出：旋转误差 < 1°，平移误差 < 1cm（兔子尺度）。

### 真实数据模式

```bash
./local_refinement \
    --map /path/to/map.ply \
    --scan /path/to/scan.ply \
    --rough-pose /path/to/rough_pose.txt \
    --output refined_pose.txt
```

#### 输入文件格式

- **`map.ply`** / **`scan.ply`**：标准 PLY，仅需要 `x y z` 三列。XYZ 单位米。
- **`rough_pose.txt`**：一行 6 个数（弧度）：
  ```
  x y z roll pitch yaw
  ```
  约定 `R = Rz(yaw) · Ry(pitch) · Rx(roll)`。

#### 关键参数

| 参数 | 默认 | 推荐范围 / 含义 |
|------|------|---------------|
| `--algo` | `Quatro` | `Quatro`（仅 yaw，IMU 已对齐时用）或 `TEASER`（完整 SO(3)） |
| `--local-radius` | 60 | 粗位姿误差 + 扫描半径 + 安全余量。误差 3m 建议 50，误差 10m 建议 80 |
| `--voxel-map` | 0.3 | 地图降采样体素（米）。地图密集时增大可省时 |
| `--voxel-scan` | 0.2 | 扫描降采样体素 |
| `--fpfh-normal-r` | 1.0 | FPFH 法线邻域。Airy 稀疏点云推荐 0.5-1.0，密集 LiDAR 推荐 0.1-0.3 |
| `--fpfh-r` | 2.0 | FPFH 描述子半径，一般取 normal-r 的 2 倍 |
| `--noise-bound` | 0.1 | TEASER 噪声上界。**最关键参数**——多帧累积建议 0.05-0.1，单帧建议 0.1-0.2 |
| `--max-corres` | 3000 | 对应数上限。TEASER TIM 矩阵 O(N²)，N=3000 时约 100MB，N=5000 时约 300MB |

#### 输出文件 `refined_pose.txt`

```
# refined pose: x y z roll pitch yaw (radians)
1.234 2.345 0.123 0.012 -0.034 1.567
# 4x4 matrix:
[R t]
[0 1]
# quality: inlier_ratio rmse_meters elapsed_seconds
0.45 0.067 0.834
```

## 工作流

```
rough_pose ──┐
             ▼
  ┌──────────────────────┐
  │ 1. extractLocalMap   │  以 rough_pose 为中心裁地图
  │    (radius)          │
  └──────────┬───────────┘
             ▼
  ┌──────────────────────┐
  │ 2. applyTransform    │  scan_in_map = T_rough · scan
  │    scan -> map frame │
  └──────────┬───────────┘
             ▼
  ┌──────────────────────┐
  │ 3. voxelDownsample   │  两片各自降采样
  └──────────┬───────────┘
             ▼
  ┌──────────────────────┐
  │ 4. FPFH 特征         │  两片各算 FPFH
  └──────────┬───────────┘
             ▼
  ┌──────────────────────┐
  │ 5. Matcher           │  cross-check + tuple test
  │    (correspondences) │
  └──────────┬───────────┘
             ▼
  ┌──────────────────────┐
  │ 6. cap at N_max      │  防止 TEASER O(N²) 爆内存
  └──────────┬───────────┘
             ▼
  ┌──────────────────────┐
  │ 7. TEASER solve      │  Quatro 或 GNC_TLS
  │    -> T_correction   │
  └──────────┬───────────┘
             ▼
  T_refined = T_correction · T_rough
```

## 关于 Quatro / IMU 预对齐

Quatro 假设 **roll = pitch = 0**——它只在 XY 平面估 yaw，不动 Z 旋转。如果场景有显著的 roll/pitch（如煤矿带坡度的巷道），**Quatro 会失真**。

正确用法：
- 用 IMU 估出 `R_rp = Ry(pitch_imu) · Rx(roll_imu)`；
- 把 `rough_pose` 写成 `T_rough.R = R_yaw · R_rp`；
- 程序里 scan 被 `T_rough` 预乘后，**roll/pitch 已经对齐到地图坐标系**；
- Quatro 估的「correction」只剩 yaw 微调；
- 最终 `T_refined = T_correction · T_rough` 自动复合上 R_rp。

如果对 IMU 不放心，直接 `--algo TEASER`（完整 SO(3)），代价是慢 2-3 倍且 PMC_EXACT 在大 N 时显著变慢。

## 已知限制

1. **Airy 稀疏点云**：单帧扫描点数少（~10k-50k），FPFH 描述子稳定性差。**强烈建议** 1-2 秒多帧累积后再调用本程序。
2. **正态分布假设**：FPFH 在动态目标（人、车）上特征不稳定。预处理时建议移除运动目标。
3. **大地图加载**：`PLYReader` 全量加载到内存。地图 > 几亿点时需要预切分。
4. **`max-corres` 取舍**：截断到 3000 时极少数情况会漏掉真 inlier。如果质量指标差，先调大 `--max-corres` 再考虑别的参数。

## 质量评判

- **`inlier_ratio` > 30% 且 `rmse` < 0.15m**：通常可用；
- **`inlier_ratio` < 15%**：粗位姿太差或场景缺特征，建议改用全局定位（需求 1，STD + TEASER）；
- **`rmse` > 0.5m 但 `inlier_ratio` 高**：FPFH 匹配出现系统偏差（特征半径不当），调 `--fpfh-r`。
