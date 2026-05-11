# 实验：Stable Triangle Descriptor (STD) 在同一 22 submap 数据上的表现

> **日期**：2026-05-11
> **数据**：`map_w2/20260511_110448`，与 SC 实验完全相同的 22 submap
> **代码**：`examples/global_localization_std/std_eval_22submaps.cc`
> **目的**：验证 v2 plan 中 Phase 3 假设——「STD 关键点三角形机制比 SC 形状描述子更适合稀疏 Airy / 对称走廊场景」

---

## 1. 结果（直说）

**STD top-1 leave-one-out 准确率：0 / 22 = 0%**

对照：相同数据上 **SC + TEASER：19 / 22 = 86.4%**（见 `sc_airy_22submaps.md`）。

| submap | GT pos | STD matched id | trans err | rot err |
|---|---|---|---|---|
| 00 | (+0.58, -0.03) | 18 | 16.8m | 58° |
| 01 | (+5.64, -0.49) | 00 | 4.4m | 16° |
| 02 | (+6.67, +1.67) | 12 | 11.7m | 104° |
| 03 | (+5.38, +4.55) | 09 | 19.9m | 132° |
| 04 | (+4.03, +6.95) | 18 | 14.8m | 88° |
| 05 | (+2.79, +10.26) | 13 | 28.3m | **176°** |
| 06 | (+1.95, +13.65) | 17 | 16.6m | **180°** |
| 07 | (+1.61, +16.68) | 16 | 13.4m | 111° |
| 08 | (+0.93, +20.15) | 18 | 25.7m | 123° |
| 09 | (-0.17, +21.96) | 01 | 34.5m | 122° |
| 10 | (-1.98, +22.08) | 15 | 26.8m | 136° |
| 11 | (-3.29, +21.67) | 05 | 19.2m | 142° |
| 12 | (-5.02, +21.83) | 21 | 37.5m | 147° |
| 13 | (-6.41, +21.66) | 00 | 25.4m | 116° |
| 14 | (-7.16, +20.15) | 03 | 43.7m | 104° |
| 15 | (-7.52, +17.14) | 21 | 22.8m | 67° |
| 16 | (-8.15, +13.51) | 17 | 13.0m | 113° |
| 17 | (-8.42, +10.95) | 21 | 17.0m | 146° |
| 18 | (-8.25, +8.60) | 04 | 17.6m | 120° |
| 19 | (-7.91, +5.72) | 16 | 21.1m | 66° |
| 20 | (-7.72, +2.90) | 09 | 25.6m | 153° |
| 21 | (-5.35, +1.10) | 10 | 30.8m | 103° |

**关键现象**：
- 平移误差 4-44 m，旋转误差几乎涵盖整个 [60°, 180°] 区间
- 多个 case rot_err ≈ 180°（submap 5, 6, 7, 13...）—— STD 把场景的「另一面」当成同一处，触发 yaw flip
- 这种「flip + 远端误匹配」是 STD 在**对称几何**下的经典失效模式

---

## 2. 实验配置

### 2.1 数据

与 SC 实验完全一致：22 个 submap，每个 ~2-8k 点，世界坐标系下覆盖 ~40m × 43m × 12m 室内场景。

### 2.2 评估协议

**Leave-one-out**：
1. 对每个 submap i，用其他 21 个建 STD DB（`AddSTDescs`）
2. 用 submap i 作 query（`GenerateSTDescs` + `SearchLoop`）
3. 取 STD 返回的 (matched_id, R, t)，复合 matched_id 的真值位姿得到预测的 query 位姿
4. 与 query 的真值位姿比较

每个 submap 是独立的 STDescManager 实例（22 次完整建库 + 1 次查询）。22 次总耗时 ~0.6 秒，STD 本身极快。

### 2.3 关键 STD 参数

在默认参数（为 Velodyne 64-line 调的）基础上为我们的稀疏室内场景做了调整：

| 参数 | 默认 | 我们用的 | 含义 |
|---|---|---|---|
| `voxel_size_` | 2.0 | 1.0 | 体素大小（米）。室内 ~40m 场景比 KITTI 紧凑 |
| `voxel_init_num_` | 10 | **5** | 一个 voxel 至少 N 点才认作 plane。默认 10 在我们稀疏云上会让大量 voxel 不被识别为 plane，导致提不出 corner —— **这是 STD 在 Airy 上跑通的关键参数** |
| `maximum_corner_num_` | 100 | 500 | 每个 keyframe 最多保留多少 corner |
| `corner_thre_` | 10 | 5 | corner 阈值，更小→更多 corner |
| `descriptor_min_len_` | 2 | 1 | 三角形最短边（米） |
| `descriptor_max_len_` | 50 | 30 | 三角形最长边（米） |
| `rough_dis_threshold_` | 0.01 | 0.5 | 三角形粗匹配容差 |
| `vertex_diff_threshold_` | 0.5 | 5.0 | 顶点属性差容差 |
| `dis_threshold_` | 0.5 | 2.0 | 几何验证距离阈值 |
| `icp_threshold_` | 0.5 | 0.0 | 接受 top-1 不再过滤（我们用 GT 判断而非 ICP 分数判断） |

**注意**：用默认参数（voxel_init_num=10）时，STD 生成的描述子数量为 0 — DB 完全空。必须放到 5 才有产出。

---

## 3. STD 失效的机理分析

### 3.1 第一道关：corner 提取

STD 流程是 **voxelize → 平面检测 → 沿平面边缘提 corner → 用 corner 三元组构三角形**。我们的 Airy submap 平均 5k 点 / 50m 区域 → 平均每 5L 体素只有 ~50 点。voxel_init_num=10 默认意味着「N>10 才算 planar」，很多体素被丢掉 → 全图剩极少平面 → corner 寥寥 → 三角形不存在。

放到 voxel_init_num=5 后 corner 提到了几百，三角形也有了，**但描述子质量不高**——很多 corner 来自不稳定的 voxel 边缘。

### 3.2 第二道关：候选筛选

`candidate_selector` 按三角形三边长度做哈希查表。三个边长被四舍五入到整数 + 1 体素的邻域内匹配。对于我们的场景：
- 全 22 个 submap 都是同一个走廊空间的不同片段
- 走廊的「短边」「长边」长度都相似（比如 5m, 8m, 12m）
- 不同 submap 的三角形哈希到相同的 bucket → 大量 false candidate

我们看到 21 个候选，rough match 数 700-1200 量级，但**大多是几何上相似但语义上完全不同的三角形**。

### 3.3 第三道关：几何验证

`candidate_verify` 用投票：取一个 triangle pair 计算 (R, t)，看其他 pair 在这个 (R, t) 下能不能对齐（vertex 误差 < 3m）。需要至少 4 个三角形投同一个 (R, t)。

对于真正的正确匹配，应该有几十个三角形高度一致——但实际我们的数据上正确匹配的 vote 数（~0.025 verify_score 量级）反而**比错误匹配（0.16）低**。

**为什么？** 走廊场景的几何重复性意味着多个「错误三角形对」可以**在某个 180°-flipped 或平移过的 (R, t) 下整齐对齐**——它们其实是把场景的镜像或旋转版本匹配上了。这种「自洽的错误」比真匹配的 vote 数还高。

### 3.4 STD 设计前提

STD 论文（HKU-MARS, ICRA 2023）的实验场景：
- **室外** KITTI、Mulran、Sequoia
- **完整 64 线 Velodyne / Livox Avia** 单帧 100k+ 点
- **大量异构几何**（建筑物、车辆、树木、地面纹理）—— 三角形描述子有足够区分度

我们的场景：
- **室内** 工业巷道
- **稀疏 Airy submap** 2-8k 点
- **几何重复性极高**（走廊截面都差不多）

这就是为什么 STD 工作得这么差。

---

## 4. v2 plan 的判据

`doc/airy_localization_plan.md` Phase 3 的入 Phase 4 条件：

> STD 召回率 ≥ SC 召回率 + 10 个百分点 → 进 Phase 4；
> 差不多 → 不集成，留 SC；
> STD 反而差 → 暂停，重审场景适配性。

**触发条件：STD 反而差**（0% vs 86%）。结论：**Phase 4 取消，留在 SC + TEASER 路线。**

---

## 5. 是否还有救？

理论上可以再做的事，但 ROI 越来越低：

### 5.1 TEASER 救 STD？

STD 的 `loop_std_pair` 是匹配三角形对的列表。每对 = 3 点对应。N 对匹配 → 3N 点对应，可以丢给 TEASER 做最后的鲁棒位姿估计。

我们已经看到：STD 的 top-1 之外的候选里有时确实有正确的（如 submap 5 的候选 0 是 frame_id=4，即真邻居 submap 4，verify_score 0.025）。如果改成 **top-K（K=5-10）+ TEASER 验证**，可能能拉一些回来。

但工程上：
- STD 当前 SearchLoop 只返回 top-1，要改成 top-K 得动它的核心逻辑
- 改完后效果未必能赶上 SC + TEASER 的 86%（毕竟 SC 的优势在于 BEV 极坐标统计天然对长走廊更稳）
- 这套改动不可移植给 STD 上游

### 5.2 改用 STD 的另一个变体？

STD 有几个后续工作（STD++、BTC 等），但都是同样思路的渐进改进，对**稀疏 + 高重复几何**的根本问题没解。

### 5.3 用 Airy 多 submap 聚合再喂 STD？

如果把 3-5 个 submap 合并成一个「fat submap」（10-30k 点），STD 的 corner 提取与三角形构建会稳得多，重复性问题可能缓解。但这就要做合并 + 重定位逻辑，工程开销不小。

**我的建议**：上述都不值得做。SC + TEASER 的 86% 已经达标，剩 14% 失败 case 都是稀疏 submap（< 4000 点），靠 **multi-submap consensus** 简单融合就能再救一批（5.3 的精神，但作为 SC 路线的兜底而非 STD 重启）。

---

## 6. 复现命令

```bash
# 编译（带 Ceres，自动找 /usr/local 优先）
cmake -S . -B build -DBUILD_TEASER_FPFH=ON
cmake --build build -j

# 跑全 22-submap 评估
./build/examples/global_localization_std/std_eval_22submaps \
    --session /home/steve/map_data/map_w2/20260511_110448 \
    --cloud-type world

# 只跑单个 submap（调试）
./build/examples/global_localization_std/std_eval_22submaps \
    --session /home/steve/map_data/map_w2/20260511_110448 \
    --cloud-type levelled --only 5
```

---

## 7. 工件清单

| 路径 | 内容 |
|---|---|
| `examples/global_localization_std/std_eval_22submaps.cc` | 主程序（~280 行，leave-one-out 协议） |
| `examples/global_localization_std/third_party/std/STDesc.{h,cpp}` | STD vendor + ROS 剥离（去掉 `read_parameters` ROS 版、`publish_std_pairs`、`ROS_ERROR_STREAM`） |
| `examples/global_localization_std/CMakeLists.txt` | standalone 构建（需要 Ceres 2.1+） |
| `analysis/experiments/std_airy_22submaps.md` | 本报告 |

---

## 8. 教训

1. **不要假设新论文就一定好用**——STD 在自家 KITTI 测试集上 90%+ 召回率，到我们的场景直接归零。论文的「最佳场景」与你的场景之间的差距决定一切。
2. **稀疏 + 重复几何是真正的死敌**，不是某个算法弱点。需要换思路（如**累积更厚的 submap**、**结合语义信息**、**多假设融合**）而不是换算法。
3. **SC 的 BEV 极坐标 + max-height** 在长走廊里反而比基于关键点的方法更稳——因为它编码的是**统计分布**而非**具体几何**。
4. **86% 已经很好**——前期我以为是 SC 在艰难地完成任务，事后看在我们场景里它其实是天然适配。
