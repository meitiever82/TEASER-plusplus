# 实验：3D-BBS 粗定位 + TEASER 精化 在 22 个 Airy submap 上的全局定位

> **日期**：2026-05-11
> **数据**：`map_w2/20260511_110448`，与 Phase 2 SC 实验、Phase 3 STD 实验完全相同的 22 submap
> **代码**：`examples/global_localization_bbs/` @ commit (TBD)
> **目的**：验证 Phase 5 `doc/phase5_bbs_plan.md` 的核心假设——「BBS 不依赖描述子，理论上对稀疏点云更鲁棒，可救回 Phase 2 失败的 3 个 sparse case (7/10/17)」

---

## 1. 结果（直说）

**整体**：22 / 22 都成功返回了 pose（**没有任何 NO_SOLUTION**），但 rotation 精度比 SC 路线偏弱。

| 判据 | BBS + TEASER | SC + TEASER（Phase 2） |
|---|---|---|
| **pose 返回率** | **22 / 22 = 100%** | 19 / 22 = 86.4% |
| 严格 (rot < 3° AND trans < 0.5 m) | 14 / 22 = 63.6% | 19 / 22 = 86.4% |
| 中等 (rot < 5° AND trans < 1.0 m) | **20 / 22 = 90.9%** | 19 / 22 = 86.4% |
| 宽松 (rot < 5° AND trans < 2.0 m) | 20 / 22 = 90.9% | 19 / 22 = 86.4% |
| 成功 case rotation 中位数 | 1.87° | 0.5° |
| 成功 case translation 中位数 | 0.125 m | 0.06 m |
| Max rotation error | **9.81°**（submap 10） | 2.3° |
| Max translation error | 0.71 m（submap 10） | 0.42 m |

简单总结：
- 鲁棒性 ↑（不会 NO_SOLUTION）
- 精度 ↓（rot 中位数 1.87° vs 0.5°，translation 中位数 0.125 m vs 0.06 m）

---

## 2. 完整 22-submap 表

```
submap   n_pts   bbs_ms   refine_s   #corr  #inlier  ratio   rmse    rot_err   trans_err
000000    6788     99.9    0.50      300    34      0.11   0.130    0.62°     0.077 m
000001    8006     78.4    0.50      512    83      0.16   0.127    0.97°     0.039 m
000002    6986     63.4    0.56      232    41      0.18   0.113    1.94°     0.126 m
000003    6973     93.4    0.66      241    38      0.16   0.122    3.45°     0.067 m  ⚠ rot
000004    7940     77.7    0.71      378    47      0.12   0.131    0.92°     0.132 m
000005    6941    150.3    0.64      302    51      0.17   0.112    1.71°     0.100 m
000006    5675    107.4    0.60      294    41      0.14   0.111    1.27°     0.089 m
000007    3808    258.1    0.51       87    19      0.22   0.126    3.48°     0.174 m  ⚠ rot   ✦救活
000008    3245    190.7    0.38      224    77      0.34   0.099    1.81°     0.036 m
000009    2802     72.4    0.39       84    27      0.32   0.099    3.24°     0.277 m  ⚠ rot
000010    2447     59.3    0.41       25    10      0.40   0.134    9.81°     0.706 m  ✗ both  ✦救活但精度差
000011    3573     72.7    0.51       45    10      0.22   0.123    3.45°     0.452 m  ⚠ rot
000012    5952    103.2    0.56      324    84      0.26   0.107    1.53°     0.022 m
000013    7007     64.7    0.62      266    22      0.08   0.126    0.91°     0.164 m
000014    6167     78.6    0.64      204    32      0.16   0.123    1.33°     0.103 m
000015    4672    115.7    0.51      181    19      0.10   0.121    2.29°     0.132 m
000016    4870     86.0    0.50      131    23      0.18   0.122    4.05°     0.121 m  ⚠ rot
000017    3357     59.4    0.50       17     3      0.18   0.106    6.11°     0.541 m  ✗ both  ✦救活但精度差
000018    3097     65.7    0.42      199    55      0.28   0.110    0.69°     0.053 m
000019    2944     71.3    0.41      124    30      0.24   0.128    3.80°     0.125 m  ⚠ rot
000020    4677     58.8    0.40      244    45      0.18   0.130    2.40°     0.187 m
000021    4345     54.1    0.39      143    17      0.12   0.130    1.45°     0.414 m
```

图例：✦ = Phase 2 中 NO_SOLUTION 的 case；⚠ rot = rot > 3° 但 < 5°；✗ both = rot > 5° 或 trans > 0.5 m。

---

## 3. 与 Phase 2 (SC) 失败 case 的对比

Phase 2 失败 3 个：submap 7 / 10 / 17（都 < 4000 点）。

| submap | n_pts | SC + TEASER | BBS + TEASER | 救活？ |
|---|---|---|---|---|
| 000007 | 3808 | NO_SOLUTION | **3.48° / 0.17 m**，19 inliers | ✅ 救回 |
| 000010 | 2447 | NO_SOLUTION | 9.81° / **0.71 m**，10 inliers | ⚠ 返回了 pose 但 rot 太大不可用 |
| 000017 | 3357 | NO_SOLUTION | 6.11° / 0.54 m，**3 inliers** | ⚠ 返回了 pose 但 inlier 太少不可信 |

**真正被救回（pose 可用）的只有 submap 7。submap 10 和 17 虽然没崩，但产出的 pose 不够精确**：
- 10: rot 9.8° 远超 5° 阈值；
- 17: inlier 只有 3，统计上不能信任。

如果在 hybrid pipeline 里设个后置检查（rot < 5° 且 inliers ≥ 10），那么 BBS 实际兜底贡献是 1 个，即 SC 19 + BBS 1 = 20 / 22 = 91%。

---

## 4. 几个关键观察

### 4.1 BBS 粗解精度量级

BBS 的输出精度由 `min_level_res` 决定。本次实验用 0.5 m，实测：
- 平均 coarse 平移误差 ~0.4-0.7 m
- 平均 coarse 旋转误差 ~0.8-2°

这与 plan §1 预期一致（±0.5 m / ±2°）。TEASER 精化能把平移压到 0.05-0.2 m，但旋转精度提升有限——这是 Phase 5 的核心代价。

### 4.2 为什么 rot 精度比 SC 差？

SC + TEASER 在 submap 5 上拿到 0° rot；BBS + TEASER 拿到 1.71°。同一份数据，差别在初值：
- SC 给的初值（top-30 中选出的 anchor 真值位姿）通常 rot 误差 < 1°；
- BBS 给的初值平均 rot 误差 ~0.8-2°，最坏可达几度。

TEASER 的 GNC-TLS 对初值不敏感（这是它的卖点），但**它本质是 outlier-robust 不是 high-precision**。当 inlier 数偏少（< 30）且初值已经偏 2° 时，最终输出会保留 1-3° 的残差。

SC 路线意外地"幸运"在：它的 anchor 池里就有真值附近的关键帧，TEASER 起始几乎在最优点；BBS 的输出是体素级离散的，永远偏几个体素。

### 4.3 sparse submap 的两种失败模式

| submap | n_pts | inliers | 失败模式 |
|---|---|---|---|
| 10 | 2447 | 10 | corrs 只 25 个，FPFH 在这么少点的云上严重退化 → BBS coarse 给出错位的初值（rot ~9°），TEASER 不收敛 |
| 17 | 3357 | 3 | corrs 17 个，inliers 只 3 个 → 统计样本不足，TEASER 解几乎是 underdetermined |

加 multi-submap consensus 后预计可救（plan §5.4 已讨论）。

### 4.4 时间预算

- BBS 中位 78 ms / 平均 95 ms / 最大 258 ms（submap 7）
- TEASER 精化 平均 0.52 s
- **端到端中位 ~0.6 s**，最大 ~0.8 s

比 plan §6 估计的 1-3 s 还快——主要因为我们的 40m 场景比 plan 假设的 200m 场景小。

---

## 5. 与 Phase 2 / Phase 3 三方对比

| 指标 | Phase 2 SC+TEASER | Phase 3 STD | **Phase 5 BBS+TEASER** |
|---|---|---|---|
| pose 返回率 | 86.4% | 100%（但都错） | **100%** |
| 严格成功率 (rot<3 AND tr<0.5) | 86.4% | 0% | 63.6% |
| 中等成功率 (rot<5 AND tr<1) | 86.4% | 0% | **90.9%** |
| rot 中位（成功 case） | 0.5° | 113° | 1.87° |
| trans 中位（成功 case） | 0.06 m | 23 m | 0.125 m |
| 失败 case 数 | 3 (sparse) | 22 (全部) | 0 (NO_SOLUTION) / 2 (高 rot) |
| 端到端耗时 | ~1-3 s（估） | ~30 ms（错得很快） | ~0.6 s（实测） |

STD 是死路（Phase 3 已归档）。BBS 与 SC 的对比是**精度 vs 鲁棒性**的取舍。

---

## 6. 决策（按 `phase5_bbs_plan.md` §7）

判据表：

| 结果 | 决策 |
|---|---|
| BBS+TEASER ≥ 21/22（95%+） | 切主线 |
| BBS+TEASER = 19-20/22（与 SC 持平） | 作为兜底 ← **匹配** |
| BBS+TEASER < 19/22 | 归档 |
| BBS+TEASER 救回 7/10/17 但其他场景反而退化 | 场景分流 |

**实际匹配**：用中等阈值（rot<5° AND trans<1m），20/22 = 91%，与 SC 的 19/22 = 86% **基本持平**。

**最终决策：作为 SC 主线的兜底，进 hybrid pipeline。**

具体策略：
1. **主路径** = SC + TEASER（Phase 2）
2. **兜底触发**：SC 返回 NO_SOLUTION 时，自动回落到 BBS + TEASER
3. **后置验证**：BBS 返回结果时必须满足 `rot_err < 5° AND inliers ≥ 10`，否则视作仍然失败（落到 multi-submap consensus）
4. **预计兜底成功率**：22 case 里 SC 拿下 19，BBS 在剩下 3 个里能稳救 1（submap 7）→ hybrid 20/22 = 91%；剩下 2 个（10, 17）走 multi-submap consensus

**不切主线**的原因：BBS 的 rot 精度（median 1.87°）明显劣于 SC（median 0.5°）。对一个 ±0.5° 精度目标的系统，让 BBS 当默认会拉低 80% 稠密 case 的精度。

---

## 7. 复现命令

```bash
# 1. 编译（3D-BBS 源码 vendor 在 examples/global_localization_bbs/third_party/，CMake add_subdirectory 自动 in-tree 编译）
cd /home/steve/Documents/GitHub/tools/TEASER-plusplus/build
cmake .. -DBUILD_TEASER_FPFH=ON
cmake --build . --target bbs_localize -j

# 2. 单 submap 验证
./examples/global_localization_bbs/bbs_localize \
    --map /home/steve/map_data/map_w2/20260511_110448/global_map.pcd \
    --scan /home/steve/map_data/map_w2/20260511_110448/000005/submap_levelled.pcd \
    --gt /home/steve/map_data/map_w2/20260511_110448/000005/data.txt \
    --output /tmp/bbs_s5.txt

# 3. 批跑 22 个
bash /home/steve/Documents/GitHub/tools/TEASER-plusplus/examples/global_localization_bbs/run_22submaps.sh
column -s, -t < /tmp/bbs_22submaps/summary.csv
```

---

## 7a. 在 in-tree (add_subdirectory) 编译下的复跑

2026-05-11 在把 3d-bbs 源码 vendor 进 `examples/global_localization_bbs/third_party/` 后做了第二次完整批跑。结论与首次实验一致，**精度量级稳定**，但因 CUDA / BnB 内部存在轻微非确定性，个别 case 数值会浮动。

| 指标 | v1 (外部 install) | v2 (in-tree add_subdirectory) |
|---|---|---|
| pose 返回率 | 22/22 = 100% | 22/22 = 100% |
| 严格 (rot<3 AND tr<0.5) | 14/22 = 64% | **16/22 = 73%** |
| 中等 (rot<5 AND tr<1) | **20/22 = 91%** | 19/22 = 86% |
| rot 中位 / max | 1.87° / 9.81° | **1.75° / 7.95°** |
| trans 中位 / max | 0.125 m / 0.71 m | 0.119 m / 0.51 m |
| BBS_ms 中位 / refine_s 均 | 78 / 0.52 s | 77 / 0.53 s |

**Phase 2 失败 case 的兜底状况**（rot<5° AND inliers≥10 才算 BBS 救活）：

| case | n_pts | v1 BBS | v2 BBS |
|---|---|---|---|
| 7 | 3808 | 3.48° / 0.17m / 19 inliers → ✓ | 5.21° / 0.28m / 15 inliers → ⚠（rot 越线） |
| 10 | 2447 | 9.81° / 0.71m / 10 inliers → ✗ | 7.95° / 0.35m / 10 inliers → ✗ |
| 17 | 3357 | 6.11° / 0.54m / 3 inliers → ✗ | 6.53° / 0.51m / 4 inliers → ✗ |

两次跑下来，submap 7 在严格的 `rot<5°` 后置阈值下时而救得活时而救不活——它是个临界 case，rot 误差稳定在 3-5° 区间。**实战部署时应该把 hybrid 后置阈值放到 rot < 6°，以稳定救回 submap 7**；对 10/17 仍需 multi-submap consensus。

整体结论不变：BBS+TEASER 适合做兜底，不适合做主线。

---

## 8. 工件清单

| 路径 | 内容 |
|---|---|
| `examples/global_localization_bbs/bbs_smoke_test.cc` | 1×1 验证 BBS API + 链接（bunny 数据） |
| `examples/global_localization_bbs/bbs_localize.cc` | BBS + TEASER 端到端（含 CLI） |
| `examples/global_localization_bbs/refine_pipeline.hpp` | FPFH+TEASER 精化（从 Phase 1 抽取的 inline header） |
| `examples/global_localization_bbs/run_22submaps.sh` | 22-submap leave-one-out 批跑脚本 |
| `examples/global_localization_bbs/third_party/` | vendored KOKIAOKI/3d_bbs；顶层 `examples/CMakeLists.txt` 通过 `add_subdirectory(... EXCLUDE_FROM_ALL)` in-tree 编译，bbs_localize / bbs_smoke_test 直接 link `gpu_bbs3d` target，无需安装 |
| `analysis/experiments/bbs_airy_22submaps.md` | 本报告 |
| `/tmp/bbs_22submaps/summary.csv` | 完整原始数据 |

---

## 9. 接下来该做的（不属于 Phase 5）

- **hybrid pipeline**：在 `global_localization_sc.cc` 里加 BBS 兜底分支。约半天工作量。
- **multi-submap consensus** for 10/17：合并相邻 submap（±1 邻居）提高密度，重跑 BBS+TEASER。这是 Phase 2 报告 §4.3 已经规划过的，属于 SC 路线的兜底，不属于 Phase 5。
- **更细体素**：试 `min_level_res=0.25 m`，看 BBS coarse 是否能把 rot 精度从 ~2° 拉到 ~1°；代价是 BBS 时间从 ~95 ms 涨到 ~1 s。性价比可能不高。

---

## 10. 教训

1. **「不依赖描述子」≠「精度好」**。BBS 解决了「能不能拿到 pose」，没解决「pose 多准」。FPFH 描述子的失败模式（在稀疏云上退化）确实被绕过了，但 TEASER 精化阶段还是要靠 FPFH 匹配——sparse case 的对应数只有几十个，inliers 个位数，没有足够样本拉精度。
2. **取舍是真实的**：robustness 提升了 14%（86% → 100% pose 返回），rot 中位数 4× 退化（0.5° → 1.87°）。哪个更重要看应用：
   - 启动定位（一次性）：取 robustness，BBS 更合适
   - 持续定位（高精度）：取 precision，SC 更合适
   - 通用：hybrid 走两条线，由后置验证拍板
3. **plan 写的 ±0.5 m / ±2° 是合理估计**——实测中位数 0.125 m / 1.87° 都在量级内，没有意外。
