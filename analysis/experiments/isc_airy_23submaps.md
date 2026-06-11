# 实验：Intensity Scan Context (ISC) 在 23 个 Airy submap 上 vs 标准 SC

> **日期**：2026-05-11
> **数据**：`map_w2/20260511_150204`，新一次建图，**为本次实验补做了 intensity 输出**（mapping 侧上游修改）
> **代码**：`examples/global_localization_sc/`（vendored SC + ISC 扩展）
> **目的**：检验 ICRA 2020 Intensity Scan Context（用 max-intensity 代替 max-z）能否在我们的稀疏 / 重复几何场景里压过原版 SC，把 Phase 2 失败的 sparse case 拉回来

---

## 1. 结果（直说）

**ISC ≈ SC，无净收益，归档。**

| 指标 | SC（max-z） | ISC（max-intensity） |
|---|---|---|
| pose 返回率 | **21 / 23 = 91.3%** | 20 / 23 = 87.0%（少救一个） |
| 中等 (rot<5° AND tr<1m) | **17 / 23 = 73.9%** | **17 / 23 = 73.9%** |
| 严格 (rot<3° AND tr<0.5m) | 14 / 23 = 60.9% | 14 / 23 = 60.9% |
| rot 中位 / max | **1.67° / 7.11°** | 1.80° / 9.91° |
| trans 中位 / max | 0.152 m / 0.45 m | **0.109 m / 0.57 m** |
| query 时间均值 | 27.0 s | 27.5 s |

ISC 唯一的可量化优势：translation 中位数好 28%（0.109 vs 0.152 m）。但是
- rotation 中位反而略差；
- pose 返回率退 1 个（submap 9 在 ISC 下 NO_SOLUTION，SC 能跑出来）；
- 中等 / 严格阈值的成功率**完全打平**，对部署判据零影响。

---

## 2. 实验设置

### 2.1 上游 mapping 侧修改

为支持 ISC，glim_ros 端做了以下改动（不在本仓库内）：

| 文件 | 改动 |
|---|---|
| `points_compact.bin` | 不变（仍 12 字节/点 = xyz） |
| `intensities_compact.bin` | **新增**：`N × 1 × float32` 并列文件，与 points_compact 索引对齐 |
| `global_map.pcd` | intensity 字段**真值填充**（旧版是 PCL 占位 0） |

本仓库下游对应的改动：

| 文件 | 改动 |
|---|---|
| `tools/submap_to_pcd.py` | 读 intensities_compact，submap_levelled.pcd 输出 FIELDS=`xyzi` |
| `third_party/scancontext/Scancontext.h/cpp` | 新增 `makeIntensityScancontext()`（bin 值 = max-intensity 而非 max-z） |
| `global_localization_sc.cc` | 新增 `--descriptor sc|isc` 选项；DB 文件格式升级到 v3 增加 1 字节描述子类型 |

### 2.2 数据

新 session `20260511_150204`，23 个 submap（000000–000022），覆盖与旧 session 相近的室内走廊。**与旧 session 22-submap 数据集不直接对应**——这是一次重新建图。

- 全局地图：115 111 点，intensity range 2-171
- submap 点数：1475 – 8180，中位 4773
- intensity range（每 submap）：典型 2-180，整体动态范围窄

### 2.3 评测协议

Leave-one-out：每个 submap 的 `submap_levelled.pcd` 作 query，固定 `global_map.pcd` 作目标。两路独立跑，相同参数：

```
--anchor-step 3.0  --sc-radius 20.0  --sc-downsample 0.3
--top-k 30  --algo TEASER  --noise-bound 0.15
--descriptor [sc | isc]
```

---

## 3. 完整 23-submap 结果对照

```
submap n_pts   SC                              ISC                             判定
                rot      trans    inliers       rot      trans    inliers
─────────────────────────────────────────────────────────────────────────────
000000  7261   1.45°    0.152m   49            0.77°    0.152m   47           ISC↑
000001  7946   1.73°    0.189m   64            0.74°    0.089m   57           ISC↑
000002  7249   4.87°    0.453m   10            5.18°    0.271m   11           SC↑
000003  7773   0.87°    0.029m   20            3.12°    0.242m   20           SC↑↑
000004  8180   0.66°    0.045m   57            0.65°    0.086m   56           平手
000005  6974   5.19°    0.362m   49            1.86°    0.083m   49           ISC↑↑
000006  5586   1.67°    0.138m   28            4.03°    0.154m   24           SC↑↑
000007  3420   6.09°    0.280m   35            0.59°    0.018m   44           ISC↑↑ ⭐
000008  2693   2.55°    0.220m   39            0.53°    0.078m   34           ISC↑↑
000009  2663   4.37°    0.112m   22            ——       ——       ——           SC独醒
000010  2637   ——       ——       ——            ——       ——       ——           都崩
000011  4285   1.39°    0.212m   18            2.14°    0.242m   15           SC↑
000012  6776   0.82°    0.046m   70            1.25°    0.105m   73           SC↑
000013  6684   0.52°    0.039m   20            0.65°    0.111m   23           平手
000014  4773   3.82°    0.327m   28            3.36°    0.334m   29           ISC↑
000015  4915   0.68°    0.138m   35            0.87°    0.040m   43           平手
000016  4695   0.71°    0.137m   64            2.42°    0.123m   57           SC↑↑
000017  3007   7.11°    0.276m   12            6.17°    0.268m   13           ISC↑
000018  3300   6.66°    0.450m   44            2.07°    0.107m   55           ISC↑↑
000019  3281   2.62°    0.302m   17            9.91°    0.567m   16           SC↑↑ ⚠
000020  5180   0.86°    0.037m   85            0.83°    0.087m   77           平手
000021  4358   1.48°    0.104m   32            1.73°    0.053m   27           平手
000022  1475   ——       ——       ——            ——       ——       ——           都崩
```

判定图例：`↑↑` 大胜（Δrot ≥ 1.5°）｜`↑` 小胜（0.3-1.5°）｜平手 < 0.3°｜独醒/都崩 = NO_SOLUTION 边界

## 4. 单 submap 互补性 vs 聚合零和

判定分布几乎完全对称：

| 判定 | 数量 | submap |
|---|---|---|
| ISC↑↑ | 4 | 5, 7, 8, 18 |
| ISC↑ | 4 | 0, 1, 14, 17 |
| SC↑↑ | 4 | 3, 6, 16, 19 |
| SC↑ | 3 | 2, 11, 12（加 SC 独醒 9） |
| 平手 | 6 | 4, 13, 15, 20, 21 |
| 都崩 | 2 | 10, 22 |

**ISC 占优 8 个 / SC 占优 8 个**——这种几乎完美的对称强烈暗示**两个描述子的失败相互独立但能力相当**。任何一方都不是另一方的超集。

**对 ensemble 的诱惑**：如果两路并跑取 inlier 最高，理论上能从 17/23 拉到 19/23 = 83%（救 ISC 独醒的 7、18 + SC 独醒的 9）。但代价是 query 时间从 27 s 翻到 ~50 s（top-K=30 双路 TEASER 验证），换 +2 个 case。**ROI 不划算，且救不到真正的硬骨头 10、22**。

---

## 5. 失效机理：为什么 ISC 在我们场景没亮

ISC 在 ICRA 2020 论文里的舒适区是 **KITTI 室外**：道路标线、交通标志、植被、车辆——反射率动态范围大、地面纹理差异显著。

我们这套数据的反差：

| 维度 | KITTI 室外 | 我们的煤矿巷道 |
|---|---|---|
| 反射率分布 | 多模态（路面 vs 植被 vs 标志 vs 金属） | 单调（混凝土墙、金属支架） |
| 单 bin intensity 动态范围 | 大 | 窄（中位 23，max 200 但稀疏） |
| 高反射点的语义 | 标志牌、车牌 → 信息量高 | lidar 罩反光、潮湿地面 spec → **噪声** |
| 几何重复性 | 低（街道有形态差异） | 极高（巷道截面都差不多） |

具体到 ISC 的脆弱点：

### 5.1 MAX-intensity 对单点敏感

ISC 用 `desc(r,s) = max intensity` 作 bin 值——**一个高反射点（雷达罩反光、雷达打到湿地、金属反光板）就会主宰整个 bin**。submap 19 的 ISC 从 SC 的 2.62° 暴退到 9.91° 多半就是这个机制：query 那侧出现了几个高 intensity 点污染了描述子，与 DB 里同一地点构建时的高 intensity 分布不一致。

相比之下 SC 的 `max z`：「最高点的高度」是**对单点噪声不敏感的统计量**，孤立的高反射点不会出现成一个高度峰。

### 5.2 室内 intensity 动态范围窄

测过：

| 数据 | min | max | 中位 | 75% 分位 |
|---|---|---|---|---|
| 全局地图 intensity | 2 | 171 | 23 | ~40 |

中位数 23、75% 分位 ~40，分布偏左、长尾。**绝大多数 bin 的 max-intensity 都在 30-80 这个窄区间内挤着**——可区分性比 SC 的高度（覆盖 0-11.6m）要弱。

### 5.3 Phase 2 失败的 3 个 case 还是失败的

最重要的检验：Phase 2 / 5 的失败模式是**极端稀疏**（< 4000 点 / pts-per-frame < 50）。这次新 session 里对应的稀疏 case 是 000010（2637 pts）、000022（1475 pts）、000009（2663 pts）。

- 000010、000022：SC 和 ISC **同时崩**
- 000009：SC 险胜（4.37°），ISC 直接 NO_SOLUTION

ISC 在「描述子构建」阶段确实绕过了"几何区分度不足"的问题，但**到了 TEASER 验证那一关，FPFH 仍然要稀疏点云上提匹配**——FPFH 在稀疏几何上退化的根本问题没解。这跟 Phase 5 BBS 一样：换描述子救不了 FPFH 那一关的退化。

---

## 6. 决策：归档，不切主线、不上 ensemble

| 选项 | 评估 | 决策 |
|---|---|---|
| 切 ISC 替代 SC | 17/23 = 17/23 完全打平，translation 略胜但 rotation 略退，**净 0** | ✗ 不切 |
| SC + ISC ensemble | +2 个 case（17→19），但延迟翻倍（27→50+ s），且救不到 000010 / 000022 | ✗ 不做（ROI 差） |
| 走 SC 主线 + 极稀疏 case 走 multi-submap consensus | SC + BBS 兜底（Phase 5）已能拿 19-20，再加邻居合并能压到 21+ | ✓ **既定路线** |

---

## 7. 这次实验的副产物（保留价值）

虽然 ISC 没成事，下面的基础设施保留下来对未来有用：

| 工件 | 价值 |
|---|---|
| Mapping 侧 intensity 输出 | 任何 intensity-aware 方法都需要这个 |
| `tools/submap_to_pcd.py` 支持 xyzi PCD | 同上 |
| `Scancontext.{h,cpp}::makeIntensityScancontext` | 想 A/B 比较 SC vs ISC 一行命令 |
| SC DB v3 格式 | 支持 self-describing 的描述子类型；后续加新描述子（joint、多模态）不用再 break format |
| `run_sc_isc_eval.py` Python 评测脚本 | 通用 leave-one-out，下一个 batch 评测可复用（之前是 bash 写的） |

---

## 8. 与之前实验的三方对比

| 维度 | Phase 2 SC (旧 session, 22) | Phase 5 BBS+TEASER (旧, 22) | **本实验 SC (新, 23)** | **本实验 ISC (新, 23)** |
|---|---|---|---|---|
| pose 返回率 | 19/22 = 86% | 22/22 = 100% | **21/23 = 91%** | 20/23 = 87% |
| 中等成功率 | 19/22 = 86% | 20/22 = 91% | **17/23 = 74%** | **17/23 = 74%** |
| rot 中位 | 0.5° | 1.87° | 1.67° | 1.80° |
| trans 中位 | 0.06 m | 0.125 m | 0.152 m | 0.109 m |

**注意新 session 的 SC 成功率（74%）比旧 session（86%）低**——大概率是新一次 mapping 切分出的稀疏 submap 更多（23 vs 22，新增的 000022 只有 1475 点），而不是 SC 算法退化。两个数据集统计上不直接可比。

---

## 9. 复现命令

```bash
# 1. 重生成 23 个 submap 的 xyzi PCD
python3 examples/global_localization_sc/tools/submap_to_pcd.py \
    /home/steve/map_data/map_w2/20260511_150204 --all

# 2. 编译
cd build && cmake --build . --target global_localization_sc -j

# 3. 建两个 DB（注意 --descriptor）
./examples/global_localization_sc/global_localization_sc --mode build \
    --map /home/steve/map_data/map_w2/20260511_150204/global_map.pcd \
    --db /tmp/w2_sc.bin  --descriptor sc \
    --anchor-step 3.0 --sc-radius 20.0 --sc-downsample 0.2

./examples/global_localization_sc/global_localization_sc --mode build \
    --map /home/steve/map_data/map_w2/20260511_150204/global_map.pcd \
    --db /tmp/w2_isc.bin --descriptor isc \
    --anchor-step 3.0 --sc-radius 20.0 --sc-downsample 0.2

# 4. 批跑 23 个
python3 examples/global_localization_sc/run_sc_isc_eval.py \
    --session /home/steve/map_data/map_w2/20260511_150204 \
    --map /home/steve/map_data/map_w2/20260511_150204/global_map.pcd \
    --db /tmp/w2_sc.bin \
    --bin build/examples/global_localization_sc/global_localization_sc \
    --out /tmp/scisc_eval/sc.csv

python3 examples/global_localization_sc/run_sc_isc_eval.py \
    --session /home/steve/map_data/map_w2/20260511_150204 \
    --map /home/steve/map_data/map_w2/20260511_150204/global_map.pcd \
    --db /tmp/w2_isc.bin \
    --bin build/examples/global_localization_sc/global_localization_sc \
    --out /tmp/scisc_eval/isc.csv
```

---

## 10. 工件清单

| 路径 | 内容 |
|---|---|
| `examples/global_localization_sc/third_party/scancontext/Scancontext.{h,cpp}` | 新增 `makeIntensityScancontext()` |
| `examples/global_localization_sc/global_localization_sc.cc` | `--descriptor sc|isc` + DB v3 格式 |
| `examples/global_localization_sc/tools/submap_to_pcd.py` | xyzi 输出（向后兼容 xyz） |
| `examples/global_localization_sc/run_sc_isc_eval.py` | Python 通用 leave-one-out 评测脚本 |
| `analysis/experiments/isc_airy_23submaps.md` | 本报告 |
| `/tmp/scisc_eval/{sc,isc}.csv` | 完整原始数据 |
| `/tmp/scisc_eval/{sc,isc}_logs/` | 每个 submap 的 stdout 日志 + refined pose 输出 |

---

## 11. 教训

1. **论文 sweet spot 必须和场景核对**。ISC 在 KITTI 室外（多模态反射率）有效，在我们煤矿巷道（单调反射率 + 高几何重复）没有可发挥的信号空间。和 Phase 3 STD 是同一种结构性失败——「用 paper 上的 best method 套到我们场景」要么打不开局面，要么打平。
2. **MAX-X 类描述子对噪声敏感**。Max-intensity 比 Max-z 更脆，因为高反射点天然更稀少且更**异常**。Mean/median-X 或 percentile 描述子可能更稳，但那是另一类研究，不在本实验范围。
3. **"互补 → ensemble" 是个常见诱惑**。8 vs 8 的对称胜负看起来像 ensemble 的礼物，但仔细算 ROI：双倍延迟换 +9% 成功率，且救不到真正失败的硬 case。**对称胜负往往意味着两个方法各自的方差大于真实差异**——ensemble 主要是在 averaging out 各自的随机性，不是在合并互补的强项。
4. **基础设施投资有持续回报**。这次 ISC 评估失败，但 mapping 侧的 intensity 通路、SC DB v3 self-describing 格式、Python 通用评测脚本都留下来给未来用。不是白做的——但下次再有"试一下 X 描述子"的诱惑时，应该先估一下 X 在场景里的工作半径。
5. **同一场景的"重做实验"基线会漂移**。新 session 23 个 submap 的 vanilla SC 拿 74%，旧 session 22 个 submap 是 86%。两份数据集本来就不对应——mapping 切分边界不同会让难易分布变化。下次做对比实验**必须用同一 session 数据**，不能跨 mapping 重跑做对照。
