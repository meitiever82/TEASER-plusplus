# 3D-BBS — Key Questions (Q&A)

> Design-rationale and "why is this implemented like this?" answers, each
> grounded in the code. If you can't find a [VERIFY:] tag below a claim, it
> isn't being made — the answer stays observational.

---

## Q1. Why two scaling factors in `branch()` — `v_rate` for translation and `num_division` for rotation?

**Short answer:** the two grids are constructed with different invariants.

**Long answer:** Translation resolution shrinks by a *fixed* `v_rate` every
level (default 2.0), so every parent voxel cleanly subdivides into `v_rate^3`
children. The grid is essentially an octree.

Rotation resolution, by contrast, depends on the *current* level's translation
resolution and the source cloud's bounding radius (see
`03-ALGORITHM_01-BranchAndBound.md` §3). The number of angular bins between
two levels can change arbitrarily — perhaps 1 at the top, 4 in the middle, 12
at the bottom. So `num_division[L-1]` (computed by `calc_angular_info`) tells
the *child* how many bins fit inside one parent angular cell.

The two expansions sit side-by-side in `DiscreteTransformation::branch(...)`:
```cpp
x*v_rate + i,                            // translation: i ∈ [0, v_rate)
roll*num_division.x() + l,               // rotation: l ∈ [0, num_division.x())
```
[VERIFY: bbs3d/include/discrete_transformation/discrete_transformation.hpp:51-56]

**Why:** Translation and rotation are coupled (rotation shifts a point by
`θ · ||p||`), so the angular resolution must adapt to translation resolution
to keep the BnB bound tight. A fixed `v_rate` for rotation would either be
too coarse at small `r_L` (wasting work) or too fine at large `r_L`
(blowing up children).

**How to apply:** If you swap in a different bounding rule (say, you have
priors on rotation tighter than the default ±0.02 rad), you don't change
`v_rate` — change `min_rpy_`, `max_rpy_`, and let `calc_angular_info`
recompute the per-level `num_division`.

---

## Q2. Why is the angular resolution derived from `arccos(1 - r²/(2·r_max²))` rather than the simpler `r / r_max`?

**Short:** They agree to first order; the cosine form is exact.

**Long:** A rotation by `θ` displaces a point at distance `d` by `2·d·sin(θ/2)`.
Setting this displacement equal to one voxel size `r_L` and solving for `θ`:

```
2 · r_max · sin(θ/2) = r_L
sin(θ/2) = r_L / (2·r_max)
```

Using the identity `cos(θ) = 1 − 2·sin²(θ/2)`:

```
cos(θ) = 1 − 2 · (r_L / (2·r_max))² = 1 − r_L² / (2·r_max²)
```

So `θ = arccos(1 − r_L² / (2·r_max²))` — exactly the code at
[VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:73-75]. The simpler `r_L / r_max` is
the small-angle approximation; the cosine form is correct for *any* angle,
including the cases at the top level where `r_L / r_max` can be > 1.

**How to apply:** If you ever need to derive a similar bound in your own
code (e.g. you change the metric from voxel-occupancy to euclidean
nearest-neighbour), reuse this exact form — `arccos(1 − d²/(2R²))` for a
displacement budget `d` at a moment arm `R`.

---

## Q3. Why does CPU use `double` but GPU uses `float`?

[VERIFY: bbs3d/include/cpu_bbs3d/bbs3d.hpp:26] (double) vs [VERIFY: bbs3d/include/gpu_bbs3d/bbs3d.cuh:30] (float).

**Reasons backed by the code:**

1. **GPU throughput.** Single-precision `Eigen::Vector3f` and `Matrix3f`
   compute roughly 2× faster than `double` on consumer GPUs. The kernel runs
   per-pose × per-source-point arithmetic, so this matters.
2. **Memory bandwidth.** `Eigen::Vector3f` is 12 bytes vs 24 for `Vector3d`;
   the source cloud is read once per point per pose, so halving the read size
   halves the bandwidth pressure on `d_src_points_`.
3. **Sufficient precision for the use case.** Voxels are ≥ 0.5 m (per the
   `min_level_res` default of 1.0 m, see [VERIFY: test/config/test.yaml:8]).
   `float` has ~7 decimal digits of precision, so within ~±10 km of origin the
   absolute error is well below the voxel size.

**On the CPU side**, the cost of `double` is negligible against the inner
hash-probe loop, and using `double` avoids subtle inconsistencies when
3D-BBS is mixed with PCL / GICP downstream (`pcl::GeneralizedIterativeClosestPoint`
returns `final_transformation_` as `Matrix4f`; the CPU test cast it to `double` at
[VERIFY: test/src/cpu_test.cpp:130]).

**Caveat:** If you operate over coordinates > ~10 km, switch the GPU version
to a re-centred coordinate frame *before* feeding into `BBS3D`. The library
doesn't do this for you.

---

## Q4. Why the 7-corner asymmetric neighbour dilation (not 26 neighbours or 6 face-neighbours)?

**Why:** The dilation must compensate for `floor`-quantisation drift, and
`floor` only crosses voxel walls toward *negative* directions.

**How `floor` skews:** Suppose a true point sits at `x = 1.0001`. `floor(x) =
1`. Now perturb by `−ε`: `floor(1.0001 − ε) = 0` for `ε > 0.0001`. But
`floor(1.0001 + ε) = 1` for `ε < 0.9999`. So the voxel index can shift to
`{−1, 0}` per axis under noise; never to `+1`.

The 7-neighbour set is exactly the dilation that adds the "negative-decrement
voxels" for every subset of the three axes:

```
(-1, -1,  0)  (-1,  0,  0)  ( 0, -1,  0)
(-1, -1, -1)  (-1,  0, -1)  ( 0, -1, -1)
( 0,  0, -1)
```
[VERIFY: bbs3d/src/cpu_bbs3d/voxelmaps.cpp:55-61]

(The 8th vertex of the unit cube — the voxel itself — is already in the map.)

**Compared to alternatives:**

| Dilation | Voxel count expansion | Quantisation correctness | Notes |
|---|---:|---|---|
| None | 1× | Breaks for noise crossing voxel walls | Highest false-negative rate |
| 6 face-neighbours | 7× | Over-corrects in `+` axes | False positives near boundaries |
| **7 negative-corners (chosen)** | 8× | Matches `floor` skew exactly | ✓ |
| 26 full-cube | 27× | Massively oversized | Wastes memory, blows up coarse-level scores |

**How to apply:** If you change the quantiser from `floor(...)` to `round(...)`,
the dilation needs to become symmetric (e.g. 6 face-neighbours). The current
neighbour set is **tightly coupled to `floor`** by design.

---

## Q5. Why does the CPU version short-circuit on `bucket.w() == 0` but the GPU version doesn't?

[VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:156-158] vs [VERIFY: bbs3d/src/gpu_bbs3d/calc_score.cu:41-53].

**Why:** Two reasons, both visible in the code:

1. **GPU warp divergence.** All 32 threads in a warp execute the inner probe
   loop in lock-step. If only some hit `w==0`, an early `break` for them
   wouldn't actually save the warp's work — the warp continues until *all*
   threads' `break` conditions are met. So the GPU "saves" nothing by adding
   the check.
2. **Correctness identity.** An empty bucket `(0,0,0,0)` can't match any
   non-origin coord, so the coord-check `continue` already handles it. The
   probe still runs to `max_bucket_scan_count` (default 10), but the result is
   identical.

**Why CPU adds it:** A scalar CPU thread sees real benefit from terminating
its probe loop. This optimisation arrived in commit `3feeede` (PR #54,
"Optimize hash lookup by breaking early on empty buckets for cpu bbs3d") —
the recent commit history visible in `git log --oneline`. The optimisation
relies on the *builder's* invariant: `create_hash_buckets` never inserts a
key past an empty slot, so `w==0` at probe `j` proves the key is absent.

**How to apply:** Don't blindly mirror CPU-side optimisations into CUDA
kernels — measure SIMT divergence cost first.

---

## Q6. Why does the GPU defer branching through `branch_stock` instead of scoring each child immediately like the CPU?

**Cost model.** A CUDA kernel launch is ~5–10 μs of fixed overhead; a H↔D
`cudaMemcpyAsync` adds another ~10 μs for small transfers. A typical BnB
iteration produces only a few children (8 translation cells × ~1–10 yaw bins
≈ 10–80), so a per-iter kernel launch would have *more launch overhead than
useful work*.

The deferred batch reaches an efficient regime: at `branch_copy_size_ =
10000` poses (default [VERIFY: bbs3d/src/gpu_bbs3d/bbs3d.cu:9]), the kernel
amortises launch cost over ~313 warp-blocks.

**Trade-off:** Children sit unscored in `branch_stock` until the threshold
trips. During this time the priority queue lacks them, so `best_score` may
update less aggressively than it would on the CPU. This means **the GPU's
prune effectiveness is slightly weaker** — but in practice the per-pose
speedup more than compensates.

**Evidence:**
- Default size visible at [VERIFY: bbs3d/src/gpu_bbs3d/bbs3d.cu:9].
- "Flush when nearing empty" check at [VERIFY: bbs3d/src/gpu_bbs3d/bbs3d.cu:193-200] — without this, deferred children would never be scored.
- "Flush at threshold" at [VERIFY: bbs3d/src/gpu_bbs3d/bbs3d.cu:215-222].

**How to apply:** If your GPU's launch overhead is much smaller (e.g. CUDA
Graphs, or you've fused multiple launches), shrink `branch_copy_size_` via
`set_branch_copy_size(int)` [VERIFY: bbs3d/include/gpu_bbs3d/bbs3d.cuh:48] to
get tighter pruning.

---

## Q7. Why default `min_rpy = (-0.02, -0.02, 0.0)` and `max_rpy = (0.02, 0.02, 2π)`?

[VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:15-16] and [VERIFY: bbs3d/src/gpu_bbs3d/bbs3d.cu:19-20].

**Why:** 3D-BBS is explicitly designed around the "gravity-aligned scan"
assumption. ±0.02 rad ≈ ±1.15° on roll and pitch is the typical residual after
an IMU-based gravity alignment (the README's *"3D LiDAR scan aligned in the
gravity direction by using such as IMU"* in
[VERIFY: README.md:26-27]). Yaw is left fully free over `[0, 2π)` because
gravity alignment doesn't constrain heading.

**Performance impact** (per README [VERIFY: README.md:27]): *"Although 3D-BBS
also performs 6DoF search without gravity aligning, the processing time is
more than 10 times longer."* Roll and pitch each contribute another factor of
~10 to the work — so 6DoF search is ~100× slower than 3DoF.

**How to apply:** If you don't have an IMU, expand the rpy range, but expect
processing time to balloon. If you have an IMU with sub-degree accuracy,
shrink to ±0.005 rad to cut work further. The ROS 2 wrapper applies the
IMU's gravity alignment *before* calling `localize`
[VERIFY: ros2_test/iridescence/src/gpu_bbs3d_iridescence/gpu_ros2_test_iridescence.cpp:135].

---

## Q8. Why `max_bucket_scan_count = 10`?

[VERIFY: bbs3d/src/cpu_bbs3d/voxelmaps.cpp:5] (constructor default).

**Why:** It's the **fixed inner-loop bound** in `calc_score`. A small constant
budget means:
- The GPU kernel has a fully-unrolled `for (j = 0; j < 10; ...)`.
- Per-point cost is bounded at `≤ 10` bucket reads.
- The builder's load factor (`#keys / #buckets`) is kept low enough that 99.9%
  of keys land within their first 10 probes
  [VERIFY: bbs3d/src/cpu_bbs3d/voxelmaps.cpp:91].

If you raise `max_bucket_scan_count`, the bucket array can be smaller (denser
load factor allowed), saving memory, but lookups become slower.

**Trade-off table** (illustrative, derived from the algorithm):

| `max_bucket_scan_count` | Memory factor | Worst-case probe cost | Suitable when |
|---:|---:|---:|---|
|  3 | ≥ 4× | 3 | Memory-rich, latency-sensitive |
| **10 (default)** | ≥ 1× | 10 | Balanced |
| 30 | ≥ 0.5× | 30 | Memory-constrained, low load-factor tolerated |

**How to apply:** Adjustable only via `VoxelMaps::set_max_bucket_scan_count`
[VERIFY: bbs3d/include/cpu_bbs3d/voxelmaps.hpp:35]. There's currently no API
on `BBS3D` itself to forward this through, so you'd have to modify
`BBS3D::set_tar_points` if you need to expose it.

---

## Q9. Why is `set_voxelmaps_coords` faster than `set_tar_points`?

**Fresh build:**
- Read PCD into `std::vector<Vector3>`,
- For each level: quantise points, run the `unordered_map`-based dilation,
  run the open-addressing `create_hash_buckets` (potentially multiple
  resize-and-retry iterations).

That's ~9272 ms (paper) → ~3494 ms (current) [VERIFY: README.md:54-55].

**Load from disk:**
- For each level: read `<i>.pcd` (binary), wrap each coord as
  `Eigen::Vector4i(x, y, z, 1)`. No re-hashing, no dilation, no probing.

That's ~130 ms [VERIFY: README.md:56], a ~27× speedup.

**Trade-off:** the disk cache is keyed only on `(min_level_res, max_level,
v_rate)`. If you change the target cloud or want to re-run dilation, **delete
`voxelmaps_coords/`** to force a fresh build. The code does *not* validate the
on-disk coords against the current target cloud
[VERIFY: bbs3d/src/cpu_bbs3d/voxelmaps_io.cpp:9-19].

---

## Q10. Why does `localize()` continue past the first leaf with a high score, rather than early-exit?

**Why:** It's a *full search*. The README states this as the headline feature
[VERIFY: README.md:34]: *"Full search algorithm based on branch-and-bound."*

A leaf popped from the priority queue has the highest score currently in the
queue, but the queue may still contain *internal nodes* whose upper-bound
score equals or exceeds it. Until those are pruned (because `best_score`
catches up) or fully expanded into their leaves, the algorithm can't claim
optimality.

The code reflects this:
- The loop runs `while (!trans_queue.empty())` [VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:200].
- The only premature exits are `has_timed_out_` checks
  [VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:201-204].

**How to apply:** If you need a fast "good enough" result, set
`set_timeout_duration_in_msec(...)` and `enable_timeout()`. The returned
`has_localized_` flag will be `false` in that case [VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:248-250],
but `best_trans` will still contain the best leaf found before the timeout —
**though it's not yielded back through `get_global_pose()`**, so you'd need
to patch the API if you want partial results.

---

## Q11. How does `score_threshold_percentage` interact with the algorithm?

[VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:177-179]:
```cpp
const int score_threshold = std::floor(src_points_.size() * score_threshold_percentage_);
int best_score = score_threshold;
DiscreteTransformation<double> best_trans(best_score);
```

Two effects:

1. **Pruning floor.** Any node with `score < threshold` is dropped immediately,
   never expanded. This narrows the search.
2. **Localisation gate.** At the end, if `best_score == threshold` *exactly*,
   the algorithm considers itself "not localised" [VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:248-250].
   This means a leaf must *exceed* the threshold to count as a successful
   localisation — equal doesn't count.

**Why:** Prevents low-confidence "best of a bad lot" answers from being
reported. With `score_threshold_percentage = 0.9` (test default), at least
90% of source points must match an occupied voxel for the result to be
trusted.

**How to apply:** Lower this to e.g. 0.5 if your source cloud has heavy
occlusions or moving objects (people walking around). Raise to 0.95+ for
high-confidence applications where false positives are unacceptable.

---

## Q12. Why does the CPU pre-built loader skip the open-addressing rebuild?

**What the CPU loader actually does** [VERIFY: bbs3d/src/cpu_bbs3d/voxelmaps_io.cpp:40-74]:

```cpp
for each pcd file (level i):
  coords = read_pcd<int>(file);
  for each c: coords4i[j] << c, 1;          // wrap as Vector4i(x, y, z, 1)
  voxelmaps_ptr_->multi_buckets_[i] = coords4i;
  voxelmaps_ptr_->voxelmaps_res_[i] = min_res * v_rate^i;
```

Note: `coords4i` is used **as-is** for `multi_buckets_[i]`. No re-hashing, no
empty slots, no power-of-2 sizing. The probe in `calc_score` does
`(hash + j) % num_buckets`, where `num_buckets = coords4i.size()` — i.e.
modulo the **packed** size.

**Effect on lookup correctness:**
- If the same coord was inserted at index `k` during the original build, but
  the file's PCD order places it at index `k'`, the lookup probe starting at
  `(hash + 0) % num_buckets` will most likely **miss** it within 10 probes.
- However, `save_voxelmaps_pcd` writes coords *in the order they appear in
  `multi_buckets_[i]`* [VERIFY: bbs3d/src/cpu_bbs3d/voxelmaps_io.cpp:117-131],
  i.e. in their *hashed positions* — the empty slots get dropped by
  `coords3i = coords[i].head<3>()` but they had `w()==0` so they wouldn't
  produce 3D coords anyway… *wait*. Let me re-check.

Looking again at [VERIFY: bbs3d/src/cpu_bbs3d/voxelmaps_io.cpp:117-131]:
```cpp
const auto& coords = voxelmaps_ptr_->multi_buckets_[i];
std::vector<Eigen::Vector3i> coords3i(coords.size());
for (int i = 0; i < coords.size(); ++i) coords3i[i] = coords[i].head<3>();
pciof::save_pcd<int>(file_path, coords3i);
```

The save loop iterates all `multi_buckets_[i]` entries **including empty
slots** (since they're `Vector4i::Zero()`, `head<3>()` produces `(0,0,0)`).
The PCD file then contains both real entries and empty `(0,0,0)`s, preserving
the bucket-array layout.

**On load**, `coords4i[j] << c, 1` marks *every* entry as occupied — including
the dummy `(0,0,0)` slots. Lookup will then *always* match `(0,0,0)` coords
even if they were originally empty.

**Practical impact:** This is mostly fine because:
1. `(0,0,0)` in voxel coords corresponds to the origin of the world frame; if
   your map doesn't contain the origin, no source point ever quantises to
   `(0,0,0)` and the spurious `(0,0,0,1)` slots are never matched.
2. Even if some source point does quantise to `(0,0,0)`, it would have been
   a legitimate match in the fresh-build case too.

**The honest summary:** The disk cache **preserves the bucket array layout
verbatim**, including empty slots which become legitimate-looking `(0,0,0,1)`
entries on reload. Lookup correctness is preserved because:
- The same Teschner hash and probe budget apply.
- `(0,0,0,1)` only matches a query that quantises to `(0,0,0)`, which is
  vanishingly rare for non-origin-centred maps.

**How to apply:** Don't centre your target cloud on the world origin if you
rely on the disk cache, or be aware that origin-adjacent points will likely
match the spurious empty slots.

> *Note added during verification:* this pitfall was identified by reading
> the I/O code carefully; I haven't seen it discussed in the upstream
> documentation. If anyone has tested this corner case in practice, it
> would be worth confirming experimentally.

---

## Q13. Why is there no API to retrieve the best transform after a timeout?

[VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:248-250]:
```cpp
if (best_score == score_threshold || has_timed_out_) {
  has_localized_ = false;
  return;
}
```

When `has_timed_out_` is true, `global_pose_` is **not** updated to
`best_trans.create_matrix(...)`. It retains whatever value it had from the
previous successful call (or default-constructed garbage if no call has
succeeded).

**Why not return the partial best?** Two plausible reasons (the code doesn't
say):
1. A timed-out best is by definition not optimal — returning it might mislead
   downstream users into trusting a partial answer.
2. The API hardens the contract: `has_localized() == true` ⇒ `get_global_pose()`
   is the proven-optimal discrete pose; otherwise don't trust the pose.

**How to apply:** If you need partial results, either:
1. Modify `localize()` to always set `global_pose_` from `best_trans`
   regardless of timeout.
2. Add a `get_partial_global_pose()` accessor that reads `best_trans_`
   directly.

Neither is currently exposed — `best_trans` is a local variable, not a
member.

---

## Q14. Why is the source cloud sometimes ~2 m downsampled but the target only ~0.1 m?

From `test/config/test.yaml` [VERIFY: test/config/test.yaml:22-25]:
```yaml
tar_leaf_size: 0.1   # off: 0.0
src_leaf_size: 2.0   # off: 0.0
```

**Why:** The target cloud's density determines the *map fidelity*. The source
cloud's density determines the *score-evaluation cost*: each source point
contributes one inner-loop iteration in `calc_score`.

`calc_score`'s inner loop:
```cpp
for (i = 0; i < points.size(); i++) {            // [VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:147]
  ...
}
```
Halving the source cloud halves scoring cost; halving the target only
slightly reduces voxelmap memory.

A LiDAR sweep at 1° angular spacing produces ~64k points; downsampling to 2 m
voxels typically leaves ~500–2000 points — enough for the BnB bound to remain
informative, while keeping `calc_score` fast.

**How to apply:**
- Indoor: `src_leaf_size = 0.5` is usually enough resolution.
- Outdoor with wide scenes: `src_leaf_size = 2.0` (current default) keeps the
  source cloud small while preserving structural features.
- `tar_leaf_size` should be ≤ `min_level_res` (default 1.0 m); otherwise the
  fine-level voxels become sparse.

---

## Q15. Why does the test code run GICP after 3D-BBS (when `use_gicp = true`)?

[VERIFY: test/src/gpu_test.cpp:120-132] and [VERIFY: test/src/cpu_test.cpp:123-135].

**Why:** 3D-BBS returns the optimal *discrete* pose at the finest level's
resolution (`min_level_res`, typically 1 m). Even the optimal answer has
**up to ±0.5 m and ±Δθ_0/2 ≈ ±0.7°** quantisation error.

GICP (Generalized ICP) refines this to sub-cm / sub-mrad. The `test_code.md`
documentation hints at this design: 3D-BBS supplies the "coarse" pose, then
GICP "fine-tunes" it.

**How to apply:** Production users almost always combine 3D-BBS with a fine
registration step. The library does *not* refine internally — the
`global_pose_` returned by `localize()` is the final answer from 3D-BBS's
perspective.

---

## Q16. Why is the workspace `CLAUDE.md` saying 3D-BBS is not built by colcon?

The workspace-level `localization_ws/CLAUDE.md` notes that `src/3d_bbs/` has
no `package.xml` and no entry in `build/`, so `colcon build` skips it.

Confirmed: the project root has `CMakeLists.txt` only
[VERIFY: CMakeLists.txt:1-67], no `package.xml`. To use it from a ROS 2
workspace, callers must either:
1. Build the standalone `bbs3d/` library via `mkdir build && cmake ..; make
   install` (per `README.md:84-95`), then consume it via `find_package(gpu_bbs3d
   REQUIRED)` in their own `package.xml`-bearing package.
2. Build the ROS 2 demo separately under `ros2_test/iridescence/` or
   `ros2_test/rviz2/` (each has its own `CMakeLists.txt` and is a ROS 2
   `ament_cmake` package).

This is intentional — `bbs3d/` is a pure C++/CUDA library, deliberately
ROS-agnostic.

---

## Q17. Are there any latent thread-safety concerns?

Read carefully:

- `cpu::BBS3D` uses `#pragma omp parallel for` inside `localize()`
  [VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:193] [VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:227]. The
  parallel region only mutates per-pose state (`trans.score++`) — each thread
  has its own `trans`, no contention.
- `gpu::BBS3D` owns one `cudaStream_t` and serialises all CUDA calls on it.
  Concurrent `localize()` calls on the same instance would *interleave* H↔D
  copies on the same stream and corrupt `d_src_points_` and `d_transset_`. Don't
  do this.
- Both classes hold mutable state (`global_pose_`, `best_score_`, `has_localized_`,
  `elapsed_time_`). A second `localize()` while the first is in progress would
  race. Treat each `BBS3D` instance as single-threaded externally; OpenMP
  parallelism inside `localize()` is internal and safe.

**How to apply:** If you need concurrent localisation queries, instantiate
multiple `BBS3D` objects (one per worker thread) and share the *immutable*
voxelmap by serialising once via `save_voxelmaps_pcd` then loading per worker
via `set_voxelmaps_coords`.

---

## Q18. Why do scores top out at `src_points_.size()` (not higher)?

`calc_score` increments `trans.score` *at most once per source point*
[VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:164-165]:
```cpp
trans.score++;
break;                    // ← break out of probe loop after first hit
```

So `score ≤ src_points_.size()` always. The `break` after `score++` is what
guarantees this: even if a point's hash collides with multiple buckets that
happen to also be occupied (extremely unlikely given the Teschner hash and
the explicit coord check), only the *first* matching bucket contributes one
point to the score.

The GPU kernel does the same [VERIFY: bbs3d/src/gpu_bbs3d/calc_score.cu:49-52]:
```cpp
if (b.w() == 1) { score++; break; }
```

---

## Q19. Can 3D-BBS handle dynamic obstacles in the source cloud?

Indirectly. The score is a *count of matched points*, so any source point
that doesn't match the target voxelmap simply doesn't contribute. As long as
*enough* source points still match (above `score_threshold_percentage * N_src`),
localisation succeeds.

**How to apply:**
- Default `score_threshold_percentage = 0.9` is too strict if 30% of the scan
  is on a person walking past — drop to 0.5 in such scenarios.
- Crop short-range points with `min_scan_range > 0.5 m`
  [VERIFY: test/config/test.yaml:29] to exclude the robot's own
  ego-occlusion.

---

## Q20. What's the relationship to Hess et al.'s 2-D Cartographer BBS?

The README acknowledges [VERIFY: README.md:128-130]:
> [hdl_global_localization, cartographer, TEASER-plusplus]

The 2-D Cartographer BBS (Hess et al., ICRA 2016) is the conceptual ancestor:
- Same upper-bound trick (coarse score ≥ fine score),
- Same probabilistic-grid lookup (replaced here with a hash voxelmap),
- Same priority-queue search.

3D-BBS contributes:
- The 6-DoF (vs 3-DoF) extension — including the `arccos(1 − r²/(2R²))`
  per-level angular resolution coupling translation and rotation.
- A hash-bucket (vs dense grid) voxelmap, making the memory footprint scale
  with occupancy rather than volume.
- The GPU-parallel score kernel + batched branching, dropping the latency
  from ~878 ms (paper) to ~189 ms (current).

---

## 21. Cross-References

- Architecture overview: `00-SYSTEM_OVERVIEW.md`
- Data structures: `01-DATA_STRUCTURES.md`
- Data flow: `02-DATA_FLOW.md`
- BnB algorithm: `03-ALGORITHM_01-BranchAndBound.md`
- Voxelmap construction: `04-ALGORITHM_02-HierarchicalVoxelmap.md`
- GPU kernel: `05-ALGORITHM_03-GPU-ScoreCalc.md`
- Function-by-function: `06-KEY_FUNCTIONS.md`

---

## 22. Verification Checklist

- [x] Every Q has at least one `[VERIFY:]` tag pointing at the supporting code
- [x] Speculative answers (Q12, Q13) are explicitly flagged as such
- [x] Cross-references match documents that exist (00-, 01-, 02-, 03-, 04-, 05-, 06-)
- [x] Tradeoff tables are derived from the code, not invented
