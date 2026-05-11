# 3D-BBS — Key Function Analysis

> Line-by-line walk-through of the most important functions. Each section is
> self-contained — you should be able to cross-check every claim against the
> referenced source line without leaving this document.

Functions analysed below:

1. `cpu::BBS3D::localize()` — main BnB loop (CPU)
2. `gpu::BBS3D::localize()` — main BnB loop (GPU)
3. `cpu::BBS3D::calc_score()` — single-pose score evaluation (CPU)
4. `gpu::calc_scores_kernel(...)` — single-pose score evaluation (GPU)
5. `cpu::BBS3D::calc_angular_info()` — per-level Δθ
6. `cpu::BBS3D::create_init_transset()` — initial transset
7. `cpu::VoxelMaps::create_voxelmaps()` — pyramid construction
8. `cpu::VoxelMaps::create_hash_buckets()` — open-addressing builder
9. `DiscreteTransformation::branch()` — child expansion
10. `DiscreteTransformation::create_matrix()` — discrete → continuous SE(3)

---

## 1. `cpu::BBS3D::localize()`

File: `bbs3d/src/cpu_bbs3d/bbs3d.cpp`, lines **170–258**
[VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:170-258]

### 1.1 Signature

```cpp
void BBS3D::localize();
```
No parameters; results read through `get_global_pose() / get_best_score() / has_localized() / has_timed_out() / get_elapsed_time()`.

### 1.2 Annotated body

| Line | Code | Effect |
|---:|---|---|
| 172–174 | start clock, compute `time_limit`, reset `has_timed_out_` | Start of frame |
| 176–179 | `best_score = floor(N_src · pct)`; `best_trans(best_score)` | Initial best is the threshold "dummy" — any leaf must beat it to be considered localised |
| 182–185 | `max_bucket_scan_count`, `max_level`; populate `ang_info_vec` | Cache scalars |
| 186 | `init_transset = create_init_transset(ang_info_vec[max_level])` | Enumerate top-level translation × rotation grid |
| 189–192 | Read `top_buckets`, `init_trans_res`, `rpy_res`, `min_rpy` for top level | Hoisted outside the parallel for |
| 193–196 | `#pragma omp parallel for num_threads(num_threads_)` over `init_transset`, calling `calc_score(...)` | First parallel scoring pass — fills `.score` for every top-level pose |
| 198 | `std::priority_queue<DT<double>> trans_queue(init_transset.begin(), init_transset.end())` | Range constructor; `O(N)` heapify |
| 200 | `while (!trans_queue.empty())` | BnB loop |
| 201–204 | If `use_timeout_` and clock exceeded `time_limit`, set `has_timed_out_=true; break` | Optional bail |
| 206–207 | `auto trans = trans_queue.top(); trans_queue.pop();` | Pop the highest-score node |
| 210–212 | `if (trans.score < best_score) continue;` | Prune (Invariant 2; see `03-ALGORITHM_01-BranchAndBound.md` §5) |
| 214–216 | If `is_leaf()`: `best_trans = trans; best_score = trans.score;` | Record best leaf |
| 218–222 | Else: read child-level `num_division/rpy_res/min_rpy` from `ang_info_vec`, get `children = trans.branch(child_level, v_rate_, num_division)` | Branch the node |
| 224–225 | Fetch `buckets` and `trans_res` for child level | Hoisted refs |
| 227–230 | `#pragma omp parallel for` over `children`, calling `calc_score(...)` | Parallel scoring of children |
| 232–239 | For each child: if `score < best_score` skip, else `push` | Push survivors back to heap |
| 244–245 | After loop ends: compute `elapsed_time_` in ms | Wall time bookkeeping |
| 248–251 | If `best_score == score_threshold` or `has_timed_out_`: `has_localized_ = false; return;` | Fail path |
| 253–254 | `global_pose_ = best_trans.create_matrix(min_res, ang_info_vec[0].rpy_res, ang_info_vec[0].min_rpy)` | Build 4×4 transform at level 0 |
| 255–257 | `best_score_ = best_score; has_timed_out_ = false; has_localized_ = true;` | Success bookkeeping |

### 1.3 Subtle points

- **Reset of `has_timed_out_` (line 174):** A leftover `true` from a previous
  call would prevent subsequent successful runs from being marked
  `has_localized_=true`. The reset is mandatory if the same `BBS3D` instance
  is used across frames.
- **`score_threshold` semantics:** When `score_threshold_percentage_ == 0`,
  `score_threshold == 0`, and `best_score` starts at 0. The `==` check at line
  248 therefore catches the "no leaf was ever popped" case as a failure.
- **Reference lifetime in the parallel for:** `rpy_res` and `min_rpy` are
  bound to specific entries of `ang_info_vec` *before* the parallel region;
  since `ang_info_vec` isn't mutated inside the loop, all threads see the same
  values without further locking.

---

## 2. `gpu::BBS3D::localize()`

File: `bbs3d/src/gpu_bbs3d/bbs3d.cu`, lines **150–240**
[VERIFY: bbs3d/src/gpu_bbs3d/bbs3d.cu:150-240]

### 2.1 Notable structural differences from the CPU version

1. **No `num_threads_`** — parallelism is the kernel's grid.
2. **No per-iter `#pragma omp parallel for`** — instead a host-side
   `branch_stock` defers branching, then a single batched kernel via
   `calc_scores(...)`.
3. **`ang_info_vec` is uploaded once per frame** via
   `cudaMemcpyAsync(d_ang_info_vec, ...)` on lines 167–172.
4. **Initial transset is scored in a single batched call**: `init_transset_output = calc_scores(init_transset, d_ang_info_vec);` (line 177).
5. **End-of-loop flush**: when `trans_queue.empty() && !branch_stock.empty()`,
  a final `calc_scores(branch_stock, ...)` is executed (lines 193–200). Without
  this, deferred children would never be scored.

### 2.2 Annotated body

| Line | Code | Effect |
|---:|---|---|
| 152–154 | start clock, compute `time_limit`, reset `has_timed_out_` | Start of frame |
| 157–160 | `best_score = floor(N_src · pct)`; `best_trans(best_score)` | Same as CPU |
| 163–165 | Compute `ang_info_vec`; convert to `thrust::device_vector` | Host build of the per-level Δθ |
| 167–172 | `cudaMemcpyAsync(d_ang_info_vec, ang_info_vec, ...)` on `stream` | One H→D copy |
| 174 | `init_transset = create_init_transset(...)` | Enumerate top-level grid |
| 177 | `init_transset_output = calc_scores(init_transset, d_ang_info_vec)` | First batched score pass |
| 179 | `priority_queue trans_queue(init_transset_output.begin(), init_transset_output.end())` | Seed heap |
| 181–182 | `branch_stock.reserve(branch_copy_size_)` | Defer-branch buffer (default 10000) |
| 183 | `while (!trans_queue.empty())` | BnB loop |
| 184–187 | Timeout check | Same as CPU |
| 189–190 | Pop top | Same as CPU |
| 193–200 | If queue would empty out and stock has unscored items: flush stock, push survivors | **End-of-loop flush** (GPU-only) |
| 202–205 | Prune | Same as CPU |
| 207–209 | Leaf: update best | Same as CPU |
| 210–213 | Non-leaf: `trans.branch(branch_stock, child_level, v_rate, num_division)` — appends to stock | Deferred branching |
| 215–222 | If `branch_stock.size() >= branch_copy_size_`: flush, push survivors, clear stock | **Size-triggered flush** (GPU-only) |
| 226–227 | Compute `elapsed_time_` | Wall time |
| 230–233 | Failure path | Same as CPU |
| 235–239 | Success path: `global_pose_ = best_trans.create_matrix(min_res, ang_info_vec[0]...)` | Build 4×4 transform |

---

## 3. `cpu::BBS3D::calc_score(...)`

File: `bbs3d/src/cpu_bbs3d/bbs3d.cpp`, lines **134–168**
[VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:134-168]

### 3.1 Signature

```cpp
void BBS3D::calc_score(
    DiscreteTransformation<double>& trans,
    const double trans_res,
    const Eigen::Vector3d& rpy_res,
    const Eigen::Vector3d& min_rpy,
    const std::vector<Eigen::Vector4i>& buckets,
    const int max_bucket_scan_count,
    const std::vector<Eigen::Vector3d>& points);
```
[VERIFY: bbs3d/include/cpu_bbs3d/bbs3d.hpp:97-104]

### 3.2 Annotated body

| Line | Code | Effect |
|---:|---|---|
| 142–143 | `num_buckets = buckets.size(); inv_res = 1.0/trans_res` | Hoisted |
| 144–145 | `transform = trans.create_matrix(trans_res, rpy_res, min_rpy)` | Build 4×4 transform once per pose |
| 147 | `for (i = 0; i < points.size(); i++)` | Loop over source points |
| 148 | `transed_point = transform * points[i]` | Eigen does T·Rz·Ry·Rx by composition order |
| 149 | `coord = (transed_point * inv_res).floor().cast<int>()` | Quantise |
| 150 | `hash = (coord[0]*73856093) ^ (coord[1]*19349669) ^ (coord[2]*83492791)` | Teschner |
| 152 | `for (j = 0; j < max_bucket_scan_count; j++)` | Probe loop (default budget 10) |
| 153 | `bucket_index = (hash + j) % num_buckets` | Linear probing |
| 154 | `bucket = buckets[bucket_index]` | Read |
| 156–158 | `if (bucket.w() == 0) break;` | **Empty-slot fast exit (CPU-only optimisation)** |
| 160–162 | `if (bucket.xyz() != coord) continue;` | Hash collision; probe next |
| 164–165 | `trans.score++; break;` | Hit |

### 3.3 Notes

- `trans` is passed by reference and mutated in place (`trans.score++`).
- The lookup probe budget is *exactly* `max_bucket_scan_count` — no extension
  beyond that even if the bucket array is half-full of collisions.
- Note that `bucket.w() == 0` is consulted *before* the coord check on each
  probe. This means a `w==0` slot mid-probe terminates the probe even if the
  next-probe slot might match — which is correct because the *builder*
  guarantees no key skips a `w==0` slot (it always inserts into the first
  empty slot it finds, line 82 of `voxelmaps.cpp`).

---

## 4. `gpu::calc_scores_kernel(...)`

Already analysed in detail in `05-ALGORITHM_03-GPU-ScoreCalc.md` §1.
Source: [VERIFY: bbs3d/src/gpu_bbs3d/calc_score.cu:7-56].

Key facts not repeated:

- Per-thread output is `trans.score` (set on line 55).
- Both `coord` and `w` are read on every probe — no early break on empty slot.
- Rotation composition matches `DiscreteTransformation::create_matrix()`.

---

## 5. `cpu::BBS3D::calc_angular_info(...)`

File: `bbs3d/src/cpu_bbs3d/bbs3d.cpp`, lines **63–107**
[VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:63-107]

### 5.1 Outline

```cpp
void calc_angular_info(std::vector<AngularInfo>& ang_info_vec) {
  max_norm = max ||p|| over src_points;            // line 64-70
  for (i = max_level; i >= 0; --i) {               // top-down
      cosine  = 1 - res[i]² / (2·max_norm²);
      ori_res = floor(arccos(max(cosine,-1)) * 1e4) / 1e4;     // line 73-76
      rpy_res_temp.{x,y,z} = (ori_res ≤ |max-min|) ? ori_res : 0;
      max_rpy_piece = (i == max_level) ? max_rpy - min_rpy
                                       : parent.rpy_res or full range;
      num_division.{x,y,z} = ceil(max_rpy_piece / rpy_res_temp);
      ang_info[i].rpy_res = num_division == 1 ? 0 : max_rpy_piece / num_division;
      ang_info[i].min_rpy = ang_info[i].rpy_res == 0 ? 0 : min_rpy;
  }
}
```

### 5.2 Annotated body

| Line | Code | Effect |
|---:|---|---|
| 64–70 | First pass over `src_points_` to find `max_norm` | Worst-case rotation arm |
| 73 | `cosine = 1 - (res[i]²/max_norm²)·0.5` | Cosine of natural Δθ |
| 74 | `ori_res = arccos(max(cosine,-1))` | Bounded to avoid NaN |
| 75 | `ori_res = floor(ori_res * 10000) / 10000` | Truncate to 4 decimal places (deterministic across runs) |
| 78–80 | Per-axis: keep `ori_res` if it fits in the user's range, else 0 | "Axis disabled" flag |
| 83–89 | Compute `max_rpy_piece` based on parent's `rpy_res` (or full range if top level) | Range a single child must cover |
| 92–95 | `num_division = ceil(max_rpy_piece / rpy_res_temp)` (or 1 if 0) | Number of bins per axis |
| 99–101 | `ang_info[i].rpy_res = max_rpy_piece / num_division` (or 0 if num_division == 1) | Actual per-bin width — note this is *not* exactly `ori_res`! |
| 103–105 | `ang_info[i].min_rpy = (rpy_res == 0) ? 0 : min_rpy` | Carry the user's lower bound only when this axis is searched |

### 5.3 Subtleties

- **Why `floor(x · 1e4) / 1e4`?** To get a stable, deterministic angular
  resolution that doesn't drift across runs because of FP rounding (e.g.
  `arccos` slightly different on different math libs). Truncating to 4
  decimals discards sub-mrad jitter.
- **`ang_info[i].rpy_res = max_rpy_piece / num_division`** (not `ori_res`) — the
  bin width is *slightly smaller than* `ori_res` because `num_division =
  ceil(...)` overshoots a bit. This guarantees the angular grid fits within
  the user's range without truncation, at the cost of slightly tighter
  rotations.
- **GPU mirror** at [VERIFY: bbs3d/src/gpu_bbs3d/bbs3d.cu:79-123] uses
  `float`s, `std::max` replaced by `max` (CUDA-friendly), but otherwise the
  same algorithm.

---

## 6. `cpu::BBS3D::create_init_transset(...)`

File: `bbs3d/src/cpu_bbs3d/bbs3d.cpp`, lines **109–132**
[VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:109-132]

```cpp
std::vector<DT<double>> create_init_transset(const AngularInfo& init_ang_info) {
  init_transset_size = (tx_max - tx_min + 1) * (ty_max - ty_min + 1)
                     * (tz_max - tz_min + 1)
                     * init_ang_info.num_division.x()
                     * init_ang_info.num_division.y()
                     * init_ang_info.num_division.z();
  transset.reserve(init_transset_size);
  for (tx ∈ [tx_min, tx_max]) for (ty ...) for (tz ...)
    for (roll < num_division.x()) for (pitch < num_division.y()) for (yaw < num_division.z())
      transset.emplace_back(DT<double>(0, max_level, tx, ty, tz, roll, pitch, yaw));
  return transset;
}
```

### 6.1 Why include both endpoints (`<= tx_max`)?

`init_tx_range_` was constructed as
`pair(floor(min_xyz.x / top_res), ceil(max_xyz.x / top_res))`
[VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:58], so the *integer* endpoints lie
inclusive of the bounding box. Looping `tx <= tx_max` covers the entire
box including its `max` corner.

### 6.2 Why start every score at 0?

Children's scores are written by `calc_score` before they are pushed into the
queue, so the initial value is overwritten. But this initial 0 still appears
briefly in the heap-pre-population step: every initial transset entry passes
through `calc_score` (line 195) *before* the `priority_queue` is constructed
(line 198), so the heap sees real scores.

---

## 7. `cpu::VoxelMaps::create_voxelmaps(...)`

File: `bbs3d/src/cpu_bbs3d/voxelmaps.cpp`, lines **9–50**
[VERIFY: bbs3d/src/cpu_bbs3d/voxelmaps.cpp:9-50]

Already analysed in detail in `04-ALGORITHM_02-HierarchicalVoxelmap.md` §2.

Single-pass per-level loop:

```cpp
for (i = 0; i < max_level_ + 1; i++) {
  unordered_voxelmap.clear();
  coords = points |> floor(/ resolution);           // line 17-20
  for c in coords:                                  // line 23-42
    if c not in map:        map[c] = 1;
    elif map[c] == -1:      map[c] = 1;
    elif map[c] ==  1:      continue;
    for n in 7_neighbours(c):
      if n not in map:      map[n] = -1;
  multi_buckets_[i]  = create_hash_buckets(map);   // line 44
  voxelmaps_res_[i]  = resolution;
  resolution        *= v_rate;                      // line 48
}
```

### 7.1 Memory churn note

`unordered_voxelmap` is freshly allocated per iteration. For a 6-level pyramid
on a 100k-point cloud, that's up to 6 large `unordered_map` allocations.
A minor optimisation would `clear()` and reuse, but the current code is fine
because construction happens once per session.

---

## 8. `cpu::VoxelMaps::create_hash_buckets(...)`

File: `bbs3d/src/cpu_bbs3d/voxelmaps.cpp`, lines **65–97**
[VERIFY: bbs3d/src/cpu_bbs3d/voxelmaps.cpp:65-97]

Already analysed in detail in `04-ALGORITHM_02-HierarchicalVoxelmap.md` §4.
Key facts:

- Outer loop: `num_buckets = |m| · 2^k`, `k = 0, 1, ..., 4` (until `16·|m|`).
- Inner loop: per-key linear probe up to `max_bucket_scan_count` (default 10).
- Exit condition: `success_rate > 0.999` or upper bound reached.
- Output buckets use `Eigen::Vector4i = (x, y, z, 1)` for occupied, `(0,0,0,0)`
  for empty.

---

## 9. `DiscreteTransformation::branch(...)`

File: `bbs3d/include/discrete_transformation/discrete_transformation.hpp`
[VERIFY: bbs3d/include/discrete_transformation/discrete_transformation.hpp:41-90]

Two overloads. Both emit `v_rate^3 · num_division.x · num_division.y · num_division.z`
children:

```cpp
for (i = 0; i < v_rate; ++i)
 for (j = 0; j < v_rate; ++j)
  for (k = 0; k < v_rate; ++k)
   for (l = 0; l < num_division.x(); ++l)
    for (m = 0; m < num_division.y(); ++m)
     for (n = 0; n < num_division.z(); ++n)
       out.emplace_back(DT(
         0,
         child_level,
         x*v_rate + i, y*v_rate + j, z*v_rate + k,
         roll*num_division.x() + l, pitch*num_division.y() + m, yaw*num_division.z() + n));
```

### 9.1 Translation index transformation

`x_child = x_parent · v_rate + i`, `i ∈ [0, v_rate)`. This is the standard
**octree** decomposition: a parent index `x` at resolution `r_parent` covers
exactly `[x, x+1) · r_parent`, which equals `[x · v_rate, (x+1) · v_rate) · r_child`.
Sub-dividing into `v_rate` intervals along each axis yields the `v_rate^3`
children.

### 9.2 Angle index transformation

`roll_child = roll_parent · num_division.x() + l`, `l ∈ [0, num_division.x())`.
This is *not* an octree split — the rotation grid is recomputed per level by
`calc_angular_info` and may have a different `num_division` per level. The
multiplication `roll_parent · num_division_child` means the parent's
*coarse-resolution* index maps to the `num_division.x()` consecutive child
indices that fall within its angular range.

---

## 10. `DiscreteTransformation::create_matrix(...)`

File: `bbs3d/include/discrete_transformation/discrete_transformation.hpp`, lines **33–39**
[VERIFY: bbs3d/include/discrete_transformation/discrete_transformation.hpp:33-39]

```cpp
Matrix4<T> create_matrix(T trans_res,
                         const Vector3<T>& rpy_res,
                         const Vector3<T>& min_rpy) {
  Eigen::Translation<T,3> translation(x*trans_res, y*trans_res, z*trans_res);
  Eigen::AngleAxis<T> rollAngle (roll *rpy_res.x() + min_rpy.x(), Vector3<T>::UnitX());
  Eigen::AngleAxis<T> pitchAngle(pitch*rpy_res.y() + min_rpy.y(), Vector3<T>::UnitY());
  Eigen::AngleAxis<T> yawAngle  (yaw  *rpy_res.z() + min_rpy.z(), Vector3<T>::UnitZ());
  return (translation * yawAngle * pitchAngle * rollAngle).matrix();
}
```

### 10.1 Convention

The resulting matrix `M = T · Rz · Ry · Rx`, applied as `p_world = M · p_body`.
In Eigen's `Affine` math, multiplying transforms in order
`A * B` produces `A(B(x))`, so this is intrinsic-XYZ Euler (rotate
body-frame X first, then Y, then Z), with translation applied last in world frame.

The GPU kernel constructs the same rotation manually
[VERIFY: bbs3d/src/gpu_bbs3d/calc_score.cu:30-34]:
```cpp
R = Rz · Ry · Rx;
p_world = R · p + t;
```
which is equivalent: `p_world = R · p + t = T · R · p`.

### 10.2 What `rpy_res = 0` means in this formula

When `rpy_res.{x,y,z}() == 0` for a given axis, both `roll * rpy_res.x()` and
`min_rpy.x()` are 0 (the latter is forced to 0 by `calc_angular_info` line 103).
So `rollAngle = AngleAxis<T>(0, X)` = identity rotation. This is how
"unsearched" axes are handled — they contribute no rotation at any level.

---

## 11. Function Call Graph

```
caller (test/ros2)
    |
    +--> BBS3D::set_tar_points
    |       \--> VoxelMaps::create_voxelmaps
    |               +--> create_neighbor_coords
    |               \--> create_hash_buckets
    |
    +--> BBS3D::set_src_points (host + (GPU) cudaMemcpyAsync)
    |
    +--> BBS3D::set_trans_search_range
    |
    +--> BBS3D::localize
            +--> calc_angular_info
            +--> create_init_transset
            +--> calc_score  (CPU, OpenMP parallel_for)
            |   OR
            |   calc_scores --> calc_scores_kernel  (GPU)
            +--> [BnB loop]
            |     +--> DiscreteTransformation::branch
            |     +--> calc_score / calc_scores
            \--> best_trans.create_matrix
```

---

## 12. Verification Checklist

- [x] Every function in this document has its lines and file references verified
- [x] Both CPU and GPU variants of `localize` and the scoring functions covered
- [x] Subtle behaviours (deterministic rounding, axis-disabled handling) called out
- [x] Call graph matches what is actually invoked in the tested entry points
  (test/src/cpu_test.cpp, test/src/gpu_test.cpp, ros2_test/...)
