# 3D-BBS — Data Structures (数据结构详解)

> Every struct/class field below was read from source. Speculation is marked
> *(inferred from usage)* and otherwise the field carries a `[VERIFY:]` tag to its
> declaration. Companion: `00-SYSTEM_OVERVIEW.md`.

---

## 1. `DiscreteTransformation<T>` — BnB Search-Tree Node

> Header-only template, shared verbatim by CPU (`T=double`) and GPU (`T=float`).
> Declared at [VERIFY: bbs3d/include/discrete_transformation/discrete_transformation.hpp:11].

```cpp
template <typename T>
class DiscreteTransformation {
public:
  int score;   // # of source points that hit an occupied bucket
  int level;   // current pyramid level; 0 == finest (leaf)
  int x, y, z; // integer voxel index in translation grid at this level
  int roll, pitch, yaw; // integer index in angular grid at this level
};
```

[VERIFY: bbs3d/include/discrete_transformation/discrete_transformation.hpp:92-100]

### 1.1 Field semantics

| Field | Type | Meaning |
|---|---|---|
| `score` | `int` | Number of source points whose transformed coord hits an occupied voxel at `level`. Larger = better. [VERIFY: bbs3d/include/discrete_transformation/discrete_transformation.hpp:93] |
| `level` | `int` | Pyramid level: `max_level` (coarsest) → `0` (finest). `is_leaf()` returns `level == 0` [VERIFY: bbs3d/include/discrete_transformation/discrete_transformation.hpp:31]. |
| `x, y, z` | `int` | Quantised translation. Continuous translation is `(x,y,z) * trans_res[level]` [VERIFY: bbs3d/include/discrete_transformation/discrete_transformation.hpp:34]. |
| `roll, pitch, yaw` | `int` | Quantised orientation. Continuous angle is `idx * rpy_res[level] + min_rpy[level]` [VERIFY: bbs3d/include/discrete_transformation/discrete_transformation.hpp:35-37]. |

### 1.2 Ordering — used by `std::priority_queue`

```cpp
bool operator<(const DiscreteTransformation& rhs) const {
  return score < rhs.score;   // [VERIFY: bbs3d/include/discrete_transformation/discrete_transformation.hpp:29]
}
```

`std::priority_queue` is a max-heap by default, so the node with the **highest
score is popped first** — the upper-bound property required for branch-and-bound.

### 1.3 `create_matrix(...)` — discrete → continuous SE(3)

```cpp
Matrix4<T> create_matrix(T trans_res, const Vector3<T>& rpy_res, const Vector3<T>& min_rpy) {
  Eigen::Translation<T,3> translation(x*trans_res, y*trans_res, z*trans_res);
  Eigen::AngleAxis<T>     rollAngle (roll *rpy_res.x() + min_rpy.x(), Vector3<T>::UnitX());
  Eigen::AngleAxis<T>     pitchAngle(pitch*rpy_res.y() + min_rpy.y(), Vector3<T>::UnitY());
  Eigen::AngleAxis<T>     yawAngle  (yaw  *rpy_res.z() + min_rpy.z(), Vector3<T>::UnitZ());
  return (translation * yawAngle * pitchAngle * rollAngle).matrix();
}
```
[VERIFY: bbs3d/include/discrete_transformation/discrete_transformation.hpp:33-39]

Composition order is **T · Rz · Ry · Rx**, applied as `p_world = T · Rz · Ry · Rx · p_body`.
This matches the rotation order used inside the GPU kernel
[VERIFY: bbs3d/src/gpu_bbs3d/calc_score.cu:31-33] (`Rz * Ry * Rx`).

### 1.4 `branch(...)` — 6-DoF child expansion

Two overloads exist:

- In-place: `void branch(std::vector<DT>& out, child_level, v_rate, num_division)`
  [VERIFY: bbs3d/include/discrete_transformation/discrete_transformation.hpp:41-63]
- Return-by-value: `std::vector<DT> branch(child_level, v_rate, num_division)`
  [VERIFY: bbs3d/include/discrete_transformation/discrete_transformation.hpp:65-90]

Both emit a 6-nested loop producing
`v_rate^3 × num_division.x * num_division.y * num_division.z` children:

```cpp
b.emplace_back(DiscreteTransformation(
    0,                              // child score starts at 0
    child_level,                    // = parent.level - 1
    x * v_rate + i,                 // i in [0, v_rate)
    y * v_rate + j,
    z * v_rate + k,
    roll  * num_division.x() + l,   // l in [0, num_division.x())
    pitch * num_division.y() + m,
    yaw   * num_division.z() + n));
```
[VERIFY: bbs3d/include/discrete_transformation/discrete_transformation.hpp:48-57]

Translation index multiplied by `v_rate` (= 2.0 default → integer `2`); angle
indices multiplied by their *level-specific* division factor.

**Why two scaling factors?** Translation resolution shrinks by a fixed
`v_rate` each step, so each parent voxel splits into `v_rate^3` children.
Angular resolution is recomputed every level from `calc_angular_info`
([VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:63-107]) because it depends on the
current translation resolution and the source-cloud radius (see
`03-ALGORITHM_01-BranchAndBound.md` §3).

---

## 2. `cpu::AngularInfo` / `gpu::AngularInfo`

Per-level angular discretisation cache. Computed once at the top of every
`localize()` call.

```cpp
// CPU [VERIFY: bbs3d/include/cpu_bbs3d/bbs3d.hpp:15-19]
struct AngularInfo {
  Eigen::Vector3i num_division; // # angle bins per axis at this level
  Eigen::Vector3d rpy_res;      // angular step [rad] per axis at this level
  Eigen::Vector3d min_rpy;      // lower bound [rad] per axis at this level
};
```

```cpp
// GPU [VERIFY: bbs3d/include/gpu_bbs3d/bbs3d.cuh:19-23]
struct AngularInfo {
  Eigen::Vector3i num_division;
  Eigen::Vector3f rpy_res;
  Eigen::Vector3f min_rpy;
};
```

A `std::vector<AngularInfo>` of size `max_level+1` is filled by
`calc_angular_info`. Notable invariants enforced there:

- If the per-axis range is smaller than the natural resolution, `rpy_res` for
  that axis is set to `0` and `num_division = 1` (axis is *not* subdivided)
  [VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:78-95].
- When `rpy_res == 0`, `min_rpy` for that axis is forced to `0` to avoid an
  uninitialised offset in `create_matrix` [VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:103-105].

---

## 3. `cpu::VoxelMaps` / `gpu::VoxelMaps`

### 3.1 Hash bucket = `Eigen::Vector4i`

A bucket is **one int4** storing the integer voxel coord (x,y,z) plus a 1-bit
occupied flag in the 4th lane:

```cpp
Eigen::Vector4i coord;
coord << voxel.first.x(), voxel.first.y(), voxel.first.z(), 1;
```
[VERIFY: bbs3d/src/cpu_bbs3d/voxelmaps.cpp:74-75] (identical pattern at [VERIFY: bbs3d/src/gpu_bbs3d/voxelmaps.cu:76-77]).

| Lane | Meaning | Notes |
|---|---|---|
| `[0]` (`x()`) | Voxel index along X | Floor-quantised at `resolution` |
| `[1]` (`y()`) | Voxel index along Y |  |
| `[2]` (`z()`) | Voxel index along Z |  |
| `[3]` (`w()`) | Occupancy flag (0 empty / 1 set) | Used as empty-slot probe sentinel |

Used in three places:

- Construction: empty buckets are zero-initialised, then probed via open addressing [VERIFY: bbs3d/src/cpu_bbs3d/voxelmaps.cpp:70-90].
- CPU lookup: an empty bucket terminates the probe early [VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:156-158].
- GPU lookup: the (0,0,0,0) bucket cannot match any non-origin coord, so the kernel just `continue`s [VERIFY: bbs3d/src/gpu_bbs3d/calc_score.cu:45-47].

### 3.2 `cpu::VoxelMaps`

```cpp
class VoxelMaps {
public:
  using UnorderedVoxelMap = std::unordered_map<Eigen::Vector3i, int, VectorHash, VctorEqual>;
  using Buckets           = std::vector<Eigen::Vector4i>;

  std::vector<Buckets> multi_buckets_;   // per-level hash buckets
  std::vector<double>  voxelmaps_res_;   // per-level resolution [m]

private:
  double min_level_res_;       // default 1.0  [m]
  int    max_level_;           // default 6
  int    max_bucket_scan_count_; // default 10 (linear-probe budget)
};
```
[VERIFY: bbs3d/include/cpu_bbs3d/voxelmaps.hpp:28-58]

Default values are set in the constructor:

```cpp
VoxelMaps::VoxelMaps() : min_level_res_(1.0), max_level_(6), max_bucket_scan_count_(10) {}
```
[VERIFY: bbs3d/src/cpu_bbs3d/voxelmaps.cpp:5]

### 3.3 `gpu::VoxelMaps`

Adds a device-side mirror plus per-level descriptor `VoxelMapInfo`:

```cpp
struct VoxelMapInfo {
  int   num_buckets;
  int   max_bucket_scan_count;
  float res;
  float inv_res;
};
```
[VERIFY: bbs3d/include/gpu_bbs3d/voxelmaps.cuh:12-17]

```cpp
class VoxelMaps {
public:
  std::vector<DeviceBuckets> d_multi_buckets_;            // one device vector per level
  thrust::device_vector<Eigen::Vector4i*> d_multi_buckets_ptrs_; // raw pointers, indexed by level
  std::vector<VoxelMapInfo>  voxelmaps_info_;             // host copy of per-level info
  thrust::device_vector<VoxelMapInfo> d_voxelmaps_info_;  // device copy
};
```
[VERIFY: bbs3d/include/gpu_bbs3d/voxelmaps.cuh:63-69]

The `d_multi_buckets_ptrs_` indirection lets the kernel index into the right
bucket array using only `trans.level` and a single global-memory read
[VERIFY: bbs3d/src/gpu_bbs3d/calc_score.cu:23].

### 3.4 The `UnorderedVoxelMap` int-flag legend (construction-only)

During build, `UnorderedVoxelMap` (a `std::unordered_map<Vector3i,int,...>`)
stores a **tri-state** flag for each candidate voxel:

| Flag | Meaning |
|---|---|
| `1`  | Voxel is occupied by ≥ 1 point |
| `-1` | Voxel is a neighbour of an occupied voxel (will be promoted to `1` if a point ever lands there) |
| (absent) | Not in the map |

Coverage logic [VERIFY: bbs3d/src/cpu_bbs3d/voxelmaps.cpp:23-42]:

```cpp
if (unordered_voxelmap.count(coord) == 0) unordered_voxelmap[coord] = 1;
else if (unordered_voxelmap[coord] == -1) unordered_voxelmap[coord] = 1;
else if (unordered_voxelmap[coord] ==  1) continue;
// emit 7 neighbour coords with flag -1 if not present
```

Once construction is done **every entry** (flag `1` or `-1`) is hashed into the
output buckets — the `-1` neighbours are what makes 3D-BBS robust to small
discretisation errors near voxel boundaries (see `04-ALGORITHM_02-HierarchicalVoxelmap.md` §3).

### 3.5 `create_neighbor_coords(...)` — exactly 7 neighbours

```cpp
neighbor_coord.emplace_back(vec.x()-1, vec.y()-1, vec.z()  );
neighbor_coord.emplace_back(vec.x()-1, vec.y(),   vec.z()  );
neighbor_coord.emplace_back(vec.x(),   vec.y()-1, vec.z()  );
neighbor_coord.emplace_back(vec.x()-1, vec.y()-1, vec.z()-1);
neighbor_coord.emplace_back(vec.x()-1, vec.y(),   vec.z()-1);
neighbor_coord.emplace_back(vec.x(),   vec.y()-1, vec.z()-1);
neighbor_coord.emplace_back(vec.x(),   vec.y(),   vec.z()-1);
```
[VERIFY: bbs3d/src/cpu_bbs3d/voxelmaps.cpp:52-63] and identical pattern at [VERIFY: bbs3d/src/gpu_bbs3d/voxelmaps.cu:54-65].

These are exactly the **"-1" corners of the surrounding 2×2×2 block** —
covering the 7 voxels reached by simultaneously *decreasing* zero or more of
x/y/z (the 8th, the voxel itself, is already in the map). This gives a one-sided
dilation that is matched on the lookup side by `(point/res).floor()`
quantisation (`floor` makes a point land at the negative-most cell of the
neighbourhood; the dilation extends the map toward the negative direction).

---

## 4. `cpu::BBS3D` — Class Layout

```cpp
class BBS3D {
private:
  Eigen::Matrix4d global_pose_;
  bool   has_timed_out_, has_localized_;
  double elapsed_time_;

  std::vector<Eigen::Vector3d> src_points_;
  std::unique_ptr<VoxelMaps>   voxelmaps_ptr_;
  std::string                  voxelmaps_folder_name_;   // "voxelmaps_coords"

  double v_rate_, inv_v_rate_;   // voxel expansion rate (default 2.0)
  int    num_threads_;           // OpenMP threads (default 4)

  int    best_score_;
  double score_threshold_percentage_;   // best_score must exceed this fraction
  bool   use_timeout_;
  std::chrono::milliseconds timeout_duration_;   // default 10 s

  Eigen::Vector3d min_xyz_, max_xyz_;            // translation search box [m]
  Eigen::Vector3d min_rpy_, max_rpy_;            // angle search box [rad]
  std::pair<int,int> init_tx_range_, init_ty_range_, init_tz_range_;
};
```
[VERIFY: bbs3d/include/cpu_bbs3d/bbs3d.hpp:107-130]

Default values (from constructor) [VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:5-17]:

| Member | Default |
|---|---|
| `v_rate_` | `2.0` |
| `num_threads_` | `4` |
| `score_threshold_percentage_` | `0.0` |
| `use_timeout_` | `false` |
| `timeout_duration_` | `10000` ms |
| `voxelmaps_folder_name_` | `"voxelmaps_coords"` |
| `min_rpy_` | `(-0.02, -0.02, 0.0)` |
| `max_rpy_` | `(0.02, 0.02, 2π)` |

The default RPY box (±0.02 rad ≈ ±1.15° on roll/pitch, 0..2π on yaw) is the
"IMU-aligned scan" assumption.

---

## 5. `gpu::BBS3D` — Class Layout

```cpp
class BBS3D {
private:
  Eigen::Matrix4f global_pose_;
  bool has_timed_out_, has_localized_;
  double elapsed_time_;

  cudaStream_t stream;                         // dedicated stream for all H↔D copies
  std::vector<Eigen::Vector3f>             src_points_;
  thrust::device_vector<Eigen::Vector3f>   d_src_points_;
  std::unique_ptr<VoxelMaps> voxelmaps_ptr_;
  std::string voxelmaps_folder_name_;

  float v_rate_, inv_v_rate_;                  // default 2.0

  int   branch_copy_size_, best_score_;        // branch_copy_size_ default 10000
  double score_threshold_percentage_;
  bool  use_timeout_;
  std::chrono::milliseconds timeout_duration_;
  Eigen::Vector3f min_xyz_, max_xyz_, min_rpy_, max_rpy_;
  std::pair<int,int> init_tx_range_, init_ty_range_, init_tz_range_;
};
```
[VERIFY: bbs3d/include/gpu_bbs3d/bbs3d.cuh:101-128]

Defaults [VERIFY: bbs3d/src/gpu_bbs3d/bbs3d.cu:7-21]:

| Member | Default |
|---|---|
| `v_rate_` | `2.0f` |
| `branch_copy_size_` | `10000` |
| `score_threshold_percentage_` | `0.0` |
| `use_timeout_` | `false` |
| `timeout_duration_` | `10000` ms |
| `min_rpy_` | `(-0.02f, -0.02f, 0.0f)` |
| `max_rpy_` | `(0.02f, 0.02f, 2π)` |
| `stream` | Created by `cudaStreamCreate` in constructor [VERIFY: bbs3d/src/gpu_bbs3d/bbs3d.cu:16] |

Note: The GPU class has **no `num_threads_` knob** — parallelism is managed by
the CUDA kernel's grid (`block_size = 32`, `num_blocks =
ceil(N/32)`) [VERIFY: bbs3d/src/gpu_bbs3d/calc_score.cu:70-71].

---

## 6. Spatial Hash Function (used everywhere)

```cpp
const std::uint32_t hash =
    (coord[0] * 73856093) ^ (coord[1] * 19349669) ^ (coord[2] * 83492791);
```
Same triple of primes used at:

- CPU build: [VERIFY: bbs3d/src/cpu_bbs3d/voxelmaps.cpp:76]
- CPU lookup: [VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:150]
- GPU build: [VERIFY: bbs3d/src/gpu_bbs3d/voxelmaps.cu:78]
- GPU lookup: [VERIFY: bbs3d/src/gpu_bbs3d/calc_score.cu:38]

This is the well-known **Teschner (2003) spatial hash**.
The primes `73856093`, `19349669`, `83492791` are chosen so that the XOR of
their products produces low collision rate on integer voxel coords.

`hash` is `uint32_t` and the bucket index is `(hash + j) % num_buckets` with `j
∈ [0, max_bucket_scan_count)` — **linear open addressing** with a bounded probe
budget.

---

## 7. `pciof::UnorderedVoxelMap<T>` (filter-side)

Different from `VoxelMaps`! This one lives in
`pointcloud_iof/filter.hpp` and stores **a list of points per voxel** rather
than an occupancy flag:

```cpp
template <typename T>
using UnorderedVoxelMap =
  std::unordered_map<Eigen::Vector3i, std::vector<Vector3<T>>, VectorHash<T>, VctorEqual<T>>;
```
[VERIFY: bbs3d/include/pointcloud_iof/filter.hpp:29]

It backs `pciof::filter(points, voxel_width)` — a centroid-style voxel-grid
downsampler [VERIFY: bbs3d/include/pointcloud_iof/filter.hpp:43-71]. To avoid
catastrophic cancellation when averaging large coordinates, the implementation
subtracts the first point as an offset before accumulating
[VERIFY: bbs3d/include/pointcloud_iof/filter.hpp:55-67].

---

## 8. Relationship Diagram

```
              ┌──────────────────┐
              │   BBS3D (CPU)    │
              │   BBS3D (GPU)    │
              └──────────────────┘
              owns │            uses (transient)
                   ▼                       ▼
        ┌──────────────────┐   ┌──────────────────────────┐
        │  VoxelMaps       │   │ std::priority_queue<     │
        │   multi_buckets_ │   │   DiscreteTransformation>│
        │   voxelmaps_res_ │   │ (max-heap by score)      │
        │   (+ device dups │   └──────────────────────────┘
        │    on GPU)       │
        └──────┬───────────┘
               │ vector of
               ▼
        ┌──────────────────────┐
        │ Buckets =            │   per level L:
        │   std::vector<       │   num_buckets ≈ 2^k × |occupied|
        │     Eigen::Vector4i> │   slot = (x, y, z, occupied)
        └──────────────────────┘

  Per BnB iteration:
        DiscreteTransformation
            │  level → VoxelMapInfo[level], buckets[level], AngularInfo[level]
            ▼
        score(trans) := Σ_{p ∈ src_points} 1[(trans · p) hashes to occupied bucket]
```

---

## 9. Verification Checklist

- [x] Every struct/class definition opened and read in full
- [x] All public fields listed with `[VERIFY:]` tag to declaration site
- [x] Default values traced to constructors (no "assumed" defaults)
- [x] CPU/GPU twin structures cross-compared
- [x] Hash function verified identical across all 4 call sites
- [x] Neighbour expansion (7 entries) cross-checked between CPU and GPU
