# 3D-BBS — Data Flow (数据流分析)

> Tracks every byte of point-cloud / pose data from raw PCD on disk to the
> 4×4 transform returned by `BBS3D::get_global_pose()`. All flow steps were
> derived by reading actual call sites listed below.

Companion docs: `00-SYSTEM_OVERVIEW.md`, `01-DATA_STRUCTURES.md`.

---

## 1. Data Pipeline Diagram

```
                       ┌──────────────────────────────┐
                       │  target PCDs    source PCDs  │   (on disk)
                       └──────┬──────────────┬────────┘
                              │              │
                              ▼              ▼
                   ┌─────────────────────────────────┐
                   │  pciof::load_tar_clouds /        │
                   │  load_src_points_with_filename   │  (PCL or pure)
                   └──────┬───────────────────┬───────┘
                          │  PCL clouds       │  PCL clouds
                          ▼                   ▼
                   ┌──────────────────┐  ┌──────────────────────┐
                   │ pcl_to_eigen     │  │ optional VoxelGrid    │
                   │ Vector3{d,f}     │  │ + crop scan range     │
                   └────────┬─────────┘  └──────────┬───────────┘
                            │                       │
                            ▼                       ▼
            ┌─────────────────────────┐   ┌──────────────────────────┐
            │  BBS3D::set_tar_points  │   │  BBS3D::set_src_points   │
            │  (or set_voxelmaps_     │   │  (CPU: copy host         │
            │   coords from disk)     │   │   GPU: + cudaMemcpyAsync)│
            └────────────┬────────────┘   └──────────┬───────────────┘
                         │ build pyramid             │
                         ▼                           │
            ┌──────────────────────────┐             │
            │ VoxelMaps::               │             │
            │   create_voxelmaps()      │             │
            │   (or load from disk)     │             │
            │                           │             │
            │ multi_buckets_[0..L]      │             │
            │ voxelmaps_res_[0..L]      │             │
            │ d_multi_buckets_ptrs_ (G) │             │
            │ d_voxelmaps_info_ (G)     │             │
            └────────────┬──────────────┘             │
                         │                            │
                         ▼                            ▼
                   ┌─────────────────────────────────────┐
                   │     BBS3D::localize()                │
                   │  1. calc_angular_info()              │
                   │  2. create_init_transset(level=max)  │
                   │  3. parallel scoring (OpenMP / CUDA) │
                   │  4. priority_queue BnB loop          │
                   │  5. best_trans.create_matrix()       │
                   └────────────┬─────────────────────────┘
                                │ Eigen::Matrix4{f,d}
                                ▼
                  ┌─────────────────────────────────────┐
                  │ get_global_pose() / get_best_score()│
                  │ has_localized() / has_timed_out()   │
                  └─────────────────────────────────────┘
```

---

## 2. Stage 1 — Loading Target Cloud (one-time)

### 2.1 PCL path (`test/`, `ros2_test/`)

```cpp
pciof::load_tar_clouds(tar_path, tar_leaf_size, tar_cloud_ptr);
```
[VERIFY: test/src/gpu_test.cpp:41], implementation [VERIFY: bbs3d/include/pointcloud_iof/pcd_loader.hpp:14-52].

Internal flow:
1. Iterate the folder; for each `.pcd` accumulate into one `pcl::PointCloud<PointXYZ>` [VERIFY: bbs3d/include/pointcloud_iof/pcd_loader.hpp:23-37].
2. If `tar_leaf_size > 0`, apply `pcl::ApproximateVoxelGrid` [VERIFY: bbs3d/include/pointcloud_iof/pcd_loader.hpp:40-46].
3. Convert PCL → Eigen via `pciof::pcl_to_eigen` [VERIFY: bbs3d/include/pointcloud_iof/pcl_eigen_converter.hpp:9-17].

### 2.2 PCL-free path (used by `voxelmaps_saver` and BBS3D itself when loading saved coords)

`pciof::load_tar_points<T>(folder, voxel_width, points)`
[VERIFY: bbs3d/include/pointcloud_iof/pcd_loader_without_pcl.hpp:11-38] does the
same job using the hand-rolled `read_pcd<T>` parser
([VERIFY: bbs3d/include/pointcloud_iof/pcd_io.hpp:12-127]).

The parser:
- Reads the first 11 ASCII header lines, extracting `FIELDS`, `SIZE`, `TYPE`, `COUNT`, `POINTS`.
- Computes byte offsets so it can `ignore` non-xyz fields without copying.
- Reads `(x, y, z)` as binary `T` (must match the field `SIZE`).
- Returns `std::vector<Eigen::Matrix<T,3,1>>`.

---

## 3. Stage 2 — Building the Hierarchical Voxelmap

Two entry points, both convergent on `VoxelMaps::create_voxelmaps(...)`:

### 3.1 Fresh build

`BBS3D::set_tar_points(points, min_level_res, max_level)`
- CPU: [VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:25-30]
- GPU: [VERIFY: bbs3d/src/gpu_bbs3d/bbs3d.cu:31-36]

calls `VoxelMaps::create_voxelmaps(points, v_rate)` which **iterates levels
`L = 0 .. max_level`** (= `max_level+1` total) doing:

```cpp
// CPU [VERIFY: bbs3d/src/cpu_bbs3d/voxelmaps.cpp:9-50]
for (i = 0; i < max_level+1; i++) {
  coords = quantise(points, resolution);            // floor(p / res)
  for c in coords:
    if !in_map(c): map[c] = 1;
    elif map[c] == -1: map[c] = 1;
    elif map[c] == 1: continue;
    for n in 7_neighbours(c):                       // -1, -1, 0 etc.
      if !in_map(n): map[n] = -1;
  buckets = create_hash_buckets(map);               // open-addressing
  multi_buckets_[i] = buckets;
  voxelmaps_res_[i] = resolution;
  resolution *= v_rate;                             // 2.0 default
}
```

For the GPU, after the host build [VERIFY: bbs3d/src/gpu_bbs3d/voxelmaps.cu:10-48] the buckets are uploaded:

```cpp
set_buckets_on_device(multi_buckets, v_rate, stream);
```
[VERIFY: bbs3d/src/gpu_bbs3d/voxelmaps.cu:51]

Inside `set_buckets_on_device` [VERIFY: bbs3d/src/gpu_bbs3d/voxelmaps.cu:101-144]:
- Each level's buckets are `cudaMemcpyAsync`'d to `d_multi_buckets_[i]`.
- A `VoxelMapInfo{res, inv_res, max_bucket_scan_count, num_buckets}` is built per level and copied to `d_voxelmaps_info_`.
- A small vector of raw `Eigen::Vector4i*` pointers (one per level) is copied to `d_multi_buckets_ptrs_`, which the kernel reads via `trans.level` indexing [VERIFY: bbs3d/src/gpu_bbs3d/calc_score.cu:23].
- `cudaStreamSynchronize(stream)` blocks until all transfers complete [VERIFY: bbs3d/src/gpu_bbs3d/voxelmaps.cu:143].

### 3.2 Pre-built load

`BBS3D::set_voxelmaps_coords(folder_path)` — fast path that skips construction.
- CPU: [VERIFY: bbs3d/src/cpu_bbs3d/voxelmaps_io.cpp:9-19]
- GPU: [VERIFY: bbs3d/src/gpu_bbs3d/voxelmaps_io.cu:8-22]

Reads `<folder>/voxelmaps_coords/voxel_params.txt` (key/value pairs:
`min_level_res`, `max_level`, `v_rate`) and per-level `.pcd` files containing
the integer voxel coords:

```cpp
voxelmaps_ptr_->multi_buckets_[i] = coords4i;                       // CPU
voxelmaps_ptr_->voxelmaps_res_[i] = min_level_res * v_rate^i;
```
[VERIFY: bbs3d/src/cpu_bbs3d/voxelmaps_io.cpp:50-56]

> ⚠️ **CPU pre-built path skips the open-addressing rebuild** — `coords4i` is
> used as-is, i.e. the buckets aren't re-hashed/probed. Lookup still works
> because `calc_score`'s probe is consistent with any layout, but the original
> probe-friendly distribution from `create_hash_buckets` is *not* preserved.
> See [VERIFY: bbs3d/src/cpu_bbs3d/voxelmaps_io.cpp:45-56]; cross-ref to GPU
> equivalent [VERIFY: bbs3d/src/gpu_bbs3d/voxelmaps_io.cu:43-58] which then
> still routes through `set_buckets_on_device` (which only does device copy, no
> re-probing either). This is intentional but worth noting: **pre-built
> coordinates are searched linearly within the on-disk ordering**, so the
> success rate is whatever the saver achieved when it built the map.

Both readers also derive `init_t{x,y,z}_range_` from the bbox of the
top-level bucket set [VERIFY: bbs3d/src/cpu_bbs3d/voxelmaps_io.cpp:59-71].

---

## 4. Stage 3 — Loading a Source Frame

Per-frame, callers feed `set_src_points`:

```cpp
// CPU [VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:32-37]
src_points_.clear(); src_points_.shrink_to_fit();
src_points_.resize(points.size());
std::copy(points.begin(), points.end(), src_points_.begin());

// GPU [VERIFY: bbs3d/src/gpu_bbs3d/bbs3d.cu:38-54]
src_points_  = copy(host);
d_src_points_.resize(...);
cudaMemcpyAsync(d_src_points_.data(), points.data(), N*sizeof(float3),
                cudaMemcpyHostToDevice, stream);
```

If the caller hasn't already set the translation range, they typically call
`set_trans_search_range(tar_points)` *once* against the target cloud
[VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:39-50]; on every per-frame call, the
search range stays fixed.

---

## 5. Stage 4 — `localize()` (the BnB loop)

The main flow inside `localize()` (per-frame):

| Step | CPU line | GPU line | Effect |
|---|---|---|---|
| 1. start clock | [VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:172] | [VERIFY: bbs3d/src/gpu_bbs3d/bbs3d.cu:152] | `start_time = now()` |
| 2. score_threshold = floor(N_src * pct) | [VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:177] | [VERIFY: bbs3d/src/gpu_bbs3d/bbs3d.cu:158] | early-termination floor |
| 3. compute `AngularInfo[0..L]` | [VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:184-185] | [VERIFY: bbs3d/src/gpu_bbs3d/bbs3d.cu:164-165] | per-level angular grid |
| 3b. (GPU only) upload `d_ang_info_vec` | n/a | [VERIFY: bbs3d/src/gpu_bbs3d/bbs3d.cu:166-172] | one `cudaMemcpyAsync` |
| 4. build initial transset | [VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:186] | [VERIFY: bbs3d/src/gpu_bbs3d/bbs3d.cu:174] | 6-nested loop at top level |
| 5. score initial transset | OpenMP `parallel for` over `calc_score` [VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:193-196] | single batched kernel via `calc_scores(...)` [VERIFY: bbs3d/src/gpu_bbs3d/bbs3d.cu:177] | every entry gets a `score` |
| 6. seed `priority_queue` | [VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:198] | [VERIFY: bbs3d/src/gpu_bbs3d/bbs3d.cu:179] | max-heap |
| 6b. (GPU only) allocate `branch_stock` | n/a | [VERIFY: bbs3d/src/gpu_bbs3d/bbs3d.cu:181] [VERIFY: bbs3d/src/gpu_bbs3d/bbs3d.cu:215-222] | per-frame branch buffer |
| 7. BnB loop | [VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:200-241] | [VERIFY: bbs3d/src/gpu_bbs3d/bbs3d.cu:183-223] | pop, branch, score, push |
| 8. timeout check (optional) | [VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:201-204] | [VERIFY: bbs3d/src/gpu_bbs3d/bbs3d.cu:184-187] | sets `has_timed_out_` |
| 9. build `global_pose_` at level 0 | [VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:253-254] | [VERIFY: bbs3d/src/gpu_bbs3d/bbs3d.cu:235-236] | `best_trans.create_matrix(min_res, ...)` |
| 10. set `elapsed_time_`, flags | [VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:244-258] | [VERIFY: bbs3d/src/gpu_bbs3d/bbs3d.cu:225-240] | finalise outputs |

Internal data exchanged in step 7:

- **CPU per-iter:** `auto trans = trans_queue.top(); trans_queue.pop();` →
  `children = trans.branch(...)` → score children in OpenMP → push children
  with `score ≥ best_score` back into the queue
  [VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:206-239].
- **GPU per-iter:** branches are accumulated into a host-side `branch_stock_`
  vector. When `branch_stock_.size() >= branch_copy_size_` (default 10000) or
  the queue is about to empty, the stock is scored with `calc_scores(stock)`
  and the survivors pushed back into the queue
  [VERIFY: bbs3d/src/gpu_bbs3d/bbs3d.cu:192-222].

### 5.1 Why GPU buffers branches

The kernel launch overhead and the H↔D transfer cost dominate small batches.
By collecting children until either:
- `branch_stock_.size() >= 10000`, or
- the queue would otherwise become empty,

the algorithm amortises the dispatch cost over thousands of poses per launch.
The CPU version does not need this — OpenMP can profitably parallelise even
small children sets.

---

## 6. Stage 5 — Single Score Evaluation

Both backends compute the same logical function:

```
score(T) := |{ p ∈ src_points : (T · p) quantised at res hits an occupied bucket }|
```

### 6.1 CPU [VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:134-167]

```cpp
transform = trans.create_matrix(trans_res, rpy_res, min_rpy);
for (i = 0; i < src_points.size(); i++) {
  transed = transform * src_points[i];
  coord   = floor(transed * inv_res).cast<int>();
  hash    = (coord.x * 73856093) ^ (coord.y * 19349669) ^ (coord.z * 83492791);
  for (j = 0; j < max_bucket_scan_count; j++) {
    bucket = buckets[(hash + j) % num_buckets];
    if (bucket.w() == 0) break;                           // empty → done probing
    if (bucket.xyz() != coord) continue;                  // hash collision
    trans.score++; break;                                 // hit
  }
}
```

### 6.2 GPU [VERIFY: bbs3d/src/gpu_bbs3d/calc_score.cu:7-56]

Each CUDA thread evaluates one pose:

```cpp
__global__ void calc_scores_kernel(...) {
  pose_index = threadIdx.x + blockIdx.x * blockDim.x;
  if (pose_index > index_size) return;          // index_size == transset_size - 1

  trans         = trans_ptr[pose_index];
  voxelmap_info = voxelmap_info_ptr[trans.level];   // res, inv_res, num_buckets, ...
  ang_info      = d_ang_info_vec_ptr[trans.level];
  buckets       = multi_buckets_ptrs[trans.level];

  for (i = 0; i < num_points; i++) {
    p           = points_ptr[i];
    translation = trans.{x,y,z} * voxelmap_info.res;
    rotation    = AngleAxis(yaw*rpy_res.z+min.z, Z)
                * AngleAxis(pitch*rpy_res.y+min.y, Y)
                * AngleAxis(roll*rpy_res.x+min.x, X);
    transed     = rotation * p + translation;
    coord       = floor(transed * voxelmap_info.inv_res);
    hash        = teschner(coord);
    for (j = 0; j < voxelmap_info.max_bucket_scan_count; j++) {
      bucket = buckets[(hash + j) % voxelmap_info.num_buckets];
      if (bucket.xyz() != coord) continue;
      if (bucket.w() == 1) { score++; break; }
    }
  }
  trans.score = score;
}
```

Two divergences from the CPU version (confirmed line-by-line):

1. **Empty-slot handling.** CPU `break`s on `w==0` *before* the coord check;
   GPU just `continue`s. As noted in §6 of `00-SYSTEM_OVERVIEW.md`, an empty
   `(0,0,0,0)` slot cannot match a non-origin point coord, so correctness is
   preserved either way. The CPU version (added in commit `3feeede`) saves the
   remaining probe iterations once an empty slot is hit
   [VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:156-158] vs. [VERIFY: bbs3d/src/gpu_bbs3d/calc_score.cu:41-53].

2. **`pose_index > index_size` guard.** The GPU passes `transset_size - 1` as
   `index_size`, then bounds the thread with `>` instead of `>=`. Effectively
   the kernel evaluates threads `0..index_size` inclusive, i.e. all
   `transset_size` threads. Threads beyond `transset_size - 1` (i.e. the tail
   of the last 32-wide block) early-exit
   [VERIFY: bbs3d/src/gpu_bbs3d/calc_score.cu:15-18] [VERIFY: bbs3d/src/gpu_bbs3d/calc_score.cu:78].

### 6.3 Branch in the GPU dispatcher

`gpu::BBS3D::calc_scores(...)` is the host-side wrapper that:
1. Allocates a `thrust::device_vector<DiscreteTransformation<float>>` sized to the input vector.
2. `cudaMemcpyAsync` host → device.
3. Launches `calc_scores_kernel<<<num_blocks, 32, 0, stream>>>(...)`.
4. `cudaMemcpyAsync` device → host on the **same** vector (now populated with `score`).
5. `cudaStreamSynchronize(stream)` and returns.

[VERIFY: bbs3d/src/gpu_bbs3d/calc_score.cu:58-92]

Notice that `DiscreteTransformation<float>` is **bit-for-bit memcpy-safe** — it
is a POD of 8 `int`s, which is why a straight `cudaMemcpyAsync` of the host
vector works.

---

## 7. Stage 6 — Output

After the BnB loop terminates (queue empty or timeout):

```cpp
if (best_score == score_threshold || has_timed_out_) {
  has_localized_ = false;
  return;
}
global_pose_ = best_trans.create_matrix(min_res,
                                        ang_info_vec[0].rpy_res,
                                        ang_info_vec[0].min_rpy);
best_score_     = best_score;
has_localized_  = true;
```
[VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:248-258] / [VERIFY: bbs3d/src/gpu_bbs3d/bbs3d.cu:230-240]

Note `best_trans.create_matrix(...)` is invoked with `min_res` (level 0
resolution) and the **level 0** angular info — i.e. the finest discretisation,
which is correct only because the BnB invariant guarantees the best discrete
node is a leaf (`level == 0`).

The caller pulls results via:

| API | Returns |
|---|---|
| `get_global_pose()` | `Eigen::Matrix4{d,f}` — 4×4 transform (body→world) |
| `get_best_score()`  | `int` — points that matched at the leaf level |
| `get_best_score_percentage()` | `score / num_src_points` |
| `get_elapsed_time()` | `double` — `localize()` wall time in ms |
| `has_localized()`   | `bool` |
| `has_timed_out()`   | `bool` |

---

## 8. Optional Side Channels

### 8.1 Gravity alignment (ROS 2 wrappers only)

The ROS 2 wrappers consume `sensor_msgs/Imu` and pre-rotate the source cloud
*before* `localize()`:

```cpp
imu_msg = imu_buffer[nearest_index(time)];
acc     = {ax, ay, az};                    // [VERIFY: ros2_test/iridescence/src/gpu_bbs3d_iridescence/gpu_ros2_test_iridescence.cpp:132-134]
pcl::transformPointCloud(*src_cloud, *src_cloud,
                         pciof::calc_gravity_alignment_matrix(acc.cast<float>()));
```
[VERIFY: ros2_test/iridescence/src/gpu_bbs3d_iridescence/gpu_ros2_test_iridescence.cpp:135]

`calc_gravity_alignment_matrix(acc)` builds an `Rx(atan2(ay, az)) · Ry(atan2(-ax, sqrt(ay²+az²)))`
rotation and embeds it in a 4×4 identity
[VERIFY: bbs3d/include/pointcloud_iof/gravity_alignment.hpp:5-18].

The standalone `test/` executables do **not** apply this — they assume the PCDs
are already gravity-aligned (the upstream `ros2bag_to_pcd` tool referenced in
`test/test_code.md:62-63` is expected to do that conversion).

### 8.2 Pre-built voxelmaps

`test/src/voxelmaps_saver.cpp` runs only the build phase, then calls
`save_voxel_params` + `save_voxelmaps_pcd` to write
`voxelmaps_coords/voxel_params.txt` and `<level>.pcd` per level
[VERIFY: test/src/voxelmaps_saver.cpp:7-31] /
[VERIFY: bbs3d/src/cpu_bbs3d/voxelmaps_io.cpp:76-135].

On subsequent runs, `set_voxelmaps_coords(folder)` is the first thing the
test driver tries; only on failure does it fall back to `set_tar_points()`
[VERIFY: test/src/gpu_test.cpp:76-81].

---

## 9. End-to-End Per-Frame Latency Budget (from `test/src/gpu_test.cpp`)

The flow inside the source-cloud loop is:

```
for each src PCD:
  1. pcl_to_eigen(...)                     // O(N_src) copy
  2. bbs3d.set_src_points(src_points)      // O(N_src) host + cudaMemcpyAsync
  3. bbs3d.localize()                      // dominates
  4. (optional) gicp.align(...)            // refinement, not part of 3D-BBS
  5. transformPointCloud + savePCDFile     // I/O for viewing
```
[VERIFY: test/src/gpu_test.cpp:96-137]

Only `localize()` is measured by `get_elapsed_time()` — the surrounding I/O is
not.

---

## 10. Verification Checklist

- [x] All stages traced top-to-bottom through actual call sites
- [x] Both CPU and GPU paths covered in parallel
- [x] Differences between the two backends explicitly annotated
- [x] PCL-free path covered (`pcd_loader_without_pcl.hpp`)
- [x] Pre-built voxelmap fast path covered (`voxelmaps_io.{cpp,cu}`)
- [x] ROS 2 IMU side-channel (gravity alignment) included
- [x] No claims without a code reference
