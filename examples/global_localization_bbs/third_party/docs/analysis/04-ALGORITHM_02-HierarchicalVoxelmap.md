# 3D-BBS — Algorithm 2: Hierarchical Hash Voxelmap

> The data structure that makes BnB tractable: a multi-resolution pyramid of
> hash-bucket sets, each built so its score gives a valid **upper bound** on
> all finer-level scores.
>
> Companions: `01-DATA_STRUCTURES.md` (struct details),
> `03-ALGORITHM_01-BranchAndBound.md` (why the bound matters).

---

## 1. Goals

1. Constant-time `point → occupancy?` lookup at every resolution.
2. Memory linear in *occupied* voxel count (not total grid volume).
3. Occupancy preserved upward through the levels (so coarse score upper-bounds fine score).
4. Tolerant of sub-voxel discretisation noise.

The implementation reaches all four goals with a single mechanism: a
**level-specific re-quantisation** of the original target cloud, plus an
**open-addressing hash bucket** with **7-neighbour dilation**.

---

## 2. Construction Algorithm — `create_voxelmaps`

CPU [VERIFY: bbs3d/src/cpu_bbs3d/voxelmaps.cpp:9-50] and GPU
[VERIFY: bbs3d/src/gpu_bbs3d/voxelmaps.cu:10-52] versions are
character-for-character equivalent up to scalar type.

```cpp
void create_voxelmaps(const std::vector<Vector3<T>>& points, const int v_rate) {
  const int N = max_level_ + 1;
  multi_buckets_.resize(N);
  voxelmaps_res_.resize(N);

  T resolution = min_level_res_;
  for (int i = 0; i < N; i++) {
    UnorderedVoxelMap unordered_voxelmap;

    // Step 1: quantise *original* points at this level's resolution.
    std::vector<Eigen::Vector3i> coords(points.size());
    std::transform(points.begin(), points.end(), coords.begin(),
                   [&](const auto& p) {
                     return (p.array() / resolution).floor().template cast<int>();
                   });

    // Step 2: insert + 7-neighbour dilation, tri-state flag {-1, 1}.
    for (const auto& coord : coords) {
      if (unordered_voxelmap.count(coord) == 0)      unordered_voxelmap[coord] = 1;
      else if (unordered_voxelmap[coord] == -1)      unordered_voxelmap[coord] = 1;
      else if (unordered_voxelmap[coord] == 1)       continue;
      for (const auto& n : create_neighbor_coords(coord)) {
        if (unordered_voxelmap.count(n) == 0)        unordered_voxelmap[n] = -1;
      }
    }

    // Step 3: hash the map into open-addressing buckets.
    multi_buckets_[i]  = create_hash_buckets(unordered_voxelmap);
    voxelmaps_res_[i]  = resolution;
    resolution        *= v_rate;          // 2.0 by default
  }
}
```

[VERIFY: bbs3d/src/cpu_bbs3d/voxelmaps.cpp:9-50]

Two implications worth highlighting:

- **Each level re-quantises from the original points.** It does *not* downsample
  level *i−1* to produce level *i*. This guarantees the upward-monotonicity
  property required for the BnB bound (every fine-level occupied voxel has a
  coarse-level ancestor that is also marked).
- **Both `1` and `-1` flags survive** into the hash buckets. The hash-bucket
  builder treats them identically (`coord << x, y, z, 1`) — i.e. once a voxel
  enters the buckets, the algorithm forgets whether it was a "real" hit or
  just a neighbour. This is what gives 3D-BBS its sub-voxel robustness.

---

## 3. The 7-Neighbour Dilation

For an occupied coord `(x, y, z)`, exactly seven neighbours are added:

```
(-1, -1,  0)   (-1,  0,  0)   ( 0, -1,  0)
(-1, -1, -1)   (-1,  0, -1)   ( 0, -1, -1)
( 0,  0, -1)
```
[VERIFY: bbs3d/src/cpu_bbs3d/voxelmaps.cpp:52-63] / [VERIFY: bbs3d/src/gpu_bbs3d/voxelmaps.cu:54-65]

Geometrically these are *exactly* the 7 voxels of the unit cube with one or
more coordinates decremented — i.e. the 8-neighbourhood `{0, -1}³` minus the
voxel itself.

### 3.1 Why these specific 7?

Look at how scoring quantises a query point:

```cpp
const Eigen::Vector3i coord = (transed_point.array() * inv_res).floor().cast<int>();
```
[VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:149], [VERIFY: bbs3d/src/gpu_bbs3d/calc_score.cu:37]

`floor` quantisation maps a continuous point to the **negative-most lattice
corner** of its containing voxel. So a small noise of size `ε > 0` in any axis
can only push the query into the voxel at the **negative direction** from the
true one (if the true point sits near the negative wall, an `−ε` offset moves
it into the neighbour at index `−1`).

The dilation deposits the source voxel into the seven cells that the
floor-quantiser would reach for *negative* perturbations along any subset of
axes — exactly compensating for that asymmetry. The "positive-side"
neighbours are unnecessary because `floor` never crosses them given small
perturbations.

This is asymmetric on purpose; a full 26-neighbour or 6-face dilation would
oversize the occupancy and create false-positive bound estimates.

### 3.2 Effect on occupancy mass

Let `|M_L^raw|` denote the count of *strictly occupied* voxels (flag 1 only),
and `|M_L^dil|` the dilated count (flags 1 + −1). In a worst case where each
occupied voxel sits in isolation, `|M_L^dil| ≤ 8 · |M_L^raw|` (one own + 7
neighbours). In practice clustered occupancy makes the ratio much smaller
(neighbours collide with other occupied voxels and merge).

---

## 4. Hash Bucket Construction — `create_hash_buckets`

```cpp
Buckets create_hash_buckets(const UnorderedVoxelMap& m) {
  std::vector<Eigen::Vector4i> buckets;
  for (int N = m.size(); N <= m.size() * 16; N *= 2) {
    buckets.assign(N, Vector4i::Zero());
    int success = 0;
    for (const auto& v : m) {
      Vector4i coord;  coord << v.first.x(), v.first.y(), v.first.z(), 1;
      uint32_t hash = (coord[0] * 73856093) ^ (coord[1] * 19349669)
                                            ^ (coord[2] * 83492791);
      for (int j = 0; j < max_bucket_scan_count_; j++) {
        uint32_t bi = (hash + j) % N;
        if (buckets[bi].w() == 0) { buckets[bi] = coord; success++; break; }
      }
    }
    if ((double)success / m.size() > 0.999) break;
  }
  return buckets;
}
```
[VERIFY: bbs3d/src/cpu_bbs3d/voxelmaps.cpp:65-97] / [VERIFY: bbs3d/src/gpu_bbs3d/voxelmaps.cu:67-99]

Properties:

| Property | Value |
|---|---|
| Open-addressing scheme | Linear probing |
| Probe limit | `max_bucket_scan_count_` (default `10`) [VERIFY: bbs3d/src/cpu_bbs3d/voxelmaps.cpp:5] |
| Empty-slot sentinel | `Vector4i::Zero()` ⇒ `w() == 0` |
| Size growth | Powers of 2, from `|m|` up to `16·|m|` |
| Acceptance threshold | `success rate > 99.9%` |

If even at `16·|m|` the probe budget is exhausted on too many keys, the loop
exits with the last (largest) attempt — meaning *some* voxels may be silently
dropped. This is rare with the `max_bucket_scan_count = 10` budget against a
load factor `≤ 1/16`, but worth knowing for very pathological inputs.

The deliberately *small* `max_bucket_scan_count` (10) is chosen so that the
**lookup-side** cost is bounded. Every score evaluation probes at most 10
buckets per point, which translates directly to a fixed inner loop on the GPU.

### 4.1 Bucket layout = `Eigen::Vector4i`

A 16-byte slot storing the integer voxel coord plus a 1-bit "occupied" flag in
the 4th lane:

```cpp
coord << v.first.x(), v.first.y(), v.first.z(), 1;
```
[VERIFY: bbs3d/src/cpu_bbs3d/voxelmaps.cpp:74-75]

`Vector4i::Zero()` thus serves a dual purpose: it both indicates "empty"
(via `w==0`) **and** is guaranteed not to match any real coord with
`w==1` (since real entries carry `w=1`).

The GPU mirrors this layout, simply uploading the same `Eigen::Vector4i`
array into a `thrust::device_vector<Eigen::Vector4i>`
[VERIFY: bbs3d/src/gpu_bbs3d/voxelmaps.cu:106-112].

---

## 5. Score Lookup — `(coord, hash) → bucket`

The lookup-side code mirrors the build-side, with the same hash:

CPU [VERIFY: bbs3d/src/cpu_bbs3d/bbs3d.cpp:147-166]:
```cpp
const Eigen::Vector3d transed_point = transform * points[i];
const Eigen::Vector3i coord  = (transed_point.array() * inv_res).floor().cast<int>();
const uint32_t hash = (coord[0]*73856093) ^ (coord[1]*19349669) ^ (coord[2]*83492791);

for (int j = 0; j < max_bucket_scan_count; j++) {
  uint32_t bi = (hash + j) % num_buckets;
  const Vector4i& b = buckets[bi];
  if (b.w() == 0) break;                                    // empty → done
  if (b.x() != coord.x() ||
      b.y() != coord.y() ||
      b.z() != coord.z()) continue;                         // collision
  trans.score++; break;                                     // hit
}
```

GPU [VERIFY: bbs3d/src/gpu_bbs3d/calc_score.cu:36-53]:
```cpp
const Vector3i coord = (transed_point.array() * inv_res).floor().cast<int>();
const uint32_t hash  = teschner(coord);
for (int j = 0; j < max_bucket_scan_count; j++) {
  uint32_t bi = (hash + j) % num_buckets;
  const Vector4i b = buckets[bi];
  if (b.x() != coord.x() ||
      b.y() != coord.y() ||
      b.z() != coord.z()) continue;
  if (b.w() == 1) { score++; break; }                       // hit
}
```

Key differences:
- CPU breaks early on `w==0`. GPU does not (an empty `(0,0,0,0)` bucket can't
  match any non-origin coord, so the loop simply runs to its probe budget).
  Per the git log this CPU optimisation was added in commit `3feeede`
  ("Optimize hash lookup by breaking early on empty buckets for cpu bbs3d
  (#54)").
- GPU checks `w == 1` *after* the coord match, mostly for symmetry with the
  build. CPU implicitly relies on `w!=0 ⇒ w==1` (the only nonzero w value the
  builder ever writes).

---

## 6. Memory Footprint

Per level `L` with `|M_L^dil|` dilated voxels and `N_buckets_L` slots:

- Bucket array: `16 · N_buckets_L` bytes (one `Vector4i` per slot).
- `N_buckets_L` grows as a power of two until success rate > 99.9%, so
  `N_buckets_L ∈ [|M_L^dil|, 16·|M_L^dil|]`. In practice the very first
  power-of-two value usually suffices [VERIFY: bbs3d/src/cpu_bbs3d/voxelmaps.cpp:68].
- Across `L_max+1` levels, the total scales roughly with
  `Σ_L |M_L^dil|`. Since `|M_L^dil|` shrinks as `L` grows (coarser voxels
  cluster more), the pyramid is mostly bounded by `L=0`.

Per the README, total construction time is **~3494 ms** on a Core i7 + RTX
2060 [VERIFY: README.md:55], with subsequent runs reducing to **130 ms** via
the disk cache [VERIFY: README.md:56]:

```
Hierarchical voxelmap construction
  Paper:               9272 ms
  Latest implementation:    3494 ms
  Load saved voxelmaps directly: 130 ms
```

---

## 7. Persistence (Disk Cache)

Saved by `save_voxel_params` + `save_voxelmaps_pcd`:

```
<target_folder>/voxelmaps_coords/
    voxel_params.txt        # text key-value
    0.pcd                   # bucket coords at level 0
    1.pcd                   # ...
    ...
    <max_level>.pcd
```

The `.pcd` files contain *the dilated occupied voxel coords as integers*
written via the project's `pciof::save_pcd<int>(...)`
[VERIFY: bbs3d/src/cpu_bbs3d/voxelmaps_io.cpp:117-131]. Each is a
header-then-binary file with `FIELDS x y z`, `TYPE I I I`
[VERIFY: bbs3d/include/pointcloud_iof/pcd_io.hpp:138-153].

`voxel_params.txt` stores three plain-text parameters:

```
min_level_res <float>
max_level     <int>
v_rate        <float>
```
[VERIFY: bbs3d/src/cpu_bbs3d/voxelmaps_io.cpp:94-97]

These three parameters fully determine the per-level resolutions:
`voxelmaps_res_[i] = min_level_res * v_rate^i`
[VERIFY: bbs3d/src/cpu_bbs3d/voxelmaps_io.cpp:55].

> The disk cache **stores `multi_buckets_[i]` directly without re-running the
> open-addressing builder**. This is intentional: each saved entry is a
> "coord with `w=1`", and the lookup probe is consistent with whatever order
> the file produces. The downside, as noted in `02-DATA_FLOW.md` §3.2, is that
> the saved layout's probe distribution is fixed at save time.

---

## 8. Per-Level Resolution Schedule

```
i = 0 → res_0  = min_level_res                        (default 1.0 m)
i = 1 → res_1  = min_level_res · v_rate               (default 2.0 m)
i = 2 → res_2  = min_level_res · v_rate²              (default 4.0 m)
...
i = L → res_L  = min_level_res · v_rate^L
```

With defaults `min_level_res = 1.0`, `v_rate = 2.0`, `max_level = 6`, the
resolution sequence is `1, 2, 4, 8, 16, 32, 64 m`. Each step doubles the
linear resolution, *increasing volume eight-fold*. So if the level-0 map has
`|M_0^dil| ≈ N`, levels above hold roughly `N/8, N/64, ...` entries — the
**pyramid is bottom-heavy** in both memory and lookup cost.

---

## 9. Worked Example

Suppose the target cloud is a 50 m × 50 m × 5 m planar room with 100,000
points. After `load_tar_clouds` with `tar_leaf_size = 0.1 m`:

- Level 0 (1 m voxels): expect roughly `50·50·5 / 1³ = 12,500` unique voxels,
  occupied to varying degrees by the 100k points (most points fall into a
  ~2,500-cell floor + walls region). Add 7-neighbour dilation: `~2,500 · 8 ≈
  20,000` dilated voxels.
- Level 1 (2 m): `~3,125` dilated voxels.
- Level 6 (64 m): essentially a single coarse cell — `O(1)` voxels.

Bucket arrays per level (using `2^k` rule, smallest sufficient):

| Level | Dilated voxels | Buckets (≥) |
|---|---|---|
| 0 | 20,000 | 32,768 |
| 1 | 3,125  |  4,096 |
| 2 |   400  |    512 |
| 3 |    50  |     64 |
| 4 |     8  |     16 |
| 5 |     2  |      2 |
| 6 |     1  |      1 |

Total pyramid memory: `(32768 + 4096 + ...) · 16 bytes ≈ 600 KB`. Well under
1 MB for a non-trivial room-scale map.

---

## 10. Verification Checklist

- [x] Construction algorithm read line-by-line in both CPU and GPU sources
- [x] 7-neighbour dilation enumerated and matched against the `floor` quantiser
- [x] Hash function (Teschner primes) verified identical across build + lookup
- [x] Open-addressing probe budget = `max_bucket_scan_count = 10` (constructor default)
- [x] Disk-cache structure (`voxel_params.txt` + `<i>.pcd`) confirmed against the I/O code
- [x] CPU vs GPU lookup-side difference (`w==0` early break) documented
