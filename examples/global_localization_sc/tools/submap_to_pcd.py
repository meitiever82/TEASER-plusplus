#!/usr/bin/env python3
"""Convert glim_ros submap directories to PCD files for use as SC + TEASER test scans.

Submap layout:
  <session>/<NNNNNN>/
    points_compact.bin       N * 3 * float32, in submap-local frame (origin = submap pose)
    intensities_compact.bin  N * 1 * float32, per-point lidar intensity (optional;
                                              if absent, PCDs are written as xyz-only)
    covs_compact.bin         N * 6 * float32, per-point cov upper triangle
    data.txt                 metadata, contains T_world_origin (4x4)

Usage:
  submap_to_pcd.py SESSION_DIR --submap 0          # one submap
  submap_to_pcd.py SESSION_DIR --all               # all submaps under SESSION_DIR
  submap_to_pcd.py SESSION_DIR --all --in-world    # also output a copy transformed to world frame

Outputs go alongside each submap as:
  submap_local.pcd     points as stored (submap-local frame)
  submap_world.pcd     (only with --in-world) points transformed via T_world_origin
  gt_pose.txt          x y z roll pitch yaw  (radians, world frame, T_world_origin)
"""
import argparse
import os
import re
import struct
import sys
from pathlib import Path

import numpy as np


def parse_matrix4_from_lines(lines, label):
    """Extract the 4x4 matrix following a line that starts with `label:`."""
    for i, line in enumerate(lines):
        if line.strip().startswith(label):
            rows = []
            for j in range(1, 5):
                vals = re.findall(r"-?\d+\.?\d*(?:[eE][+-]?\d+)?", lines[i + j])
                if len(vals) < 4:
                    raise ValueError(f"bad matrix row near {label}: {lines[i + j]!r}")
                rows.append([float(v) for v in vals[:4]])
            return np.array(rows, dtype=np.float64)
    raise KeyError(label)


def rot_to_rpy(R):
    """ZYX (yaw-pitch-roll) extraction, matching the local_refinement convention."""
    pitch = np.arcsin(-R[2, 0])
    roll = np.arctan2(R[2, 1], R[2, 2])
    yaw = np.arctan2(R[1, 0], R[0, 0])
    return roll, pitch, yaw


def write_pcd_binary(path, pts, intensities=None):
    """Write a minimal binary PCD.

    If `intensities` is None, writes FIELDS x y z (12 bytes/point).
    If provided (1D array of length n), writes FIELDS x y z intensity (16 bytes/point).
    """
    n = len(pts)
    if intensities is None:
        header = (
            "# .PCD v0.7 - Point Cloud Data file format\n"
            "VERSION 0.7\n"
            "FIELDS x y z\n"
            "SIZE 4 4 4\n"
            "TYPE F F F\n"
            "COUNT 1 1 1\n"
            f"WIDTH {n}\n"
            "HEIGHT 1\n"
            "VIEWPOINT 0 0 0 1 0 0 0\n"
            f"POINTS {n}\n"
            "DATA binary\n"
        )
        payload = pts.astype(np.float32)
    else:
        if len(intensities) != n:
            raise ValueError(f"len(intensities)={len(intensities)} != len(pts)={n}")
        header = (
            "# .PCD v0.7 - Point Cloud Data file format\n"
            "VERSION 0.7\n"
            "FIELDS x y z intensity\n"
            "SIZE 4 4 4 4\n"
            "TYPE F F F F\n"
            "COUNT 1 1 1 1\n"
            f"WIDTH {n}\n"
            "HEIGHT 1\n"
            "VIEWPOINT 0 0 0 1 0 0 0\n"
            f"POINTS {n}\n"
            "DATA binary\n"
        )
        # Interleave xyz + i into a (n, 4) float32 array
        payload = np.empty((n, 4), dtype=np.float32)
        payload[:, :3] = pts.astype(np.float32)
        payload[:, 3] = intensities.astype(np.float32)
    with open(path, "wb") as f:
        f.write(header.encode("ascii"))
        payload.tofile(f)


def process_submap(submap_dir: Path, in_world: bool):
    points_bin = submap_dir / "points_compact.bin"
    intensities_bin = submap_dir / "intensities_compact.bin"
    data_txt = submap_dir / "data.txt"
    if not points_bin.exists():
        print(f"[skip] {submap_dir}: no points_compact.bin")
        return
    if not data_txt.exists():
        print(f"[skip] {submap_dir}: no data.txt")
        return

    pts = np.fromfile(points_bin, dtype=np.float32).reshape(-1, 3)
    intensities = None
    if intensities_bin.exists():
        intensities = np.fromfile(intensities_bin, dtype=np.float32)
        if len(intensities) != len(pts):
            print(f"[warn] {submap_dir}: intensity len {len(intensities)} != pts len {len(pts)} — dropping intensity")
            intensities = None
    with data_txt.open() as f:
        lines = f.readlines()
    T = parse_matrix4_from_lines(lines, "T_world_origin")
    R = T[:3, :3]
    t = T[:3, 3]
    roll, pitch, yaw = rot_to_rpy(R)

    # Always write the local-frame cloud.
    local_path = submap_dir / "submap_local.pcd"
    write_pcd_binary(local_path, pts, intensities)

    # ALSO write a "gravity-levelled" cloud: origin stays at scan center, but the cloud is
    # rotated so its z-axis aligns with world z (gravity). This is what SC actually wants:
    # scanner-centric + z = gravity, with yaw still unknown. In production this comes from
    # the IMU; here we synthesise it from GT by undoing R's roll and pitch but NOT its yaw.
    #
    # Decomposition: R = R_yaw(yaw) * R_pitch(pitch) * R_roll(roll).
    # We apply R_level = (R_pitch * R_roll)^-1 to local points; the residual cloud-to-world
    # rotation is just R_yaw, exactly the "perfectly-levelled scanner facing some yaw" state.
    Ry = np.array([
        [np.cos(pitch),  0, np.sin(pitch)],
        [0,              1, 0            ],
        [-np.sin(pitch), 0, np.cos(pitch)],
    ])
    Rx = np.array([
        [1, 0,              0            ],
        [0, np.cos(roll), -np.sin(roll)],
        [0, np.sin(roll),  np.cos(roll)],
    ])
    # We want: R_world @ P_local = R_yaw @ P_levelled, so P_levelled = R_pr @ P_local
    # (NOT R_pr.T - applying the inverse would also rotate the cloud by the inverse, leaving
    # the full rotation as R_yaw @ R_pr^2, with z not aligned to gravity).
    R_pr = Ry @ Rx
    R_level = R_pr
    pts_levelled = (R_level @ pts.T).T
    levelled_path = submap_dir / "submap_levelled.pcd"
    write_pcd_binary(levelled_path, pts_levelled.astype(np.float32), intensities)

    # Optionally write a world-frame cloud (useful for visualizing alignment against the global map).
    if in_world:
        world_pts = (R @ pts.T).T + t
        world_path = submap_dir / "submap_world.pcd"
        write_pcd_binary(world_path, world_pts, intensities)

    gt_path = submap_dir / "gt_pose.txt"
    with gt_path.open("w") as f:
        f.write(f"# T_world_origin as x y z roll pitch yaw (radians)\n")
        f.write(f"{t[0]} {t[1]} {t[2]} {roll} {pitch} {yaw}\n")
        f.write("# 4x4:\n")
        for row in T:
            f.write(" ".join(f"{v:+.9f}" for v in row) + "\n")

    int_tag = ""
    if intensities is not None:
        int_tag = f" i[{intensities.min():.1f},{intensities.max():.1f}]"
    print(
        f"{submap_dir.name}: n={len(pts)}{int_tag}, "
        f"pos=({t[0]:+.2f},{t[1]:+.2f},{t[2]:+.2f}), "
        f"rpy=({np.degrees(roll):+.1f},{np.degrees(pitch):+.1f},{np.degrees(yaw):+.1f})deg"
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("session_dir", type=Path)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--submap", type=int, help="single submap id (e.g. 0)")
    g.add_argument("--all", action="store_true")
    ap.add_argument("--in-world", action="store_true",
                    help="also write a copy transformed to world frame (for visualization)")
    args = ap.parse_args()

    if not args.session_dir.exists():
        print(f"no such dir: {args.session_dir}", file=sys.stderr)
        sys.exit(1)

    if args.submap is not None:
        sub = args.session_dir / f"{args.submap:06d}"
        process_submap(sub, args.in_world)
    else:
        for sub in sorted(args.session_dir.iterdir()):
            if sub.is_dir() and re.match(r"\d{6}", sub.name):
                process_submap(sub, args.in_world)


if __name__ == "__main__":
    main()
