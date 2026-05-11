// Phase 5: BBS coarse global localization + TEASER local refinement.
//
// Usage (single-submap test):
//   bbs_localize \
//       --map     /path/to/global_map.pcd \
//       --scan    /path/to/submap_levelled.pcd \
//       --gt      /path/to/data.txt        (optional, for error reporting)
//       --output  /tmp/result.txt
//
// Stage A: gpu::BBS3D runs branch-and-bound over the voxelized map with the scan
// gravity-aligned (3DOF: yaw + xy + z, roll/pitch fixed to ±0.02 rad).
//
// Stage B: refine_pipeline::refineLocalPose uses BBS's coarse 4x4 as the rough pose,
// runs FPFH + TEASER GNC_TLS to get a precise correction.
//
// Output: refined 4x4 in world frame, plus quality metrics.

#include <array>
#include <chrono>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <pcl/io/pcd_io.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

#include <teaser/ply_io.h>
#include <teaser/geometry.h>

#include <gpu_bbs3d/bbs3d.cuh>

#include "refine_pipeline.hpp"

namespace {

struct Args {
  std::string map_path;
  std::string scan_path;
  std::string gt_path;
  std::string output_path = "/tmp/bbs_result.txt";
  std::string gt_frame = "levelled";  // levelled | raw  (levelled: strip roll/pitch from GT)

  // BBS knobs
  float min_level_res = 0.5f;
  int max_level = 6;
  float score_threshold = 0.9f;
  float bbs_src_leaf = 0.2f;
  float bbs_tar_leaf = 0.1f;
  int bbs_timeout_ms = 10000;

  // Refinement knobs (subset of refine_pipeline::Config)
  std::string algo = "TEASER";
  double noise_bound = 0.15;
  double local_radius = 30.0;
  double voxel_size_map = 0.3;
  double voxel_size_scan = 0.2;
  double fpfh_normal_radius = 1.0;
  double fpfh_radius = 2.0;
  int max_corres = 3000;
};

void printUsage(const char* a0) {
  std::cout
    << "Usage: " << a0 << " --map MAP.pcd --scan SCAN.pcd [--gt GT.txt] [--output OUT.txt]\n"
    << "BBS opts:    --min-level-res 0.5  --max-level 6  --score-threshold 0.9\n"
    << "             --bbs-src-leaf 0.2  --bbs-tar-leaf 0.1  --bbs-timeout-ms 10000\n"
    << "Refine opts: --algo TEASER|Quatro  --noise-bound 0.15  --local-radius 30\n"
    << "             --voxel-map 0.3  --voxel-scan 0.2\n"
    << "             --fpfh-normal-r 1.0  --fpfh-r 2.0  --max-corres 3000\n";
}

bool parseArgs(int argc, char** argv, Args& a) {
  for (int i = 1; i < argc; ++i) {
    std::string k = argv[i];
    auto next = [&](double& v) { if (i + 1 < argc) v = std::stod(argv[++i]); };
    auto nextS = [&](std::string& v) { if (i + 1 < argc) v = argv[++i]; };
    auto nextF = [&](float& v) { if (i + 1 < argc) v = std::stof(argv[++i]); };
    auto nextI = [&](int& v) { if (i + 1 < argc) v = std::stoi(argv[++i]); };
    if      (k == "--map")               nextS(a.map_path);
    else if (k == "--scan")              nextS(a.scan_path);
    else if (k == "--gt")                nextS(a.gt_path);
    else if (k == "--gt-frame")          nextS(a.gt_frame);
    else if (k == "--output")            nextS(a.output_path);
    else if (k == "--min-level-res")     nextF(a.min_level_res);
    else if (k == "--max-level")         nextI(a.max_level);
    else if (k == "--score-threshold")   nextF(a.score_threshold);
    else if (k == "--bbs-src-leaf")      nextF(a.bbs_src_leaf);
    else if (k == "--bbs-tar-leaf")      nextF(a.bbs_tar_leaf);
    else if (k == "--bbs-timeout-ms")    nextI(a.bbs_timeout_ms);
    else if (k == "--algo")              nextS(a.algo);
    else if (k == "--noise-bound")       next(a.noise_bound);
    else if (k == "--local-radius")      next(a.local_radius);
    else if (k == "--voxel-map")         next(a.voxel_size_map);
    else if (k == "--voxel-scan")        next(a.voxel_size_scan);
    else if (k == "--fpfh-normal-r")     next(a.fpfh_normal_radius);
    else if (k == "--fpfh-r")            next(a.fpfh_radius);
    else if (k == "--max-corres")        nextI(a.max_corres);
    else if (k == "-h" || k == "--help") { printUsage(argv[0]); return false; }
    else {
      std::cerr << "Unknown arg: " << k << "\n";
      printUsage(argv[0]);
      return false;
    }
  }
  if (a.map_path.empty() || a.scan_path.empty()) {
    printUsage(argv[0]);
    return false;
  }
  return true;
}

bool endsWith(const std::string& s, const std::string& suf) {
  return s.size() >= suf.size()
      && std::equal(suf.rbegin(), suf.rend(), s.rbegin(),
                    [](char a, char b) { return std::tolower(a) == std::tolower(b); });
}

bool loadCloud(const std::string& path, teaser::PointCloud& out) {
  if (endsWith(path, ".ply")) {
    teaser::PLYReader r;
    return r.read(path, out) == 0;
  }
  if (endsWith(path, ".pcd")) {
    pcl::PointCloud<pcl::PointXYZ> tmp;
    if (pcl::io::loadPCDFile<pcl::PointXYZ>(path, tmp) < 0) return false;
    out.clear();
    out.reserve(tmp.size());
    for (const auto& p : tmp.points) {
      if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z)) continue;
      out.push_back({p.x, p.y, p.z});
    }
    return true;
  }
  std::cerr << "Unsupported extension: " << path << "\n";
  return false;
}

std::vector<Eigen::Vector3f> toEigenVec(const teaser::PointCloud& c) {
  std::vector<Eigen::Vector3f> out;
  out.reserve(c.size());
  for (const auto& p : c) out.emplace_back(p.x, p.y, p.z);
  return out;
}

// First-hit voxel downsample on a vector of Eigen::Vector3f (used for BBS src/tar).
std::vector<Eigen::Vector3f> voxelDown(const std::vector<Eigen::Vector3f>& pts, float leaf) {
  if (leaf <= 0.0f) return pts;
  struct Key { int x, y, z; bool operator==(const Key& o) const { return x == o.x && y == o.y && z == o.z; } };
  struct Hash {
    size_t operator()(const Key& k) const {
      return std::hash<long long>()(
          static_cast<long long>(k.x) * 73856093LL
          ^ static_cast<long long>(k.y) * 19349663LL
          ^ static_cast<long long>(k.z) * 83492791LL);
    }
  };
  std::unordered_map<Key, Eigen::Vector3f, Hash> seen;
  seen.reserve(pts.size());
  const float inv = 1.0f / leaf;
  for (const auto& p : pts) {
    Key k{static_cast<int>(std::floor(p.x() * inv)),
          static_cast<int>(std::floor(p.y() * inv)),
          static_cast<int>(std::floor(p.z() * inv))};
    seen.emplace(k, p);
  }
  std::vector<Eigen::Vector3f> out;
  out.reserve(seen.size());
  for (auto& kv : seen) out.push_back(kv.second);
  return out;
}

// Parse ground-truth pose. Supports two formats:
//   (A) gt_pose.txt from submap_to_pcd.py: "# ...\n<x y z roll pitch yaw>\n..."
//   (B) data.txt from glim_ros:             "id: N\nT_world_origin: \n<4 rows of 4 floats>"
bool loadGT(const std::string& path, Eigen::Matrix4d& T_gt) {
  std::ifstream f(path);
  if (!f.is_open()) return false;
  std::string line;
  std::vector<std::string> lines;
  while (std::getline(f, line)) lines.push_back(line);

  // Try (A): first non-comment line with 6 floats.
  for (const auto& l : lines) {
    if (l.empty() || l[0] == '#') continue;
    std::istringstream iss(l);
    double x, y, z, roll, pitch, yaw;
    if (iss >> x >> y >> z >> roll >> pitch >> yaw) {
      // Sanity: extra tokens after yaw mean it's a matrix row, not RPY.
      std::string tail;
      if (!(iss >> tail)) {
        Eigen::AngleAxisd rx(roll, Eigen::Vector3d::UnitX());
        Eigen::AngleAxisd ry(pitch, Eigen::Vector3d::UnitY());
        Eigen::AngleAxisd rz(yaw, Eigen::Vector3d::UnitZ());
        T_gt = Eigen::Matrix4d::Identity();
        T_gt.block<3, 3>(0, 0) = (rz * ry * rx).toRotationMatrix();
        T_gt(0, 3) = x; T_gt(1, 3) = y; T_gt(2, 3) = z;
        return true;
      }
      break;  // fall through to (B)
    }
  }

  // Try (B): collect 4 lines after "T_world_origin:" with 4 floats each.
  for (size_t i = 0; i < lines.size(); ++i) {
    if (lines[i].find("T_world_origin") != std::string::npos) {
      std::vector<std::array<double, 4>> rows;
      for (size_t j = i + 1; j < lines.size() && rows.size() < 4; ++j) {
        std::istringstream iss(lines[j]);
        double a, b, c, d;
        if (iss >> a >> b >> c >> d) rows.push_back({a, b, c, d});
      }
      if (rows.size() == 4) {
        T_gt = Eigen::Matrix4d::Identity();
        for (int r = 0; r < 4; ++r)
          for (int c = 0; c < 4; ++c)
            T_gt(r, c) = rows[r][c];
        return true;
      }
    }
  }
  return false;
}

void reportError(const Eigen::Matrix4d& T_est, const Eigen::Matrix4d& T_gt,
                 const std::string& label) {
  Eigen::Matrix3d R_err = T_gt.block<3, 3>(0, 0).transpose() * T_est.block<3, 3>(0, 0);
  double cos_th = std::fmax(-1.0, std::fmin(1.0, (R_err.trace() - 1.0) / 2.0));
  double rot_err_deg = std::acos(cos_th) * 180.0 / M_PI;
  double t_err = (T_gt.block<3, 1>(0, 3) - T_est.block<3, 1>(0, 3)).norm();
  std::cout << label << ": Δrot = " << rot_err_deg << " deg, Δtrans = " << t_err << " m\n";
}

}  // namespace

int main(int argc, char** argv) {
  Args a;
  if (!parseArgs(argc, argv, a)) return 1;

  // ----- Stage 0: load -----
  std::cout << "[load] map  " << a.map_path << "\n";
  teaser::PointCloud map_cloud;
  if (!loadCloud(a.map_path, map_cloud)) { std::cerr << "map load failed\n"; return 1; }
  std::cout << "[load] scan " << a.scan_path << "\n";
  teaser::PointCloud scan_cloud;
  if (!loadCloud(a.scan_path, scan_cloud)) { std::cerr << "scan load failed\n"; return 1; }
  std::cout << "       map=" << map_cloud.size() << " pts, scan=" << scan_cloud.size() << " pts\n";

  Eigen::Matrix4d T_gt = Eigen::Matrix4d::Identity();
  bool has_gt = false;
  if (!a.gt_path.empty()) {
    has_gt = loadGT(a.gt_path, T_gt);
    if (has_gt) {
      std::cout << "[gt] loaded from " << a.gt_path << " (frame=" << a.gt_frame << ")\n";
      // If scan is submap_levelled, the input cloud has had R_pitch·R_roll baked in,
      // so the pose-to-recover only contains R_yaw + translation. Strip roll/pitch from GT.
      if (a.gt_frame == "levelled") {
        Eigen::Matrix3d R = T_gt.block<3, 3>(0, 0);
        double yaw = std::atan2(R(1, 0), R(0, 0));
        Eigen::AngleAxisd rz(yaw, Eigen::Vector3d::UnitZ());
        T_gt.block<3, 3>(0, 0) = rz.toRotationMatrix();
      }
    } else {
      std::cout << "[gt] not parseable, skipping error reporting\n";
    }
  }

  // ----- Stage A: BBS coarse -----
  std::cout << "\n=== Stage A: 3D-BBS coarse ===\n";
  auto map_eigen = toEigenVec(map_cloud);
  auto scan_eigen = toEigenVec(scan_cloud);
  auto map_ds  = voxelDown(map_eigen, a.bbs_tar_leaf);
  auto scan_ds = voxelDown(scan_eigen, a.bbs_src_leaf);
  std::cout << "  BBS map_ds=" << map_ds.size() << "  scan_ds=" << scan_ds.size() << "\n";

  gpu::BBS3D bbs;
  bbs.set_tar_points(map_ds, a.min_level_res, a.max_level);
  bbs.set_trans_search_range(map_ds);
  bbs.set_angular_search_range(
      Eigen::Vector3f(-0.02f, -0.02f, 0.0f),
      Eigen::Vector3f(+0.02f, +0.02f, static_cast<float>(2.0 * M_PI)));
  bbs.set_score_threshold_percentage(a.score_threshold);
  bbs.enable_timeout();
  bbs.set_timeout_duration_in_msec(a.bbs_timeout_ms);
  bbs.set_src_points(scan_ds);

  auto t0 = std::chrono::steady_clock::now();
  bbs.localize();
  auto t1 = std::chrono::steady_clock::now();
  double bbs_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

  std::cout << "  elapsed: " << bbs_ms << " ms\n";
  std::cout << "  has_localized: " << bbs.has_localized()
            << "  timed_out: " << bbs.has_timed_out() << "\n";
  std::cout << "  best score: " << bbs.get_best_score()
            << " (" << bbs.get_best_score_percentage() * 100.0f << "%)\n";

  if (!bbs.has_localized()) {
    std::cerr << "BBS did not converge.\n";
    return 2;
  }

  Eigen::Matrix4d T_coarse = bbs.get_global_pose().cast<double>();
  std::cout << "  T_coarse:\n" << T_coarse << "\n";
  if (has_gt) reportError(T_coarse, T_gt, "  [coarse vs GT]");

  // ----- Stage B: TEASER refine -----
  std::cout << "\n=== Stage B: TEASER refine ===\n";
  refine_pipeline::Config rc;
  rc.algorithm = a.algo;
  rc.noise_bound = a.noise_bound;
  rc.local_radius = a.local_radius;
  rc.voxel_size_map = a.voxel_size_map;
  rc.voxel_size_scan = a.voxel_size_scan;
  rc.fpfh_normal_radius = a.fpfh_normal_radius;
  rc.fpfh_radius = a.fpfh_radius;
  rc.max_correspondences = a.max_corres;

  Eigen::Matrix4d T_refined = T_coarse;
  refine_pipeline::Quality q;
  bool refined_ok = refine_pipeline::refineLocalPose(rc, map_cloud, scan_cloud, T_coarse, T_refined, q);
  if (!refined_ok) {
    std::cerr << "Refinement failed. Falling back to coarse pose.\n";
    T_refined = T_coarse;
  }

  std::cout << "\n=== Result ===\n";
  std::cout << "T_refined:\n" << T_refined << "\n";
  std::cout << "  correspondences: " << q.n_correspondences << "\n";
  std::cout << "  inliers:         " << q.n_inliers
            << " (" << q.inlier_ratio * 100 << "%)\n";
  std::cout << "  rmse:            " << q.rmse << " m\n";
  std::cout << "  BBS+TEASER total: " << (bbs_ms / 1000.0 + q.elapsed_s) << " s\n";
  if (has_gt) {
    reportError(T_coarse,  T_gt, "[coarse  vs GT]");
    reportError(T_refined, T_gt, "[refined vs GT]");
  }

  std::ofstream out(a.output_path);
  if (out.is_open()) {
    out << "# T_refined (4x4)\n" << T_refined << "\n";
    out << "# quality: bbs_ms refine_s n_corres n_inliers inlier_ratio rmse\n";
    out << bbs_ms << " " << q.elapsed_s << " " << q.n_correspondences << " "
        << q.n_inliers << " " << q.inlier_ratio << " " << q.rmse << "\n";
    if (has_gt) {
      Eigen::Matrix3d Re = T_gt.block<3,3>(0,0).transpose() * T_refined.block<3,3>(0,0);
      double cos_th = std::fmax(-1.0, std::fmin(1.0, (Re.trace() - 1.0) / 2.0));
      out << "# error: rot_deg trans_m\n";
      out << std::acos(cos_th) * 180.0 / M_PI << " "
          << (T_gt.block<3,1>(0,3) - T_refined.block<3,1>(0,3)).norm() << "\n";
    }
    std::cout << "Wrote " << a.output_path << "\n";
  }
  return refined_ok ? 0 : 3;
}
