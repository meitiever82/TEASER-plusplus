// FPFH + TEASER refinement helpers shared by Phase 1 (local_refinement) and Phase 5 (BBS).
// Logic is a verbatim lift from examples/local_refinement/local_refinement.cc, exposed as
// inline functions so multiple translation units can include this safely.

#pragma once

#include <chrono>
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <teaser/registration.h>
#include <teaser/matcher.h>

namespace refine_pipeline {

struct Config {
  std::string algorithm = "TEASER";    // Quatro | TEASER
  double local_radius = 60.0;
  double voxel_size_map = 0.3;
  double voxel_size_scan = 0.2;
  double fpfh_normal_radius = 1.0;
  double fpfh_radius = 2.0;
  double noise_bound = 0.15;
  int max_correspondences = 3000;
};

struct Quality {
  int n_correspondences = 0;
  int n_inliers = 0;
  double inlier_ratio = 0.0;
  double rmse = 0.0;
  double elapsed_s = 0.0;
};

inline teaser::PointCloud applyTransform(const teaser::PointCloud& cloud,
                                         const Eigen::Matrix4d& T) {
  teaser::PointCloud out;
  out.reserve(cloud.size());
  const Eigen::Matrix3d R = T.block<3, 3>(0, 0);
  const Eigen::Vector3d t = T.block<3, 1>(0, 3);
  for (const auto& p : cloud) {
    Eigen::Vector3d v(p.x, p.y, p.z);
    Eigen::Vector3d w = R * v + t;
    out.push_back({static_cast<float>(w(0)),
                   static_cast<float>(w(1)),
                   static_cast<float>(w(2))});
  }
  return out;
}

inline teaser::PointCloud extractLocalMap(const teaser::PointCloud& map,
                                          const Eigen::Vector3d& center,
                                          double radius) {
  teaser::PointCloud out;
  const double r2 = radius * radius;
  for (const auto& p : map) {
    double dx = p.x - center(0);
    double dy = p.y - center(1);
    double dz = p.z - center(2);
    if (dx * dx + dy * dy + dz * dz <= r2) out.push_back(p);
  }
  return out;
}

inline teaser::PointCloud voxelDownsample(const teaser::PointCloud& cloud, double voxel) {
  if (voxel <= 0.0 || cloud.size() == 0) return cloud;
  struct Key { int x, y, z; bool operator==(const Key& o) const { return x == o.x && y == o.y && z == o.z; } };
  struct KeyHash {
    size_t operator()(const Key& k) const {
      return std::hash<long long>()(
          static_cast<long long>(k.x) * 73856093LL
          ^ static_cast<long long>(k.y) * 19349663LL
          ^ static_cast<long long>(k.z) * 83492791LL);
    }
  };
  std::unordered_map<Key, teaser::PointXYZ, KeyHash> seen;
  seen.reserve(cloud.size());
  const double inv = 1.0 / voxel;
  for (const auto& p : cloud) {
    Key k{static_cast<int>(std::floor(p.x * inv)),
          static_cast<int>(std::floor(p.y * inv)),
          static_cast<int>(std::floor(p.z * inv))};
    seen.emplace(k, p);
  }
  teaser::PointCloud out;
  out.reserve(seen.size());
  for (auto& kv : seen) out.push_back(kv.second);
  return out;
}

inline teaser::RobustRegistrationSolver::Params makeParams(const Config& cfg) {
  teaser::RobustRegistrationSolver::Params p;
  p.noise_bound = cfg.noise_bound;
  p.cbar2 = 1.0;
  p.estimate_scaling = false;
  p.rotation_max_iterations = 100;
  p.rotation_gnc_factor = 1.4;
  p.rotation_cost_threshold = 2e-4;
  p.max_clique_time_limit = 10.0;
  if (cfg.algorithm == "Quatro") {
    p.rotation_estimation_algorithm =
        teaser::RobustRegistrationSolver::ROTATION_ESTIMATION_ALGORITHM::QUATRO;
    p.inlier_selection_mode =
        teaser::RobustRegistrationSolver::INLIER_SELECTION_MODE::PMC_HEU;
  } else {
    p.rotation_estimation_algorithm =
        teaser::RobustRegistrationSolver::ROTATION_ESTIMATION_ALGORITHM::GNC_TLS;
    p.inlier_selection_mode =
        teaser::RobustRegistrationSolver::INLIER_SELECTION_MODE::PMC_EXACT;
  }
  return p;
}

// Refine (map, scan, T_rough) -> T_refined.  Returns false if pipeline fails.
inline bool refineLocalPose(const Config& cfg, const teaser::PointCloud& map,
                            const teaser::PointCloud& scan,
                            const Eigen::Matrix4d& T_rough,
                            Eigen::Matrix4d& T_refined, Quality& q) {
  const auto t0 = std::chrono::steady_clock::now();

  Eigen::Vector3d center = T_rough.block<3, 1>(0, 3);
  teaser::PointCloud local_map = extractLocalMap(map, center, cfg.local_radius);
  std::cout << "  local map: " << local_map.size() << " / " << map.size() << " pts\n";
  if (local_map.size() < 50) {
    std::cerr << "Local map too small. Increase --local-radius.\n";
    return false;
  }

  teaser::PointCloud scan_initial = applyTransform(scan, T_rough);
  teaser::PointCloud local_map_ds = voxelDownsample(local_map, cfg.voxel_size_map);
  teaser::PointCloud scan_ds = voxelDownsample(scan_initial, cfg.voxel_size_scan);
  std::cout << "  after downsample: scan " << scan_ds.size()
            << ", local_map " << local_map_ds.size() << "\n";

  teaser::FPFHEstimation fpfh;
  auto scan_feat = fpfh.computeFPFHFeatures(scan_ds, cfg.fpfh_normal_radius, cfg.fpfh_radius);
  auto map_feat  = fpfh.computeFPFHFeatures(local_map_ds, cfg.fpfh_normal_radius, cfg.fpfh_radius);

  teaser::Matcher matcher;
  auto correspondences = matcher.calculateCorrespondences(
      scan_ds, local_map_ds, *scan_feat, *map_feat,
      false, true, true, 0.95);
  q.n_correspondences = static_cast<int>(correspondences.size());
  std::cout << "  correspondences: " << correspondences.size() << "\n";
  if (correspondences.size() < 6) {
    std::cerr << "Too few correspondences (need >=6).\n";
    return false;
  }
  if (static_cast<int>(correspondences.size()) > cfg.max_correspondences) {
    correspondences.resize(cfg.max_correspondences);
    std::cout << "  truncated to " << cfg.max_correspondences << "\n";
  }

  auto params = makeParams(cfg);
  teaser::RobustRegistrationSolver solver(params);
  solver.solve(scan_ds, local_map_ds, correspondences);
  auto sol = solver.getSolution();
  if (!sol.valid) {
    std::cerr << "TEASER returned invalid solution (max-clique too small).\n";
    return false;
  }

  Eigen::Matrix4d T_correction = Eigen::Matrix4d::Identity();
  T_correction.block<3, 3>(0, 0) = sol.rotation;
  T_correction.block<3, 1>(0, 3) = sol.translation;
  T_refined = T_correction * T_rough;

  q.n_inliers = static_cast<int>(solver.getTranslationInliers().size());
  q.inlier_ratio = correspondences.empty()
                     ? 0.0
                     : static_cast<double>(q.n_inliers) / correspondences.size();

  double sq_sum = 0.0;
  const auto& inlier_idx = solver.getInputOrderedTranslationInliers();
  for (int idx : inlier_idx) {
    const auto& a = scan_ds[correspondences[idx].first];
    const auto& b = local_map_ds[correspondences[idx].second];
    Eigen::Vector3d va(a.x, a.y, a.z), vb(b.x, b.y, b.z);
    Eigen::Vector3d r = (sol.rotation * va + sol.translation) - vb;
    sq_sum += r.squaredNorm();
  }
  q.rmse = q.n_inliers > 0 ? std::sqrt(sq_sum / q.n_inliers) : 0.0;

  const auto t1 = std::chrono::steady_clock::now();
  q.elapsed_s = std::chrono::duration<double>(t1 - t0).count();
  return true;
}

}  // namespace refine_pipeline
