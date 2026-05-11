// Phase 5 smoke test: confirm gpu_bbs3d links + runs a localize() end-to-end.
//
// Test design: take bun_zipper_res3.ply as the "map" (tar), apply a known
// (yaw + translation) transform to a copy of it as the "scan" (src),
// then ask BBS to recover the transform.
//
// Pass criteria:
//   - has_localized() returns true
//   - recovered pose is within (~min_level_res) of ground truth
//
// This isolates Phase 5's first risk: does our local BBS install actually run.

#include <chrono>
#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <gpu_bbs3d/bbs3d.cuh>
#include <teaser/ply_io.h>

namespace {

std::vector<Eigen::Vector3f> load_ply_points(const std::string& path) {
  teaser::PLYReader reader;
  teaser::PointCloud cloud;
  if (reader.read(path, cloud) != 0) {
    throw std::runtime_error("PLY read failed: " + path);
  }
  std::vector<Eigen::Vector3f> out;
  out.reserve(cloud.size());
  for (size_t i = 0; i < cloud.size(); ++i) {
    const auto& p = cloud[i];
    out.emplace_back(p.x, p.y, p.z);
  }
  return out;
}

std::vector<Eigen::Vector3f> apply_transform(
    const std::vector<Eigen::Vector3f>& pts,
    const Eigen::Matrix4f& T) {
  std::vector<Eigen::Vector3f> out;
  out.reserve(pts.size());
  const Eigen::Matrix3f R = T.block<3, 3>(0, 0);
  const Eigen::Vector3f t = T.block<3, 1>(0, 3);
  for (const auto& p : pts) out.emplace_back(R * p + t);
  return out;
}

}  // namespace

int main(int argc, char** argv) {
  const std::string ply_path = (argc > 1)
      ? std::string(argv[1])
      : "./example_data/bun_zipper_res3.ply";

  std::cout << "[bbs_smoke] loading " << ply_path << "\n";
  std::vector<Eigen::Vector3f> tar = load_ply_points(ply_path);
  std::cout << "[bbs_smoke] tar points: " << tar.size() << "\n";

  // Bunny is ~0.15 m wide. Scale by 100 so search resolution makes sense.
  for (auto& p : tar) p *= 100.0f;

  // Ground-truth transform: yaw=30°, t=(2.0, 1.5, 0.3)
  const float gt_yaw = static_cast<float>(M_PI) / 6.0f;
  Eigen::Matrix4f T_gt = Eigen::Matrix4f::Identity();
  T_gt.block<3, 3>(0, 0) =
      Eigen::AngleAxisf(gt_yaw, Eigen::Vector3f::UnitZ()).toRotationMatrix();
  T_gt.block<3, 1>(0, 3) = Eigen::Vector3f(2.0f, 1.5f, 0.3f);

  // src = inv(T_gt) · tar — i.e. tar = T_gt · src (BBS estimates pose of src in tar frame).
  std::vector<Eigen::Vector3f> src = apply_transform(tar, T_gt.inverse());

  gpu::BBS3D bbs;
  const float min_level_res = 0.5f;
  const int max_level = 6;
  bbs.set_tar_points(tar, min_level_res, max_level);
  bbs.set_trans_search_range(tar);
  // 3DOF: gravity-aligned, only yaw varies.
  bbs.set_angular_search_range(
      Eigen::Vector3f(-0.02f, -0.02f, 0.0f),
      Eigen::Vector3f(+0.02f, +0.02f, static_cast<float>(2.0 * M_PI)));
  bbs.set_score_threshold_percentage(0.9f);
  bbs.enable_timeout();
  bbs.set_timeout_duration_in_msec(10000);
  bbs.set_src_points(src);

  std::cout << "[bbs_smoke] running localize()...\n";
  const auto t0 = std::chrono::steady_clock::now();
  bbs.localize();
  const auto t1 = std::chrono::steady_clock::now();
  const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

  std::cout << "[bbs_smoke] elapsed: " << ms << " ms\n";
  std::cout << "[bbs_smoke] has_localized: " << bbs.has_localized() << "\n";
  std::cout << "[bbs_smoke] has_timed_out: " << bbs.has_timed_out() << "\n";
  std::cout << "[bbs_smoke] best_score: " << bbs.get_best_score()
            << " (" << bbs.get_best_score_percentage() * 100.0f << "% of src)\n";

  if (!bbs.has_localized()) {
    std::cerr << "[bbs_smoke] FAIL: localize() did not converge\n";
    return 1;
  }

  const Eigen::Matrix4f T_est = bbs.get_global_pose();
  std::cout << "[bbs_smoke] T_gt:\n" << T_gt << "\n";
  std::cout << "[bbs_smoke] T_est:\n" << T_est << "\n";

  // Error: T_err = T_gt^-1 · T_est should be ≈ identity.
  const Eigen::Matrix4f T_err = T_gt.inverse() * T_est;
  const float trans_err = T_err.block<3, 1>(0, 3).norm();
  const float angle_err_rad =
      Eigen::AngleAxisf(T_err.block<3, 3>(0, 0)).angle();
  const float angle_err_deg = angle_err_rad * 180.0f / static_cast<float>(M_PI);

  std::cout << "[bbs_smoke] translation error: " << trans_err << " m\n";
  std::cout << "[bbs_smoke] rotation error: " << angle_err_deg << " deg\n";

  // Tolerance = a few voxels of slack, since BBS output is coarse by design.
  const bool pass = trans_err < 3.0f * min_level_res && angle_err_deg < 5.0f;
  std::cout << "[bbs_smoke] " << (pass ? "PASS" : "FAIL") << "\n";
  return pass ? 0 : 2;
}
