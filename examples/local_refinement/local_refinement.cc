// Local refinement example: given a rough pose, refine to high precision using FPFH + TEASER.
//
// Use case (per airy_localization_plan.md demand 2):
//   - Pre-built map (PLY)
//   - Single (or accumulated) scan
//   - Rough pose with 3-10m / few-degrees error
//   - Goal: refine to <=0.1m / <=0.5deg
//
// Pipeline:
//   1. Extract local map around the rough pose
//   2. Pre-transform scan by rough pose into the map frame
//   3. Downsample both clouds
//   4. Compute FPFH features
//   5. Match features with cross-check + tuple test
//   6. Truncate correspondences if too many (avoid TEASER's O(N^2) blow-up)
//   7. Run TEASER (Quatro or full GNC_TLS) to get correction transform
//   8. Compose correction with rough pose -> refined pose
//
// Demo mode (no args): the bunny PLY is used as both map and scan, with a synthetic ground-truth
// transformation and rough-pose noise applied. The pipeline must recover the ground truth.
//
// Real-data mode: pass --map, --scan, --rough-pose paths. See README.md.

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <teaser/ply_io.h>
#include <teaser/registration.h>
#include <teaser/matcher.h>

namespace {

// ----- Configuration -----------------------------------------------------------------------------

struct Config {
  std::string map_path;
  std::string scan_path;
  std::string rough_pose_path;
  std::string output_path = "refined_pose.txt";
  std::string algorithm = "Quatro";   // Quatro | TEASER

  double local_radius = 60.0;          // local map crop radius (meters)
  double voxel_size_map = 0.3;         // map downsample voxel (meters)
  double voxel_size_scan = 0.2;        // scan downsample voxel (meters)
  double fpfh_normal_radius = 1.0;     // FPFH normal estimation radius
  double fpfh_radius = 2.0;            // FPFH descriptor radius
  double noise_bound = 0.1;            // TEASER noise bound (meters)
  int max_correspondences = 3000;      // hard cap to keep TEASER's O(N^2) memory bounded
  bool demo_mode = true;               // set false when --map/--scan/--rough-pose given
};

// ----- Pose IO -----------------------------------------------------------------------------------

inline Eigen::Matrix3d eulerRPYtoRot(double roll, double pitch, double yaw) {
  // Convention: R = Rz(yaw) * Ry(pitch) * Rx(roll), applied to a column vector.
  Eigen::AngleAxisd rx(roll, Eigen::Vector3d::UnitX());
  Eigen::AngleAxisd ry(pitch, Eigen::Vector3d::UnitY());
  Eigen::AngleAxisd rz(yaw, Eigen::Vector3d::UnitZ());
  return (rz * ry * rx).toRotationMatrix();
}

// Load 6DOF pose from "x y z roll pitch yaw" (radians).
bool loadPose(const std::string& path, Eigen::Matrix4d& T) {
  std::ifstream f(path);
  if (!f.is_open()) {
    std::cerr << "Failed to open pose file: " << path << "\n";
    return false;
  }
  double x, y, z, roll, pitch, yaw;
  if (!(f >> x >> y >> z >> roll >> pitch >> yaw)) {
    std::cerr << "Pose file must contain: x y z roll pitch yaw (radians)\n";
    return false;
  }
  T = Eigen::Matrix4d::Identity();
  T.block<3, 3>(0, 0) = eulerRPYtoRot(roll, pitch, yaw);
  T(0, 3) = x;
  T(1, 3) = y;
  T(2, 3) = z;
  return true;
}

void writePose(const std::string& path, const Eigen::Matrix4d& T, double inlier_ratio,
               double rmse, double elapsed_s) {
  std::ofstream f(path);
  Eigen::Matrix3d R = T.block<3, 3>(0, 0);
  Eigen::Vector3d t = T.block<3, 1>(0, 3);
  // Euler RPY back out of R (ZYX intrinsic: same as the loader convention above).
  // pitch in [-pi/2, pi/2]; standard atan2 extraction:
  double pitch = std::asin(-R(2, 0));
  double roll = std::atan2(R(2, 1), R(2, 2));
  double yaw = std::atan2(R(1, 0), R(0, 0));
  f << "# refined pose: x y z roll pitch yaw (radians)\n";
  f << t(0) << " " << t(1) << " " << t(2) << " "
    << roll << " " << pitch << " " << yaw << "\n";
  f << "# 4x4 matrix:\n";
  f << T << "\n";
  f << "# quality: inlier_ratio rmse_meters elapsed_seconds\n";
  f << inlier_ratio << " " << rmse << " " << elapsed_s << "\n";
}

// ----- Point cloud utilities ---------------------------------------------------------------------

teaser::PointCloud applyTransform(const teaser::PointCloud& cloud, const Eigen::Matrix4d& T) {
  teaser::PointCloud out;
  out.reserve(cloud.size());
  const Eigen::Matrix3d R = T.block<3, 3>(0, 0).cast<double>();
  const Eigen::Vector3d t = T.block<3, 1>(0, 3).cast<double>();
  for (const auto& p : cloud) {
    Eigen::Vector3d v(p.x, p.y, p.z);
    Eigen::Vector3d w = R * v + t;
    out.push_back({static_cast<float>(w(0)), static_cast<float>(w(1)), static_cast<float>(w(2))});
  }
  return out;
}

// Crop points within a sphere of `radius` centered at `center`.
teaser::PointCloud extractLocalMap(const teaser::PointCloud& map, const Eigen::Vector3d& center,
                                   double radius) {
  teaser::PointCloud out;
  const double r2 = radius * radius;
  for (const auto& p : map) {
    double dx = p.x - center(0);
    double dy = p.y - center(1);
    double dz = p.z - center(2);
    if (dx * dx + dy * dy + dz * dz <= r2) {
      out.push_back(p);
    }
  }
  return out;
}

// Simple voxel-grid downsampling: keep the first point seen per integer voxel cell.
// (This is intentionally light - PCL's true VoxelGrid averages within cells, but for FPFH
// feature stability the simpler "first-hit" approach is adequate and avoids the PCL dependency
// at the example level.)
teaser::PointCloud voxelDownsample(const teaser::PointCloud& cloud, double voxel) {
  if (voxel <= 0.0 || cloud.size() == 0) return cloud;
  struct Key {
    int x, y, z;
    bool operator==(const Key& o) const { return x == o.x && y == o.y && z == o.z; }
  };
  struct KeyHash {
    size_t operator()(const Key& k) const {
      auto h = std::hash<long long>();
      return h(static_cast<long long>(k.x) * 73856093LL ^
               static_cast<long long>(k.y) * 19349663LL ^
               static_cast<long long>(k.z) * 83492791LL);
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

// ----- TEASER parameter helper -------------------------------------------------------------------

teaser::RobustRegistrationSolver::Params makeParams(const Config& cfg) {
  teaser::RobustRegistrationSolver::Params params;
  params.noise_bound = cfg.noise_bound;
  params.cbar2 = 1.0;
  params.estimate_scaling = false;
  params.rotation_max_iterations = 100;
  params.rotation_gnc_factor = 1.4;
  params.rotation_cost_threshold = 2e-4;
  params.max_clique_time_limit = 10.0;  // 10 seconds, not the 3600 default
  if (cfg.algorithm == "Quatro") {
    params.rotation_estimation_algorithm =
        teaser::RobustRegistrationSolver::ROTATION_ESTIMATION_ALGORITHM::QUATRO;
    params.inlier_selection_mode =
        teaser::RobustRegistrationSolver::INLIER_SELECTION_MODE::PMC_HEU;
  } else {
    params.rotation_estimation_algorithm =
        teaser::RobustRegistrationSolver::ROTATION_ESTIMATION_ALGORITHM::GNC_TLS;
    params.inlier_selection_mode =
        teaser::RobustRegistrationSolver::INLIER_SELECTION_MODE::PMC_EXACT;
  }
  return params;
}

// ----- CLI parsing -------------------------------------------------------------------------------

void printUsage(const char* argv0) {
  std::cout
      << "Usage:\n"
      << "  " << argv0 << "                                  # demo mode (bunny self-test)\n"
      << "  " << argv0 << " --map MAP.ply --scan SCAN.ply --rough-pose POSE.txt [options]\n"
      << "\n"
      << "Options:\n"
      << "  --output PATH            output refined pose path (default refined_pose.txt)\n"
      << "  --algo Quatro|TEASER     rotation algorithm (default Quatro)\n"
      << "  --local-radius R         local map crop radius in meters (default 60)\n"
      << "  --voxel-map V            map voxel size (default 0.3)\n"
      << "  --voxel-scan V           scan voxel size (default 0.2)\n"
      << "  --fpfh-normal-r R        FPFH normal radius (default 1.0)\n"
      << "  --fpfh-r R               FPFH descriptor radius (default 2.0)\n"
      << "  --noise-bound NB         TEASER noise bound in meters (default 0.1)\n"
      << "  --max-corres N           cap on correspondences to keep N^2 bounded (default 3000)\n";
}

bool parseArgs(int argc, char** argv, Config& cfg) {
  if (argc == 1) return true;  // demo
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    auto need = [&](const char* opt) {
      if (i + 1 >= argc) {
        std::cerr << opt << " requires a value\n";
        std::exit(1);
      }
      return std::string(argv[++i]);
    };
    if (a == "--map") { cfg.map_path = need("--map"); cfg.demo_mode = false; }
    else if (a == "--scan") { cfg.scan_path = need("--scan"); cfg.demo_mode = false; }
    else if (a == "--rough-pose") { cfg.rough_pose_path = need("--rough-pose"); cfg.demo_mode = false; }
    else if (a == "--output") { cfg.output_path = need("--output"); }
    else if (a == "--algo") { cfg.algorithm = need("--algo"); }
    else if (a == "--local-radius") { cfg.local_radius = std::stod(need("--local-radius")); }
    else if (a == "--voxel-map") { cfg.voxel_size_map = std::stod(need("--voxel-map")); }
    else if (a == "--voxel-scan") { cfg.voxel_size_scan = std::stod(need("--voxel-scan")); }
    else if (a == "--fpfh-normal-r") { cfg.fpfh_normal_radius = std::stod(need("--fpfh-normal-r")); }
    else if (a == "--fpfh-r") { cfg.fpfh_radius = std::stod(need("--fpfh-r")); }
    else if (a == "--noise-bound") { cfg.noise_bound = std::stod(need("--noise-bound")); }
    else if (a == "--max-corres") { cfg.max_correspondences = std::stoi(need("--max-corres")); }
    else if (a == "-h" || a == "--help") { printUsage(argv[0]); std::exit(0); }
    else { std::cerr << "Unknown argument: " << a << "\n"; printUsage(argv[0]); std::exit(1); }
  }
  if (!cfg.demo_mode) {
    if (cfg.map_path.empty() || cfg.scan_path.empty() || cfg.rough_pose_path.empty()) {
      std::cerr << "Real-data mode requires --map, --scan, --rough-pose\n";
      return false;
    }
  }
  return true;
}

// ----- Demo data generation ----------------------------------------------------------------------

// Loads the bunny, uses the whole thing as "map", and crops one half as "scan".
// Applies a known ground-truth SE(3), then adds noise to it to form the rough pose.
struct DemoSetup {
  teaser::PointCloud map_cloud;
  teaser::PointCloud scan_cloud;
  Eigen::Matrix4d T_gt;
  Eigen::Matrix4d T_rough;
};

bool buildDemo(DemoSetup& demo) {
  teaser::PLYReader reader;
  // Bunny is loaded relative to executable cwd; build sets up example_data/ next to it.
  std::string bunny_path = "./example_data/bun_zipper_res3.ply";
  if (reader.read(bunny_path, demo.map_cloud) != 0) {
    std::cerr << "Failed to read bunny PLY at " << bunny_path << "\n";
    return false;
  }

  // Use the entire bunny as the scan (in scan frame). The map is the bunny placed in some
  // arbitrary pose; we'll transform the bunny BACKWARDS by an inverse to fabricate a map.
  // Equivalently: scan = bunny, map = T_gt * bunny.
  demo.scan_cloud = demo.map_cloud;  // copy

  // Ground-truth SE(3) - small but non-trivial rotation + translation.
  std::random_device rd;
  std::mt19937 gen(42);
  std::uniform_real_distribution<double> ang(-15.0 * M_PI / 180.0, 15.0 * M_PI / 180.0);
  std::uniform_real_distribution<double> trn(-0.05, 0.05);

  double roll_gt = ang(gen), pitch_gt = ang(gen), yaw_gt = ang(gen);
  Eigen::Matrix4d T_gt = Eigen::Matrix4d::Identity();
  T_gt.block<3, 3>(0, 0) = eulerRPYtoRot(roll_gt, pitch_gt, yaw_gt);
  T_gt(0, 3) = trn(gen);
  T_gt(1, 3) = trn(gen);
  T_gt(2, 3) = trn(gen);
  demo.T_gt = T_gt;

  // Build map by applying T_gt to bunny: map_point = T_gt * scan_point.
  demo.map_cloud = applyTransform(demo.scan_cloud, T_gt);

  // Rough pose: ground truth perturbed by a few-cm translation (bunny is ~10cm, so
  // 3-10m translation noise would make scan disjoint from local map. Scale down to bunny size).
  std::uniform_real_distribution<double> rough_ang(-3.0 * M_PI / 180.0, 3.0 * M_PI / 180.0);
  std::uniform_real_distribution<double> rough_trn(-0.02, 0.02);
  double droll = rough_ang(gen), dpitch = rough_ang(gen), dyaw = rough_ang(gen);
  Eigen::Matrix4d T_rough = T_gt;
  T_rough.block<3, 3>(0, 0) = T_gt.block<3, 3>(0, 0) *
                              eulerRPYtoRot(droll, dpitch, dyaw);
  T_rough(0, 3) += rough_trn(gen);
  T_rough(1, 3) += rough_trn(gen);
  T_rough(2, 3) += rough_trn(gen);
  demo.T_rough = T_rough;
  return true;
}

// ----- Main pipeline -----------------------------------------------------------------------------

struct Quality {
  int n_correspondences = 0;
  int n_inliers = 0;
  double inlier_ratio = 0.0;
  double rmse = 0.0;
  double elapsed_s = 0.0;
};

// Run the full refinement on (map, scan, T_rough). Returns refined T and quality.
bool refineLocalPose(const Config& cfg, const teaser::PointCloud& map,
                     const teaser::PointCloud& scan, const Eigen::Matrix4d& T_rough,
                     Eigen::Matrix4d& T_refined, Quality& q) {
  const auto t0 = std::chrono::steady_clock::now();

  // 1. Crop local map around rough pose origin.
  Eigen::Vector3d center = T_rough.block<3, 1>(0, 3);
  teaser::PointCloud local_map = extractLocalMap(map, center, cfg.local_radius);
  std::cout << "  local map: " << local_map.size() << " / " << map.size() << " points\n";
  if (local_map.size() < 50) {
    std::cerr << "Local map too small. Increase --local-radius.\n";
    return false;
  }

  // 2. Pre-transform scan into the map frame using T_rough.
  teaser::PointCloud scan_initial = applyTransform(scan, T_rough);

  // 3. Downsample.
  teaser::PointCloud local_map_ds = voxelDownsample(local_map, cfg.voxel_size_map);
  teaser::PointCloud scan_ds = voxelDownsample(scan_initial, cfg.voxel_size_scan);
  std::cout << "  after downsample: scan " << scan_ds.size()
            << ", local_map " << local_map_ds.size() << "\n";

  // 4. FPFH features.
  teaser::FPFHEstimation fpfh;
  auto scan_features = fpfh.computeFPFHFeatures(scan_ds, cfg.fpfh_normal_radius, cfg.fpfh_radius);
  auto map_features = fpfh.computeFPFHFeatures(local_map_ds, cfg.fpfh_normal_radius, cfg.fpfh_radius);

  // 5. Match.
  teaser::Matcher matcher;
  auto correspondences = matcher.calculateCorrespondences(
      scan_ds, local_map_ds, *scan_features, *map_features,
      /*use_absolute_scale=*/false, /*use_crosscheck=*/true,
      /*use_tuple_test=*/true, /*tuple_scale=*/0.95);
  q.n_correspondences = static_cast<int>(correspondences.size());
  std::cout << "  correspondences: " << correspondences.size() << "\n";

  if (correspondences.size() < 6) {
    std::cerr << "Too few correspondences (need >=6).\n";
    return false;
  }

  // 6. Cap correspondences to bound TEASER's O(N^2) memory.
  if (static_cast<int>(correspondences.size()) > cfg.max_correspondences) {
    std::cout << "  truncating to " << cfg.max_correspondences << " correspondences\n";
    correspondences.resize(cfg.max_correspondences);
  }

  // 7. TEASER. Construct a fresh solver each call (Quatro has a static-variable hazard if reused).
  auto params = makeParams(cfg);
  teaser::RobustRegistrationSolver solver(params);
  solver.solve(scan_ds, local_map_ds, correspondences);
  auto sol = solver.getSolution();
  if (!sol.valid) {
    std::cerr << "TEASER returned invalid solution (max-clique too small).\n";
    return false;
  }

  // 8. Compose: scan_in_map = T_correction * scan_initial = T_correction * T_rough * scan.
  Eigen::Matrix4d T_correction = Eigen::Matrix4d::Identity();
  T_correction.block<3, 3>(0, 0) = sol.rotation;
  T_correction.block<3, 1>(0, 3) = sol.translation;
  T_refined = T_correction * T_rough;

  // 9. Compute quality metrics on the final correspondences.
  q.n_inliers = static_cast<int>(solver.getTranslationInliers().size());
  q.inlier_ratio = correspondences.empty()
                       ? 0.0
                       : static_cast<double>(q.n_inliers) / correspondences.size();

  // Build inlier residual RMSE.
  double sq_sum = 0.0;
  const auto& inlier_idx = solver.getInputOrderedTranslationInliers();
  const auto& corr = correspondences;
  for (int idx : inlier_idx) {
    Eigen::Vector3d a(scan_ds[corr[idx].first].x, scan_ds[corr[idx].first].y,
                      scan_ds[corr[idx].first].z);
    Eigen::Vector3d b(local_map_ds[corr[idx].second].x, local_map_ds[corr[idx].second].y,
                      local_map_ds[corr[idx].second].z);
    Eigen::Vector3d r = (sol.rotation * a + sol.translation) - b;
    sq_sum += r.squaredNorm();
  }
  q.rmse = q.n_inliers > 0 ? std::sqrt(sq_sum / q.n_inliers) : 0.0;

  const auto t1 = std::chrono::steady_clock::now();
  q.elapsed_s = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1e6;
  return true;
}

}  // namespace

int main(int argc, char** argv) {
  Config cfg;
  if (!parseArgs(argc, argv, cfg)) return 1;

  teaser::PointCloud map_cloud, scan_cloud;
  Eigen::Matrix4d T_rough = Eigen::Matrix4d::Identity();
  Eigen::Matrix4d T_gt = Eigen::Matrix4d::Identity();
  bool has_gt = false;

  if (cfg.demo_mode) {
    std::cout << "[demo mode] using bunny PLY as both map and scan with synthetic noise.\n";
    // Demo scale is tiny (bunny ~0.1m). Scale parameters down for the demo to make sense.
    cfg.local_radius = 1.0;
    cfg.voxel_size_map = 0.005;
    cfg.voxel_size_scan = 0.005;
    cfg.fpfh_normal_radius = 0.02;
    cfg.fpfh_radius = 0.04;
    cfg.noise_bound = 0.01;

    DemoSetup demo;
    if (!buildDemo(demo)) return 1;
    map_cloud = std::move(demo.map_cloud);
    scan_cloud = std::move(demo.scan_cloud);
    T_rough = demo.T_rough;
    T_gt = demo.T_gt;
    has_gt = true;
  } else {
    teaser::PLYReader reader;
    if (reader.read(cfg.map_path, map_cloud) != 0) {
      std::cerr << "Failed to read map: " << cfg.map_path << "\n";
      return 1;
    }
    if (reader.read(cfg.scan_path, scan_cloud) != 0) {
      std::cerr << "Failed to read scan: " << cfg.scan_path << "\n";
      return 1;
    }
    if (!loadPose(cfg.rough_pose_path, T_rough)) return 1;
    std::cout << "Loaded map " << map_cloud.size() << " pts, scan " << scan_cloud.size()
              << " pts.\n";
  }

  std::cout << "Algorithm: " << cfg.algorithm << "\n";
  std::cout << "Rough pose:\n" << T_rough << "\n";

  Eigen::Matrix4d T_refined;
  Quality q;
  if (!refineLocalPose(cfg, map_cloud, scan_cloud, T_rough, T_refined, q)) {
    return 1;
  }

  std::cout << "\n=====================================\n";
  std::cout << "       Local Refinement Result       \n";
  std::cout << "=====================================\n";
  std::cout << "Refined pose:\n" << T_refined << "\n";
  std::cout << "Correspondences: " << q.n_correspondences << "\n";
  std::cout << "Inliers:         " << q.n_inliers << " (" << q.inlier_ratio * 100 << "%)\n";
  std::cout << "Inlier RMSE:     " << q.rmse << " m\n";
  std::cout << "Elapsed:         " << q.elapsed_s << " s\n";

  if (has_gt) {
    Eigen::Matrix3d R_err = T_gt.block<3, 3>(0, 0).transpose() * T_refined.block<3, 3>(0, 0);
    double cos_th = std::fmax(-1.0, std::fmin(1.0, (R_err.trace() - 1.0) / 2.0));
    double rot_err_deg = std::acos(cos_th) * 180.0 / M_PI;
    double t_err = (T_gt.block<3, 1>(0, 3) - T_refined.block<3, 1>(0, 3)).norm();
    std::cout << "\nGround truth error: rotation = " << rot_err_deg << " deg, translation = "
              << t_err << " m\n";
  }

  writePose(cfg.output_path, T_refined, q.inlier_ratio, q.rmse, q.elapsed_s);
  std::cout << "\nWrote " << cfg.output_path << "\n";
  return 0;
}
