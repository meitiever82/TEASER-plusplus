// Global localization example: Scan Context retrieval + FPFH/TEASER refinement.
//
// Use case (per airy_localization_plan.md demand 1):
//   - Pre-built map (PLY)
//   - Single (or accumulated) scan
//   - No prior pose at all
//   - Goal: 6DOF pose in the map frame, ~0.1m precision
//
// Two-stage pipeline:
//   Stage A (Scan Context retrieval): grid-sample anchor positions on the map (offline mode),
//     compute one SC descriptor per anchor, store. At query time compute SC for the scan and
//     retrieve top-K anchors via ringkey kNN + pairwise SC distance.
//   Stage B (per-candidate verification): for each top-K candidate, treat the anchor position
//     (plus SC's yaw estimate) as a rough pose and run FPFH + TEASER (the same pipeline as
//     local_refinement). Score by inlier ratio + RMSE; pick the best.
//
// Subcommands (selected with --mode):
//   build      Sample anchors from --map, write database to --db.
//   query      Just run SC retrieval against --db (for debugging - no TEASER step).
//   localize   Full pipeline: --map + --db + --scan -> refined 6DOF pose.
//
// See README.md for end-to-end recipes.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>

#include <teaser/ply_io.h>
#include <teaser/registration.h>
#include <teaser/matcher.h>

#include "third_party/scancontext/Scancontext.h"

namespace {

// ----- Configuration -----------------------------------------------------------------------------

struct Config {
  std::string mode;                    // "build" | "query" | "localize"

  std::string map_path;
  std::string scan_path;
  std::string db_path = "sc_database.bin";
  std::string output_path = "refined_pose.txt";

  // Offline DB building
  double anchor_step = 10.0;           // grid sample step in meters
  double anchor_z_max_dev = 5.0;       // discard anchors whose nearest map z deviates > this
  double sc_radius = 80.0;             // matches SCManager::PC_MAX_RADIUS - keep in sync
  double sc_downsample = 0.5;          // voxel downsample for SC input (SC is shape-based)

  // Retrieval
  int top_k = 5;                       // top-K candidates from SC to verify

  // Per-candidate TEASER refinement (mirrors local_refinement defaults)
  std::string algorithm = "Quatro";
  double local_radius = 60.0;
  double voxel_size_map = 0.3;
  double voxel_size_scan = 0.2;
  double fpfh_normal_radius = 1.0;
  double fpfh_radius = 2.0;
  double noise_bound = 0.1;
  int max_correspondences = 3000;
  int verify_min_inliers = 10;         // hard floor on absolute inlier count
};

// ----- Pose IO -----------------------------------------------------------------------------------

inline Eigen::Matrix3d eulerRPYtoRot(double roll, double pitch, double yaw) {
  Eigen::AngleAxisd rx(roll, Eigen::Vector3d::UnitX());
  Eigen::AngleAxisd ry(pitch, Eigen::Vector3d::UnitY());
  Eigen::AngleAxisd rz(yaw, Eigen::Vector3d::UnitZ());
  return (rz * ry * rx).toRotationMatrix();
}

void writePose(const std::string& path, const Eigen::Matrix4d& T, int anchor_id,
               double inlier_ratio, double rmse, double elapsed_s) {
  std::ofstream f(path);
  Eigen::Matrix3d R = T.block<3, 3>(0, 0);
  Eigen::Vector3d t = T.block<3, 1>(0, 3);
  double pitch = std::asin(-R(2, 0));
  double roll = std::atan2(R(2, 1), R(2, 2));
  double yaw = std::atan2(R(1, 0), R(0, 0));
  f << "# winning anchor: " << anchor_id << "\n";
  f << "# refined pose: x y z roll pitch yaw (radians)\n";
  f << t(0) << " " << t(1) << " " << t(2) << " "
    << roll << " " << pitch << " " << yaw << "\n";
  f << "# 4x4 matrix:\n";
  f << T << "\n";
  f << "# quality: inlier_ratio rmse_meters elapsed_seconds\n";
  f << inlier_ratio << " " << rmse << " " << elapsed_s << "\n";
}

// ----- Point cloud utilities (shared with local_refinement) ------------------------------------

// Load .ply (via teaser::PLYReader) or .pcd (via PCL) based on file extension.
// Returns true on success. Drops any extra PCD fields - we only keep xyz.
bool loadCloud(const std::string& path, teaser::PointCloud& out) {
  auto ends_with = [](const std::string& s, const std::string& suf) {
    return s.size() >= suf.size() &&
           std::equal(suf.rbegin(), suf.rend(), s.rbegin(),
                      [](char a, char b) { return std::tolower(a) == std::tolower(b); });
  };
  if (ends_with(path, ".ply")) {
    teaser::PLYReader reader;
    return reader.read(path, out) == 0;
  }
  if (ends_with(path, ".pcd")) {
    pcl::PointCloud<pcl::PointXYZ> tmp;
    if (pcl::io::loadPCDFile<pcl::PointXYZ>(path, tmp) < 0) return false;
    out.clear();
    out.reserve(tmp.size());
    for (const auto& p : tmp.points) {
      if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z)) continue;
      out.push_back({p.x, p.y, p.z});
    }
    std::cout << "\tRead " << out.size() << " total vertices (from PCD)" << std::endl;
    return true;
  }
  std::cerr << "Unsupported file extension (need .ply or .pcd): " << path << "\n";
  return false;
}

teaser::PointCloud applyTransform(const teaser::PointCloud& cloud, const Eigen::Matrix4d& T) {
  teaser::PointCloud out;
  out.reserve(cloud.size());
  const Eigen::Matrix3d R = T.block<3, 3>(0, 0);
  const Eigen::Vector3d t = T.block<3, 1>(0, 3);
  for (const auto& p : cloud) {
    Eigen::Vector3d v(p.x, p.y, p.z);
    Eigen::Vector3d w = R * v + t;
    out.push_back({static_cast<float>(w(0)), static_cast<float>(w(1)), static_cast<float>(w(2))});
  }
  return out;
}

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

teaser::PointCloud voxelDownsample(const teaser::PointCloud& cloud, double voxel) {
  if (voxel <= 0.0 || cloud.size() == 0) return cloud;
  struct Key { int x, y, z;
    bool operator==(const Key& o) const { return x == o.x && y == o.y && z == o.z; } };
  struct KeyHash { size_t operator()(const Key& k) const {
    auto h = std::hash<long long>();
    return h(static_cast<long long>(k.x) * 73856093LL ^
             static_cast<long long>(k.y) * 19349663LL ^
             static_cast<long long>(k.z) * 83492791LL);
  } };
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

// Convert teaser::PointCloud (XYZ float) to pcl::PointCloud<PointXYZI> (XYZ + zero intensity)
// for SCManager input. SC only consumes XYZ; intensity is unused.
pcl::PointCloud<pcl::PointXYZI>::Ptr toPclXYZI(const teaser::PointCloud& src) {
  auto out = std::make_shared<pcl::PointCloud<pcl::PointXYZI>>();
  out->reserve(src.size());
  for (const auto& p : src) {
    pcl::PointXYZI pi;
    pi.x = p.x; pi.y = p.y; pi.z = p.z; pi.intensity = 0.0f;
    out->push_back(pi);
  }
  return out;
}

// ----- TEASER parameters (matches local_refinement) --------------------------------------------

teaser::RobustRegistrationSolver::Params makeTeaserParams(const Config& cfg) {
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

// ----- Scan Context database -------------------------------------------------------------------
//
// Built directly from the SCManager primitives (makeScancontext, makeRingkeyFromScancontext,
// distanceBtnScanContext). We don't call detectLoopClosureID because it bakes in loop-closure
// assumptions (NUM_EXCLUDE_RECENT, SC_DIST_THRES, single-best return) that don't fit global
// localization.

struct Anchor {
  int id = 0;
  Eigen::Vector3d position = Eigen::Vector3d::Zero();
  Eigen::MatrixXd sc_desc;             // 20 x 60 (matches SCManager defaults)
  std::vector<float> ringkey;          // length 20
};

struct Match {
  int anchor_id;
  Eigen::Vector3d anchor_position;
  double sc_distance;
  double yaw_rad;                      // SC's yaw alignment estimate
};

class SCDatabase {
public:
  SCDatabase() : sc_() {}

  // Build by grid-sampling anchor positions over the XY bounding box of the map.
  // At each anchor we extract a local cloud of radius sc_radius, compute SC, and store.
  // Anchors with too-few points are skipped.
  void build(const teaser::PointCloud& map, double sample_step, double sc_radius,
             double sc_downsample, double z_max_dev) {
    if (map.size() == 0) return;

    // Configure SC's internal polar-binning radius to match the geometric crop.
    // Without this, SC keeps PC_MAX_RADIUS=80m default and bins are mismatched.
    sc_.setRadius(sc_radius);
    sc_radius_ = sc_radius;
    sc_.setLidarHeight(0.0);  // anchors are already centered at the local-frame origin

    double xmin =  1e18, xmax = -1e18;
    double ymin =  1e18, ymax = -1e18;
    double zmin =  1e18, zmax = -1e18;
    double zsum = 0.0;
    for (const auto& p : map) {
      xmin = std::min(xmin, double(p.x)); xmax = std::max(xmax, double(p.x));
      ymin = std::min(ymin, double(p.y)); ymax = std::max(ymax, double(p.y));
      zmin = std::min(zmin, double(p.z)); zmax = std::max(zmax, double(p.z));
      zsum += p.z;
    }
    double z_anchor = zsum / map.size();  // anchor sits at average map height

    std::cout << "Map bounds: x[" << xmin << ", " << xmax << "] y[" << ymin << ", " << ymax
              << "] z[" << zmin << ", " << zmax << "]\n";

    int n_skipped = 0;
    int next_id = 0;
    for (double y = ymin; y <= ymax; y += sample_step) {
      for (double x = xmin; x <= xmax; x += sample_step) {
        Eigen::Vector3d c(x, y, z_anchor);

        // Crop local cloud at this anchor.
        auto local = extractLocalMap(map, c, sc_radius);
        if (local.size() < 200) {  // SC needs reasonable coverage
          ++n_skipped;
          continue;
        }

        // Downsample for SC (it's a shape descriptor; raw density adds noise to bin heights).
        auto local_ds = voxelDownsample(local, sc_downsample);

        // SC's makeScancontext assumes the cloud is in the *robot's* local frame.
        // For each anchor we translate the local crop so the anchor sits at origin.
        for (auto& p : local_ds) {
          p.x -= static_cast<float>(c(0));
          p.y -= static_cast<float>(c(1));
          // z is left as-is; SC uses LIDAR_HEIGHT to add an offset, see Scancontext.h.
        }

        // SC input expects pcl::PointXYZI.
        auto pcl_cloud = toPclXYZI(local_ds);

        Anchor a;
        a.id = next_id++;
        a.position = c;
        a.sc_desc = sc_.makeScancontext(*pcl_cloud);
        Eigen::MatrixXd ring = sc_.makeRingkeyFromScancontext(a.sc_desc);
        a.ringkey.resize(ring.rows());
        for (int i = 0; i < ring.rows(); ++i) a.ringkey[i] = static_cast<float>(ring(i, 0));

        anchors_.push_back(std::move(a));
        if (anchors_.size() % 50 == 0) {
          std::cout << "  built " << anchors_.size() << " anchors...\n";
        }
      }
    }
    std::cout << "DB built: " << anchors_.size() << " anchors (skipped " << n_skipped
              << " with sparse coverage).\n";
    (void)z_max_dev;  // reserved for future per-anchor z filtering
  }

  // Top-K retrieval: brute-force pairwise SC distance.
  // For a few thousand anchors this is fast (each anchor's SC is 20x60 doubles, dist is O(rows*cols)).
  // If needed, a nanoflann tree on ringkeys could prune candidates first - but for our scale,
  // the brute-force is fine and simpler.
  std::vector<Match> queryTopK(const teaser::PointCloud& scan, double sc_downsample, int k) {
    // Keep SC's polar-binning radius synced with what was used at build time.
    if (sc_radius_ > 0) sc_.setRadius(sc_radius_);
    sc_.setLidarHeight(0.0);

    auto scan_ds = voxelDownsample(scan, sc_downsample);
    auto pcl_scan = toPclXYZI(scan_ds);
    Eigen::MatrixXd query_sc = sc_.makeScancontext(*pcl_scan);

    std::vector<Match> all;
    all.reserve(anchors_.size());
    for (auto& a : anchors_) {
      Eigen::MatrixXd cand = a.sc_desc;  // distanceBtnScanContext takes non-const ref
      auto dist_yaw = sc_.distanceBtnScanContext(query_sc, cand);
      Match m;
      m.anchor_id = a.id;
      m.anchor_position = a.position;
      m.sc_distance = dist_yaw.first;
      // distanceBtnScanContext returns yaw as a sector index; convert to radians.
      // Each sector is 360 / PC_NUM_SECTOR (= 6) degrees.
      double yaw_deg = dist_yaw.second * (360.0 / 60.0);
      m.yaw_rad = yaw_deg * M_PI / 180.0;
      all.push_back(m);
    }
    std::sort(all.begin(), all.end(),
              [](const Match& a, const Match& b) { return a.sc_distance < b.sc_distance; });
    if (k > 0 && static_cast<int>(all.size()) > k) all.resize(k);
    return all;
  }

  // Binary serialization. Format:
  //   magic ('SCDB') | version (uint32) | n_anchors (uint64) | sc_rows (uint32) | sc_cols (uint32)
  //   per anchor: id (int32) | pos (3 * double) | sc_desc (sc_rows*sc_cols * double) | ringkey
  bool save(const std::string& path) const {
    std::ofstream f(path, std::ios::binary);
    if (!f) return false;
    const char magic[4] = {'S','C','D','B'};
    uint32_t version = 2;  // v2 adds sc_radius after sc_cols
    uint64_t n = anchors_.size();
    uint32_t sc_rows = anchors_.empty() ? 0 : static_cast<uint32_t>(anchors_[0].sc_desc.rows());
    uint32_t sc_cols = anchors_.empty() ? 0 : static_cast<uint32_t>(anchors_[0].sc_desc.cols());
    double sc_radius = sc_radius_;
    f.write(magic, 4);
    f.write(reinterpret_cast<const char*>(&version), 4);
    f.write(reinterpret_cast<const char*>(&n), 8);
    f.write(reinterpret_cast<const char*>(&sc_rows), 4);
    f.write(reinterpret_cast<const char*>(&sc_cols), 4);
    f.write(reinterpret_cast<const char*>(&sc_radius), 8);
    for (const auto& a : anchors_) {
      int32_t id = a.id;
      double pos[3] = {a.position(0), a.position(1), a.position(2)};
      f.write(reinterpret_cast<const char*>(&id), 4);
      f.write(reinterpret_cast<const char*>(pos), 24);
      // SC desc - column-major because Eigen defaults to column-major.
      for (int c = 0; c < a.sc_desc.cols(); ++c) {
        for (int r = 0; r < a.sc_desc.rows(); ++r) {
          double v = a.sc_desc(r, c);
          f.write(reinterpret_cast<const char*>(&v), 8);
        }
      }
      uint32_t rk_size = static_cast<uint32_t>(a.ringkey.size());
      f.write(reinterpret_cast<const char*>(&rk_size), 4);
      f.write(reinterpret_cast<const char*>(a.ringkey.data()), 4 * rk_size);
    }
    return true;
  }

  bool load(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return false;
    char magic[4];
    uint32_t version, sc_rows, sc_cols;
    uint64_t n;
    f.read(magic, 4);
    if (magic[0] != 'S' || magic[1] != 'C' || magic[2] != 'D' || magic[3] != 'B') return false;
    f.read(reinterpret_cast<char*>(&version), 4);
    f.read(reinterpret_cast<char*>(&n), 8);
    f.read(reinterpret_cast<char*>(&sc_rows), 4);
    f.read(reinterpret_cast<char*>(&sc_cols), 4);
    sc_radius_ = 0.0;
    if (version >= 2) {
      f.read(reinterpret_cast<char*>(&sc_radius_), 8);
    }
    anchors_.clear();
    anchors_.reserve(n);
    for (uint64_t i = 0; i < n; ++i) {
      Anchor a;
      int32_t id;
      double pos[3];
      f.read(reinterpret_cast<char*>(&id), 4);
      f.read(reinterpret_cast<char*>(pos), 24);
      a.id = id;
      a.position = Eigen::Vector3d(pos[0], pos[1], pos[2]);
      a.sc_desc.resize(sc_rows, sc_cols);
      for (uint32_t c = 0; c < sc_cols; ++c) {
        for (uint32_t r = 0; r < sc_rows; ++r) {
          double v;
          f.read(reinterpret_cast<char*>(&v), 8);
          a.sc_desc(r, c) = v;
        }
      }
      uint32_t rk_size;
      f.read(reinterpret_cast<char*>(&rk_size), 4);
      a.ringkey.resize(rk_size);
      f.read(reinterpret_cast<char*>(a.ringkey.data()), 4 * rk_size);
      anchors_.push_back(std::move(a));
    }
    return static_cast<bool>(f);
  }

  size_t size() const { return anchors_.size(); }

private:
  SCManager sc_;
  std::vector<Anchor> anchors_;
  double sc_radius_ = 0.0;  // saved + restored across save/load so query uses build-time radius
};

// ----- TEASER refinement at one anchor ---------------------------------------------------------
//
// Same pipeline as local_refinement, but we don't have a 6DOF rough_pose - just an XY anchor
// position and a yaw estimate from SC. Reconstruct a 4x4 rough_pose from those and run.

struct VerifyResult {
  bool ok = false;
  Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
  int n_correspondences = 0;
  int n_inliers = 0;
  double inlier_ratio = 0.0;
  double rmse = 0.0;
};

VerifyResult verifyOneCandidate(const Config& cfg, const teaser::PointCloud& map,
                                const teaser::PointCloud& scan, const Match& cand) {
  VerifyResult vr;

  // Build a rough pose: anchor XY/Z + SC's yaw, roll/pitch zero.
  Eigen::Matrix4d T_rough = Eigen::Matrix4d::Identity();
  T_rough.block<3, 3>(0, 0) = eulerRPYtoRot(0, 0, cand.yaw_rad);
  T_rough.block<3, 1>(0, 3) = cand.anchor_position;

  // 1. Crop local map around anchor.
  auto local = extractLocalMap(map, cand.anchor_position, cfg.local_radius);
  if (local.size() < 50) return vr;

  // 2. Pre-transform scan into the map frame using the rough pose.
  auto scan_init = applyTransform(scan, T_rough);

  // 3. Downsample.
  auto local_ds = voxelDownsample(local, cfg.voxel_size_map);
  auto scan_ds = voxelDownsample(scan_init, cfg.voxel_size_scan);

  // 4. FPFH features.
  teaser::FPFHEstimation fpfh;
  auto scan_features = fpfh.computeFPFHFeatures(scan_ds, cfg.fpfh_normal_radius, cfg.fpfh_radius);
  auto map_features = fpfh.computeFPFHFeatures(local_ds, cfg.fpfh_normal_radius, cfg.fpfh_radius);

  // 5. Match.
  teaser::Matcher matcher;
  auto correspondences = matcher.calculateCorrespondences(
      scan_ds, local_ds, *scan_features, *map_features,
      false, true, true, 0.95);
  vr.n_correspondences = static_cast<int>(correspondences.size());
  if (vr.n_correspondences < 6) return vr;

  if (static_cast<int>(correspondences.size()) > cfg.max_correspondences) {
    correspondences.resize(cfg.max_correspondences);
  }

  // 6. TEASER (fresh solver per candidate - avoids inlier_graph_ / Quatro-static state leak).
  auto params = makeTeaserParams(cfg);
  teaser::RobustRegistrationSolver solver(params);
  solver.solve(scan_ds, local_ds, correspondences);
  auto sol = solver.getSolution();
  if (!sol.valid) return vr;

  // 7. Compose.
  Eigen::Matrix4d T_corr = Eigen::Matrix4d::Identity();
  T_corr.block<3, 3>(0, 0) = sol.rotation;
  T_corr.block<3, 1>(0, 3) = sol.translation;
  vr.T = T_corr * T_rough;

  // 8. Quality.
  vr.n_inliers = static_cast<int>(solver.getTranslationInliers().size());
  vr.inlier_ratio = correspondences.empty()
                        ? 0.0
                        : static_cast<double>(vr.n_inliers) / correspondences.size();
  double sq_sum = 0.0;
  for (int idx : solver.getInputOrderedTranslationInliers()) {
    Eigen::Vector3d a(scan_ds[correspondences[idx].first].x,
                      scan_ds[correspondences[idx].first].y,
                      scan_ds[correspondences[idx].first].z);
    Eigen::Vector3d b(local_ds[correspondences[idx].second].x,
                      local_ds[correspondences[idx].second].y,
                      local_ds[correspondences[idx].second].z);
    Eigen::Vector3d r = (sol.rotation * a + sol.translation) - b;
    sq_sum += r.squaredNorm();
  }
  vr.rmse = vr.n_inliers > 0 ? std::sqrt(sq_sum / vr.n_inliers) : 0.0;
  vr.ok = vr.n_inliers >= cfg.verify_min_inliers;
  return vr;
}

// ----- CLI ---------------------------------------------------------------------------------------

void printUsage(const char* a) {
  std::cout
      << "Usage:\n"
      << "  " << a << " --mode build    --map MAP.ply [--db OUT.bin] [--anchor-step S]\n"
      << "                                [--sc-radius R] [--sc-downsample V]\n"
      << "  " << a << " --mode query    --db DB.bin --scan SCAN.ply [--top-k K]\n"
      << "  " << a << " --mode localize --map MAP.ply --db DB.bin --scan SCAN.ply\n"
      << "                                [--top-k K] [--algo Quatro|TEASER] [--output OUT.txt]\n"
      << "                                [+ all local_refinement-style refinement options]\n";
}

bool parseArgs(int argc, char** argv, Config& cfg) {
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    auto need = [&](const char* opt) {
      if (i + 1 >= argc) { std::cerr << opt << " needs a value\n"; std::exit(1); }
      return std::string(argv[++i]);
    };
    if (a == "--mode") cfg.mode = need("--mode");
    else if (a == "--map") cfg.map_path = need("--map");
    else if (a == "--scan") cfg.scan_path = need("--scan");
    else if (a == "--db") cfg.db_path = need("--db");
    else if (a == "--output") cfg.output_path = need("--output");
    else if (a == "--anchor-step") cfg.anchor_step = std::stod(need("--anchor-step"));
    else if (a == "--sc-radius") cfg.sc_radius = std::stod(need("--sc-radius"));
    else if (a == "--sc-downsample") cfg.sc_downsample = std::stod(need("--sc-downsample"));
    else if (a == "--top-k") cfg.top_k = std::stoi(need("--top-k"));
    else if (a == "--algo") cfg.algorithm = need("--algo");
    else if (a == "--local-radius") cfg.local_radius = std::stod(need("--local-radius"));
    else if (a == "--voxel-map") cfg.voxel_size_map = std::stod(need("--voxel-map"));
    else if (a == "--voxel-scan") cfg.voxel_size_scan = std::stod(need("--voxel-scan"));
    else if (a == "--fpfh-normal-r") cfg.fpfh_normal_radius = std::stod(need("--fpfh-normal-r"));
    else if (a == "--fpfh-r") cfg.fpfh_radius = std::stod(need("--fpfh-r"));
    else if (a == "--noise-bound") cfg.noise_bound = std::stod(need("--noise-bound"));
    else if (a == "--max-corres") cfg.max_correspondences = std::stoi(need("--max-corres"));
    else if (a == "--min-inliers") cfg.verify_min_inliers = std::stoi(need("--min-inliers"));
    else if (a == "-h" || a == "--help") { printUsage(argv[0]); std::exit(0); }
    else { std::cerr << "Unknown arg: " << a << "\n"; printUsage(argv[0]); std::exit(1); }
  }
  if (cfg.mode != "build" && cfg.mode != "query" && cfg.mode != "localize") {
    std::cerr << "--mode must be build / query / localize\n";
    printUsage(argv[0]);
    return false;
  }
  return true;
}

// ----- Mode dispatch -----------------------------------------------------------------------------

int doBuild(const Config& cfg) {
  if (cfg.map_path.empty()) { std::cerr << "--map required for build\n"; return 1; }
  teaser::PointCloud map_cloud;
  if (!loadCloud(cfg.map_path, map_cloud)) {
    std::cerr << "Failed to read map: " << cfg.map_path << "\n"; return 1;
  }
  std::cout << "Loaded map: " << map_cloud.size() << " points.\n";

  SCDatabase db;
  auto t0 = std::chrono::steady_clock::now();
  db.build(map_cloud, cfg.anchor_step, cfg.sc_radius, cfg.sc_downsample, cfg.anchor_z_max_dev);
  auto t1 = std::chrono::steady_clock::now();
  std::cout << "Built " << db.size() << " anchors in "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() / 1000.0
            << " s.\n";

  if (!db.save(cfg.db_path)) {
    std::cerr << "Failed to write " << cfg.db_path << "\n"; return 1;
  }
  std::cout << "Wrote " << cfg.db_path << "\n";
  return 0;
}

int doQuery(const Config& cfg) {
  if (cfg.scan_path.empty()) { std::cerr << "--scan required for query\n"; return 1; }
  SCDatabase db;
  if (!db.load(cfg.db_path)) { std::cerr << "Failed to load " << cfg.db_path << "\n"; return 1; }
  std::cout << "Loaded DB: " << db.size() << " anchors.\n";

  teaser::PointCloud scan_cloud;
  if (!loadCloud(cfg.scan_path, scan_cloud)) {
    std::cerr << "Failed to read scan: " << cfg.scan_path << "\n"; return 1;
  }
  std::cout << "Loaded scan: " << scan_cloud.size() << " points.\n";

  auto t0 = std::chrono::steady_clock::now();
  auto matches = db.queryTopK(scan_cloud, cfg.sc_downsample, cfg.top_k);
  auto t1 = std::chrono::steady_clock::now();
  std::cout << "SC retrieval: " << matches.size() << " candidates in "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << " ms.\n";

  std::cout << "\nTop-" << matches.size() << " candidates:\n";
  std::cout << "  rank | anchor | sc_dist |          position           | yaw (deg)\n";
  for (size_t i = 0; i < matches.size(); ++i) {
    const auto& m = matches[i];
    std::cout << "  " << (i + 1) << "    | " << m.anchor_id << "    | "
              << m.sc_distance << "  | ("
              << m.anchor_position(0) << ", "
              << m.anchor_position(1) << ", "
              << m.anchor_position(2) << ") | "
              << m.yaw_rad * 180.0 / M_PI << "\n";
  }
  return 0;
}

int doLocalize(const Config& cfg) {
  if (cfg.map_path.empty() || cfg.scan_path.empty()) {
    std::cerr << "localize requires --map, --scan, --db\n"; return 1;
  }
  teaser::PointCloud map_cloud, scan_cloud;
  if (!loadCloud(cfg.map_path, map_cloud)) {
    std::cerr << "Failed to read map\n"; return 1;
  }
  if (!loadCloud(cfg.scan_path, scan_cloud)) {
    std::cerr << "Failed to read scan\n"; return 1;
  }
  SCDatabase db;
  if (!db.load(cfg.db_path)) { std::cerr << "Failed to load DB\n"; return 1; }
  std::cout << "Map: " << map_cloud.size() << " pts, scan: " << scan_cloud.size()
            << " pts, DB: " << db.size() << " anchors.\n";

  auto t0 = std::chrono::steady_clock::now();

  // Stage A: SC retrieval.
  auto matches = db.queryTopK(scan_cloud, cfg.sc_downsample, cfg.top_k);
  auto t_sc = std::chrono::steady_clock::now();
  std::cout << "SC retrieval: " << matches.size() << " candidates in "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t_sc - t0).count()
            << " ms.\n";

  // Stage B: verify each candidate.
  int best_idx = -1;
  VerifyResult best_vr;
  for (size_t i = 0; i < matches.size(); ++i) {
    const auto& m = matches[i];
    std::cout << "\nCandidate " << (i + 1) << "/" << matches.size()
              << " (anchor " << m.anchor_id << ", sc_dist " << m.sc_distance
              << ", yaw " << m.yaw_rad * 180.0 / M_PI << " deg):\n";
    auto vr = verifyOneCandidate(cfg, map_cloud, scan_cloud, m);
    std::cout << "  correspondences: " << vr.n_correspondences
              << ", inliers: " << vr.n_inliers
              << " (" << vr.inlier_ratio * 100 << "%)"
              << ", rmse: " << vr.rmse << " m"
              << ", ok: " << (vr.ok ? "yes" : "no") << "\n";

    // Winner = highest inlier count among "ok" candidates; if none are ok, highest inlier_ratio.
    bool prefer = false;
    if (!best_vr.ok && vr.ok) prefer = true;
    else if (best_vr.ok == vr.ok) {
      if (vr.n_inliers > best_vr.n_inliers) prefer = true;
    }
    if (prefer || best_idx < 0) {
      best_idx = static_cast<int>(i);
      best_vr = vr;
    }
  }
  auto t1 = std::chrono::steady_clock::now();
  double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() / 1000.0;

  if (best_idx < 0 || !best_vr.ok) {
    std::cerr << "\nNo candidate passed verification (need >= " << cfg.verify_min_inliers
              << " inliers).\n";
    return 2;
  }

  std::cout << "\n=====================================\n";
  std::cout << "       Global Localization Result      \n";
  std::cout << "=====================================\n";
  std::cout << "Winning anchor: " << matches[best_idx].anchor_id
            << " (rank " << best_idx + 1 << " in SC retrieval)\n";
  std::cout << "Refined pose:\n" << best_vr.T << "\n";
  std::cout << "Correspondences: " << best_vr.n_correspondences << "\n";
  std::cout << "Inliers:         " << best_vr.n_inliers
            << " (" << best_vr.inlier_ratio * 100 << "%)\n";
  std::cout << "Inlier RMSE:     " << best_vr.rmse << " m\n";
  std::cout << "Total elapsed:   " << elapsed << " s\n";

  writePose(cfg.output_path, best_vr.T, matches[best_idx].anchor_id,
            best_vr.inlier_ratio, best_vr.rmse, elapsed);
  std::cout << "Wrote " << cfg.output_path << "\n";
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc == 1) {
    printUsage(argv[0]);
    std::cout << "\nThis example needs real data. See README.md for a recipe using the\n"
              << "3DMatch sample (cloud_bin_36.ply as map, cloud_bin_2.ply as scan).\n";
    return 0;
  }
  Config cfg;
  if (!parseArgs(argc, argv, cfg)) return 1;
  if (cfg.mode == "build") return doBuild(cfg);
  if (cfg.mode == "query") return doQuery(cfg);
  return doLocalize(cfg);
}
