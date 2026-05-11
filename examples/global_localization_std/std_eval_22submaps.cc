// Phase 3 / STD vs SC head-to-head test.
//
// Setup: load 22 glim_ros submaps from a session dir. Each submap has a known
// T_world_origin (ground truth pose) in data.txt.
//
// Protocol (leave-one-out):
//   For each submap i in [0..N-1]:
//     1. Build a fresh STDescManager DB with all submaps except i.
//     2. Generate STDescs from submap i (the held-out query).
//     3. SearchLoop -> matched index k + (R, t) between submap i and submap k.
//     4. Compare predicted T_world_query against ground truth T_world_origin[i].
//
// Reports per-submap: did STD find a match? what was the geographic error?
//
// Inputs:
//   --session DIR     e.g. /home/steve/map_data/map_w2/20260511_110448
//   --n N             number of submaps to process (default = auto-detect)
//   --cloud-type local|levelled|world   (default: levelled)
//   --skip-i ID       just run this single submap as query (for debugging)

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

#include "third_party/std/STDesc.h"

namespace fs = std::filesystem;

namespace {

struct Pose {
  Eigen::Vector3d t;
  Eigen::Matrix3d R;
};

struct Submap {
  int id;
  Pose T_world_origin;
  pcl::PointCloud<pcl::PointXYZI>::Ptr cloud;
};

Eigen::Matrix4d toMatrix4(const Pose& p) {
  Eigen::Matrix4d M = Eigen::Matrix4d::Identity();
  M.block<3, 3>(0, 0) = p.R;
  M.block<3, 1>(0, 3) = p.t;
  return M;
}

// Parse the "T_world_origin:" 4x4 block from glim_ros's data.txt.
bool parsePoseFromDataTxt(const fs::path& data_txt, Pose& pose) {
  std::ifstream f(data_txt);
  if (!f) return false;
  std::string line;
  while (std::getline(f, line)) {
    if (line.find("T_world_origin") != std::string::npos) {
      Eigen::Matrix4d M;
      for (int r = 0; r < 4; ++r) {
        std::string row;
        if (!std::getline(f, row)) return false;
        std::istringstream iss(row);
        for (int c = 0; c < 4; ++c) {
          if (!(iss >> M(r, c))) return false;
        }
      }
      pose.R = M.block<3, 3>(0, 0);
      pose.t = M.block<3, 1>(0, 3);
      return true;
    }
  }
  return false;
}

// Load N*3*float32 from points_compact.bin (intensity = 0).
pcl::PointCloud<pcl::PointXYZI>::Ptr loadPointsCompact(const fs::path& bin) {
  auto out = std::make_shared<pcl::PointCloud<pcl::PointXYZI>>();
  std::ifstream f(bin, std::ios::binary);
  if (!f) return out;
  f.seekg(0, std::ios::end);
  size_t n_bytes = f.tellg();
  f.seekg(0, std::ios::beg);
  size_t n_pts = n_bytes / (3 * sizeof(float));
  out->reserve(n_pts);
  for (size_t i = 0; i < n_pts; ++i) {
    float xyz[3];
    f.read(reinterpret_cast<char*>(xyz), 12);
    pcl::PointXYZI p;
    p.x = xyz[0]; p.y = xyz[1]; p.z = xyz[2]; p.intensity = 0.0f;
    out->push_back(p);
  }
  return out;
}

// Apply just pitch+roll from a pose (to gravity-align the local cloud) and
// translate to the origin. Z becomes vertical; yaw is preserved.
pcl::PointCloud<pcl::PointXYZI>::Ptr levelCloud(
    const pcl::PointCloud<pcl::PointXYZI>::Ptr& src, const Eigen::Matrix3d& R) {
  // Decompose R = R_yaw * R_pitch * R_roll, then apply (R_pitch * R_roll) to points.
  double pitch = std::asin(-R(2, 0));
  double roll = std::atan2(R(2, 1), R(2, 2));
  Eigen::Matrix3d Ry, Rx;
  Ry << std::cos(pitch), 0, std::sin(pitch),
        0,               1, 0,
       -std::sin(pitch), 0, std::cos(pitch);
  Rx << 1, 0,              0,
        0, std::cos(roll), -std::sin(roll),
        0, std::sin(roll),  std::cos(roll);
  Eigen::Matrix3d R_level = Ry * Rx;

  auto out = std::make_shared<pcl::PointCloud<pcl::PointXYZI>>();
  out->reserve(src->size());
  for (const auto& p : src->points) {
    Eigen::Vector3d v(p.x, p.y, p.z);
    Eigen::Vector3d w = R_level * v;
    pcl::PointXYZI q;
    q.x = w(0); q.y = w(1); q.z = w(2); q.intensity = p.intensity;
    out->push_back(q);
  }
  return out;
}

// Transform a cloud by [R | t].
pcl::PointCloud<pcl::PointXYZI>::Ptr transformCloud(
    const pcl::PointCloud<pcl::PointXYZI>::Ptr& src, const Pose& pose) {
  auto out = std::make_shared<pcl::PointCloud<pcl::PointXYZI>>();
  out->reserve(src->size());
  for (const auto& p : src->points) {
    Eigen::Vector3d v(p.x, p.y, p.z);
    Eigen::Vector3d w = pose.R * v + pose.t;
    pcl::PointXYZI q;
    q.x = w(0); q.y = w(1); q.z = w(2); q.intensity = p.intensity;
    out->push_back(q);
  }
  return out;
}

// Pose composition for SE(3): given submap i and submap k poses in world,
// the relative transform R_ik, t_ik s.t. p_k_frame = R_ik * p_i_frame + t_ik is
//   R_ik = R_k^T * R_i
//   t_ik = R_k^T * (t_i - t_k)
// STD's SearchLoop returns (R, t) mapping query frame -> matched frame's frame.
Pose relativePose(const Pose& Tw_i, const Pose& Tw_k) {
  Pose rel;
  rel.R = Tw_k.R.transpose() * Tw_i.R;
  rel.t = Tw_k.R.transpose() * (Tw_i.t - Tw_k.t);
  return rel;
}

std::vector<Submap> loadAllSubmaps(const fs::path& session_dir,
                                   const std::string& cloud_kind) {
  std::vector<Submap> out;
  std::regex submap_re("[0-9]{6}");
  std::vector<fs::path> dirs;
  for (const auto& e : fs::directory_iterator(session_dir)) {
    if (e.is_directory() && std::regex_match(e.path().filename().string(), submap_re)) {
      dirs.push_back(e.path());
    }
  }
  std::sort(dirs.begin(), dirs.end());

  for (const auto& d : dirs) {
    Submap s;
    s.id = std::stoi(d.filename().string());
    if (!parsePoseFromDataTxt(d / "data.txt", s.T_world_origin)) {
      std::cerr << "skip " << d << ": cannot parse data.txt\n";
      continue;
    }
    auto local = loadPointsCompact(d / "points_compact.bin");
    if (local->empty()) {
      std::cerr << "skip " << d << ": no points\n";
      continue;
    }
    if (cloud_kind == "local") {
      s.cloud = local;
    } else if (cloud_kind == "levelled") {
      s.cloud = levelCloud(local, s.T_world_origin.R);
    } else if (cloud_kind == "world") {
      s.cloud = transformCloud(local, s.T_world_origin);
    } else {
      std::cerr << "unknown --cloud-type " << cloud_kind << "\n";
      std::exit(1);
    }
    out.push_back(std::move(s));
  }
  return out;
}

// One leave-one-out trial: rebuild DB without submap `held`, query with it.
struct TrialResult {
  int held_id = -1;
  int matched_db_index = -1;     // index into other_submaps, -1 if no match
  int matched_submap_id = -1;    // actual submap id, -1 if no match
  double std_score = 0.0;        // STD's icp/plane score (higher is better)
  Pose predicted_T_world_query;  // GT pose of matched + (R, t) = predicted query pose
  double trans_err_m = -1.0;
  double rot_err_deg = -1.0;
  double elapsed_s = 0.0;
};

TrialResult runLeaveOneOut(const std::vector<Submap>& submaps, int held_idx,
                           ConfigSetting cfg) {
  TrialResult r;
  r.held_id = submaps[held_idx].id;

  auto t0 = std::chrono::steady_clock::now();

  // Build a fresh DB with all submaps except `held_idx`.
  // skip_near_num=0 because we're using leave-one-out, not loop-closure-style.
  cfg.skip_near_num_ = 0;
  STDescManager mgr(cfg);

  std::vector<int> db_to_submap_idx;
  int total_db_descs = 0;
  for (size_t i = 0; i < submaps.size(); ++i) {
    if (static_cast<int>(i) == held_idx) continue;
    auto cloud = submaps[i].cloud;
    std::vector<STDesc> stds;
    auto copy = std::make_shared<pcl::PointCloud<pcl::PointXYZI>>(*cloud);
    mgr.GenerateSTDescs(copy, stds);
    mgr.AddSTDescs(stds);
    total_db_descs += stds.size();
    db_to_submap_idx.push_back(i);
  }

  // Query.
  const auto& q = submaps[held_idx];
  auto q_copy = std::make_shared<pcl::PointCloud<pcl::PointXYZI>>(*q.cloud);
  std::vector<STDesc> q_stds;
  mgr.GenerateSTDescs(q_copy, q_stds);

  std::pair<int, double> loop_result{-1, 0.0};
  std::pair<Eigen::Vector3d, Eigen::Matrix3d> loop_transform;
  std::vector<std::pair<STDesc, STDesc>> loop_pairs;
  mgr.SearchLoop(q_stds, loop_result, loop_transform, loop_pairs);

  auto t1 = std::chrono::steady_clock::now();
  r.elapsed_s = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1e6;

  r.matched_db_index = loop_result.first;
  r.std_score = loop_result.second;
  if (r.matched_db_index < 0 ||
      static_cast<size_t>(r.matched_db_index) >= db_to_submap_idx.size()) {
    return r;
  }
  int matched_orig_idx = db_to_submap_idx[r.matched_db_index];
  r.matched_submap_id = submaps[matched_orig_idx].id;

  // Predict query world pose from matched keyframe's GT pose plus STD's relative transform.
  //   STD returns (R_qm, t_qm) such that p_query = R_qm * p_matched + t_qm   (their convention,
  //   based on `loop_transform.first` = t and `loop_transform.second` = R; reading SearchLoop
  //   it computes the transform that aligns query points to matched - we use it as a relative).
  // T_world_query = T_world_matched * T_matched_query
  Pose rel;
  rel.t = loop_transform.first;
  rel.R = loop_transform.second;
  const Pose& Twm = submaps[matched_orig_idx].T_world_origin;
  r.predicted_T_world_query.R = Twm.R * rel.R;
  r.predicted_T_world_query.t = Twm.R * rel.t + Twm.t;

  // Compare with GT.
  Eigen::Vector3d t_err_vec = r.predicted_T_world_query.t - q.T_world_origin.t;
  r.trans_err_m = t_err_vec.norm();
  Eigen::Matrix3d R_err = q.T_world_origin.R.transpose() * r.predicted_T_world_query.R;
  double cos_th = std::fmax(-1.0, std::fmin(1.0, (R_err.trace() - 1.0) / 2.0));
  r.rot_err_deg = std::acos(cos_th) * 180.0 / M_PI;
  return r;
}

}  // namespace

int main(int argc, char** argv) {
  std::string session_dir;
  std::string cloud_kind = "levelled";
  int only_i = -1;
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    auto need = [&](const char* k) {
      if (i + 1 >= argc) { std::cerr << k << " needs a value\n"; std::exit(1); }
      return std::string(argv[++i]);
    };
    if (a == "--session") session_dir = need("--session");
    else if (a == "--cloud-type") cloud_kind = need("--cloud-type");
    else if (a == "--only") only_i = std::stoi(need("--only"));
    else if (a == "-h" || a == "--help") {
      std::cout << "Usage: " << argv[0]
                << " --session DIR [--cloud-type local|levelled|world] [--only ID]\n";
      return 0;
    } else { std::cerr << "unknown arg: " << a << "\n"; return 1; }
  }
  if (session_dir.empty()) {
    std::cerr << "need --session\n"; return 1;
  }

  std::cout << "Loading submaps from " << session_dir << " (cloud=" << cloud_kind << ")\n";
  auto submaps = loadAllSubmaps(session_dir, cloud_kind);
  std::cout << "Loaded " << submaps.size() << " submaps.\n\n";

  ConfigSetting cfg;
  init_default_parameters(cfg);
  // Tweaks for sparse Airy submaps (~2-8k points, ~40m scene scale).
  // Default STD is calibrated for dense Velodyne 64-line (100k+ pts per keyframe).
  cfg.voxel_size_ = 1.0;
  cfg.voxel_init_num_ = 5;          // 10 default; need fewer points to call a voxel "planar"
  cfg.maximum_corner_num_ = 500;
  cfg.corner_thre_ = 5.0;
  cfg.descriptor_min_len_ = 1.0;    // 2 default
  cfg.descriptor_max_len_ = 30.0;   // 50 default
  cfg.candidate_num_ = 50;
  cfg.icp_threshold_ = 0.0;
  cfg.rough_dis_threshold_ = 0.5;
  cfg.vertex_diff_threshold_ = 5.0;
  cfg.dis_threshold_ = 2.0;
  cfg.non_max_suppression_radius_ = 1.0;  // 2.0 default; allow denser corners

  std::cout << "STD config: voxel=" << cfg.voxel_size_
            << " tri_len=[" << cfg.descriptor_min_len_ << ", " << cfg.descriptor_max_len_ << "]"
            << " candidates=" << cfg.candidate_num_ << "\n\n";

  std::cout << "id | n_pts | GT pos      | matched | trans_err | rot_err | elapsed | result\n";
  std::cout << "---|-------|-------------|---------|-----------|---------|---------|-------\n";
  int n_pass = 0, n_total = 0;
  for (size_t i = 0; i < submaps.size(); ++i) {
    if (only_i >= 0 && submaps[i].id != only_i) continue;
    auto r = runLeaveOneOut(submaps, i, cfg);
    ++n_total;
    const auto& q = submaps[i];
    bool found = r.matched_submap_id >= 0;
    bool good = found && r.trans_err_m >= 0 && r.trans_err_m < 1.0 && r.rot_err_deg < 5.0;
    if (good) ++n_pass;
    char buf[256];
    if (found) {
      std::snprintf(buf, sizeof(buf),
                    "%02d | %5zu | (%+.2f,%+.2f) | %02d      | %.2fm    | %.1f°   | %.2fs   | %s\n",
                    q.id, q.cloud->size(), q.T_world_origin.t(0), q.T_world_origin.t(1),
                    r.matched_submap_id, r.trans_err_m, r.rot_err_deg, r.elapsed_s,
                    good ? "PASS" : "WRONG");
    } else {
      std::snprintf(buf, sizeof(buf),
                    "%02d | %5zu | (%+.2f,%+.2f) | --      |     -      |    -    | %.2fs   | NO_MATCH\n",
                    q.id, q.cloud->size(), q.T_world_origin.t(0), q.T_world_origin.t(1),
                    r.elapsed_s);
    }
    std::cout << buf;
  }
  std::cout << "---\n";
  std::cout << "PASS rate: " << n_pass << " / " << n_total << "\n";
  return 0;
}
