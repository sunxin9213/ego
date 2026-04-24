#ifndef _GRID_MAP_H
#define _GRID_MAP_H

#include <Eigen/Eigen>
#include <Eigen/StdVector>
#include <cv_bridge/cv_bridge.h>
#include <geometry_msgs/PoseStamped.h>
#include <iostream>
#include <random>
#include <nav_msgs/Odometry.h>
#include <queue>
#include <ros/ros.h>
#include <tuple>
#include <visualization_msgs/Marker.h>

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl_conversions/pcl_conversions.h>

#include <message_filters/subscriber.h>
#include <message_filters/sync_policies/approximate_time.h>
#include <message_filters/sync_policies/exact_time.h>
#include <message_filters/time_synchronizer.h>

#include <plan_env/raycast.h>

/**
 * @brief 对数几率函数
 * @param x 输入概率值(0到1之间)
 * @return 对数几率值
 */
#define logit(x) (log((x) / (1 - (x))))

using namespace std;

/**
 * @brief 体素哈希函数模板
 *
 * 用于将矩阵(如Eigen矩阵)转换为哈希值,用于unordered_map等哈希表中。
 * @tparam T 矩阵类型
 */
template <typename T>
struct matrix_hash : std::unary_function<T, size_t> {
  std::size_t operator()(T const& matrix) const {
    size_t seed = 0;
    for (size_t i = 0; i < matrix.size(); ++i) {
      auto elem = *(matrix.data() + i);
      seed ^= std::hash<typename T::Scalar>()(elem) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
  }
};

/**
 * @brief 地图构建参数结构体
 *
 * 存储3D地图构建所需的所有参数,包括地图属性、相机参数、深度图像过滤、
 * 射线投射、局部地图更新和可视化等配置。
 */
struct MappingParameters {

  /* map properties */
  // 地图原点和尺寸
  Eigen::Vector3d map_origin_, map_size_;
  // 地图在位置空间的最小和最大边界
  Eigen::Vector3d map_min_boundary_, map_max_boundary_;
  // 地图在索引空间的体素数量
  Eigen::Vector3i map_voxel_num_;
  // 局部更新范围
  Eigen::Vector3d local_update_range_;
  // 地图分辨率及其倒数
  double resolution_, resolution_inv_;
  // 障碍物膨胀距离
  double obstacles_inflation_;
  // 地图坐标系ID
  string frame_id_;
  // 姿态数据类型(1:PoseStamped, 2:Odometry)
  int pose_type_;

  /* camera parameters */
  // 相机内参:主点和焦距
  double cx_, cy_, fx_, fy_;

  /* depth image projection filtering */
  // 深度图像投影过滤参数
  double depth_filter_maxdist_, depth_filter_mindist_, depth_filter_tolerance_;  // 最大/最小/容差距离
  int depth_filter_margin_;  // 深度过滤边距
  bool use_depth_filter_;    // 是否使用深度过滤
  double k_depth_scaling_factor_;  // 深度缩放因子
  int skip_pixel_;  // 跳过的像素数

  /* raycasting */
  // 占用概率参数(命中、未命中、最小、最大、占据)
  double p_hit_, p_miss_, p_min_, p_max_, p_occ_;
  // 占用概率的对数形式
  double prob_hit_log_, prob_miss_log_, clamp_min_log_, clamp_max_log_,
      min_occupancy_log_;
  // 射线长度范围
  double min_ray_length_, max_ray_length_;

  /* local map update and clear */
  // 局部地图边距
  int local_map_margin_;

  /* visualization and computation time display */
  // 可视化参数:截断高度、虚拟天花板高度、地面高度
  double visualization_truncate_height_, virtual_ceil_height_, ground_height_;
  bool show_occ_time_;  // 是否显示占用更新时间

  /* active mapping */
  // 未知区域标志
  double unknown_flag_;
};

/**
 * @brief 地图构建数据结构
 *
 * 存储地图构建过程中的中间数据,包括占用率缓冲区、相机位姿、深度图像
 * 以及用于加速射线投射的标志缓冲区等。
 */
struct MappingData {
  // main map data, occupancy of each voxel and Euclidean distance

  // 占用率缓冲区(存储每个体素的占用概率对数值)
  std::vector<double> occupancy_buffer_;
  // 膨胀占用缓冲区(标记障碍物)
  std::vector<char> occupancy_buffer_inflate_;

  // camera position and pose data

  // 相机位置和姿态(当前和上一次)
  Eigen::Vector3d camera_pos_, last_camera_pos_;
  Eigen::Quaterniond camera_q_, last_camera_q_;

  // depth image data

  // 深度图像(当前和上一次)
  cv::Mat depth_image_, last_depth_image_;
  // 图像计数
  int image_cnt_;

  // flags of map state

  // 地图状态标志:需要更新、已局部更新
  bool occ_need_update_, local_updated_;
  // 是否有第一帧深度图
  bool has_first_depth_;
  // 是否有里程计和点云数据
  bool has_odom_, has_cloud_;

  // depth image projected point cloud

  // 深度图像投影得到的点云
  vector<Eigen::Vector3d> proj_points_;
  // 投影点数量
  int proj_points_cnt;

  // flag buffers for speeding up raycasting

  // 射线投射加速标志缓冲区
  vector<short> count_hit_, count_hit_and_miss_;
  vector<char> flag_traverse_, flag_rayend_;
  // 射线编号
  char raycast_num_;
  // 体素缓存队列
  queue<Eigen::Vector3i> cache_voxel_;

  // range of updating grid

  // 局部更新边界(最小和最大索引)
  Eigen::Vector3i local_bound_min_, local_bound_max_;

  // computation time

  // 融合时间和最大融合时间
  double fuse_time_, max_fuse_time_;
  // 更新次数
  int update_num_;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

/**
 * @brief 3D栅格地图类
 *
 * 基于八叉树或体素哈希的3D栅格地图实现,用于无人机路径规划中的环境表示。
 * 主要功能包括:
 * - 基于深度图像的地图更新
 * - 障碍物膨胀
 * - 射线投射更新占用率
 * - 地图可视化
 */
class GridMap {
public:
  GridMap() {}
  ~GridMap() {}

  /**
   * @brief 姿态数据类型枚举
   */
  enum { POSE_STAMPED = 1, ODOMETRY = 2, INVALID_IDX = -10000 };

  // occupancy map management

  /**
   * @brief 重置占用缓冲区
   */
  void resetBuffer();

  /**
   * @brief 在指定范围内重置占用缓冲区
   * @param min 最小边界位置
   * @param max 最大边界位置
   */
  void resetBuffer(Eigen::Vector3d min, Eigen::Vector3d max);

  /**
   * @brief 将世界坐标转换为体素索引
   * @param pos 世界坐标位置
   * @param id 输出的体素索引
   */
  inline void posToIndex(const Eigen::Vector3d& pos, Eigen::Vector3i& id);

  /**
   * @brief 将体素索引转换为世界坐标
   * @param id 体素索引
   * @param pos 输出的世界坐标
   */
  inline void indexToPos(const Eigen::Vector3i& id, Eigen::Vector3d& pos);

  /**
   * @brief 将体素索引转换为内存地址
   * @param id 体素索引
   * @return 内存地址
   */
  inline int toAddress(const Eigen::Vector3i& id);

  /**
   * @brief 将坐标分量转换为内存地址
   * @param x X坐标
   * @param y Y坐标
   * @param z Z坐标
   * @return 内存地址
   */
  inline int toAddress(int& x, int& y, int& z);

  /**
   * @brief 检查位置是否在地图范围内
   * @param pos 世界坐标位置
   * @return 是否在地图内
   */
  inline bool isInMap(const Eigen::Vector3d& pos);

  /**
   * @brief 检查索引是否在地图范围内
   * @param idx 体素索引
   * @return 是否在地图内
   */
  inline bool isInMap(const Eigen::Vector3i& idx);

  /**
   * @brief 设置位置占用率
   * @param pos 世界坐标位置
   * @param occ 占用值(0或1)
   */
  inline void setOccupancy(Eigen::Vector3d pos, double occ = 1);

  /**
   * @brief 设置位置为已占用(障碍物)
   * @param pos 世界坐标位置
   */
  inline void setOccupied(Eigen::Vector3d pos);

  /**
   * @brief 获取位置占用状态
   * @param pos 世界坐标位置
   * @return 占用状态(1:占用, 0:空闲, -1:未知)
   */
  inline int getOccupancy(Eigen::Vector3d pos);

  /**
   * @brief 获取体素索引的占用状态
   * @param id 体素索引
   * @return 占用状态(1:占用, 0:空闲, -1:未知)
   */
  inline int getOccupancy(Eigen::Vector3i id);

  /**
   * @brief 获取膨胀后的占用状态
   * @param pos 世界坐标位置
   * @return 膨胀占用状态(1:占用, 0:空闲)
   */
  inline int getInflateOccupancy(Eigen::Vector3d pos);

  /**
   * @brief 限制索引在地图范围内
   * @param id 输入/输出的体素索引
   */
  inline void boundIndex(Eigen::Vector3i& id);

  /**
   * @brief 检查索引是否为未知区域
   * @param id 体素索引
   * @return 是否未知
   */
  inline bool isUnknown(const Eigen::Vector3i& id);

  /**
   * @brief 检查位置是否为未知区域
   * @param pos 世界坐标位置
   * @return 是否未知
   */
  inline bool isUnknown(const Eigen::Vector3d& pos);

  /**
   * @brief 检查索引是否为已知空闲区域
   * @param id 体素索引
   * @return 是否空闲
   */
  inline bool isKnownFree(const Eigen::Vector3i& id);

  /**
   * @brief 检查索引是否为已知占用区域
   * @param id 体素索引
   * @return 是否占用
   */
  inline bool isKnownOccupied(const Eigen::Vector3i& id);

  /**
   * @brief 初始化地图
   * @param nh ROS节点句柄
   */
  void initMap(ros::NodeHandle& nh);

  /**
   * @brief 发布地图(原始占用率)
   */
  void publishMap();

  /**
   * @brief 发布膨胀后的地图
   * @param all_info 是否显示完整信息
   */
  void publishMapInflate(bool all_info = false);

  /**
   * @brief 发布未知区域
   */
  void publishUnknown();

  /**
   * @brief 发布深度图像
   */
  void publishDepth();

  /**
   * @brief 是否有深度观测
   * @return 是否有深度数据
   */
  bool hasDepthObservation();

  /**
   * @brief 里程计是否有效
   * @return 里程计是否有效
   */
  bool odomValid();

  /**
   * @brief 获取地图区域
   * @param ori 输出原点
   * @param size 输出尺寸
   */
  void getRegion(Eigen::Vector3d& ori, Eigen::Vector3d& size);

  /**
   * @brief 获取地图分辨率
   * @return 分辨率值
   */
  inline double getResolution();

  /**
   * @brief 获取地图原点
   * @return 地图原点坐标
   */
  Eigen::Vector3d getOrigin();

  /**
   * @brief 获取体素总数
   * @return 体素数量
   */
  int getVoxelNum();

  typedef std::shared_ptr<GridMap> Ptr;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  // 地图构建参数
  MappingParameters mp_;
  // 地图构建数据
  MappingData md_;

  // get depth image and camera pose
  /**
   * @brief 深度图像和姿态同步回调
   * @param img 深度图像消息
   * @param pose 姿态消息
   */
  void depthPoseCallback(const sensor_msgs::ImageConstPtr& img,
                         const geometry_msgs::PoseStampedConstPtr& pose);

  /**
   * @brief 深度图像和里程计同步回调
   * @param img 深度图像消息
   * @param odom 里程计消息
   */
  void depthOdomCallback(const sensor_msgs::ImageConstPtr& img, const nav_msgs::OdometryConstPtr& odom);

  /**
   * @brief 点云回调
   * @param img 点云消息
   */
  void cloudCallback(const sensor_msgs::PointCloud2ConstPtr& img);

  /**
   * @brief 里程计回调
   * @param odom 里程计消息
   */
  void odomCallback(const nav_msgs::OdometryConstPtr& odom);

  // update occupancy by raycasting

  /**
   * @brief 占用率更新回调(定时器触发)
   * @param event 定时器事件
   */
  void updateOccupancyCallback(const ros::TimerEvent& /*event*/);

  /**
   * @brief 可视化回调(定时器触发)
   * @param event 定时器事件
   */
  void visCallback(const ros::TimerEvent& /*event*/);

  // main update process

  /**
   * @brief 投影深度图像
   */
  void projectDepthImage();

  /**
   * @brief 射线投射处理
   */
  void raycastProcess();

  /**
   * @brief 清理和膨胀局部地图
   */
  void clearAndInflateLocalMap();

  /**
   * @brief 膨胀点(在周围标记障碍物)
   * @param pt 中心点索引
   * @param step 膨胀步数
   * @param pts 输出的膨胀点列表
   */
  inline void inflatePoint(const Eigen::Vector3i& pt, int step, vector<Eigen::Vector3i>& pts);

  /**
   * @brief 设置缓存占用
   * @param pos 位置
   * @param occ 占用值
   * @return 设置结果
   */
  int setCacheOccupancy(Eigen::Vector3d pos, int occ);

  /**
   * @brief 获取地图中最近的点
   * @param pt 目标点
   * @param camera_pt 相机位置
   * @return 地图中最近的点
   */
  Eigen::Vector3d closetPointInMap(const Eigen::Vector3d& pt, const Eigen::Vector3d& camera_pt);

  // 消息过滤器同步策略别名
  typedef message_filters::sync_policies::ApproximateTime<sensor_msgs::Image, nav_msgs::Odometry>
      SyncPolicyImageOdom;
  typedef message_filters::sync_policies::ApproximateTime<sensor_msgs::Image, geometry_msgs::PoseStamped>
      SyncPolicyImagePose;
  typedef shared_ptr<message_filters::Synchronizer<SyncPolicyImagePose>> SynchronizerImagePose;
  typedef shared_ptr<message_filters::Synchronizer<SyncPolicyImageOdom>> SynchronizerImageOdom;

  // ROS节点句柄
  ros::NodeHandle node_;
  // 消息订阅者
  shared_ptr<message_filters::Subscriber<sensor_msgs::Image>> depth_sub_;
  shared_ptr<message_filters::Subscriber<geometry_msgs::PoseStamped>> pose_sub_;
  shared_ptr<message_filters::Subscriber<nav_msgs::Odometry>> odom_sub_;
  // 消息同步器
  SynchronizerImagePose sync_image_pose_;
  SynchronizerImageOdom sync_image_odom_;

  // 独立订阅者和发布者
  ros::Subscriber indep_cloud_sub_, indep_odom_sub_;
  ros::Publisher map_pub_, map_inf_pub_;
  ros::Publisher unknown_pub_;
  // 定时器
  ros::Timer occ_timer_, vis_timer_;

  // 随机噪声生成器
  uniform_real_distribution<double> rand_noise_;
  normal_distribution<double> rand_noise2_;
  default_random_engine eng_;
};

/* ============================== definition of inline function
 * ============================== */

inline int GridMap::toAddress(const Eigen::Vector3i& id) {
  return id(0) * mp_.map_voxel_num_(1) * mp_.map_voxel_num_(2) + id(1) * mp_.map_voxel_num_(2) + id(2);
}

inline int GridMap::toAddress(int& x, int& y, int& z) {
  return x * mp_.map_voxel_num_(1) * mp_.map_voxel_num_(2) + y * mp_.map_voxel_num_(2) + z;
}

inline void GridMap::boundIndex(Eigen::Vector3i& id) {
  Eigen::Vector3i id1;
  id1(0) = max(min(id(0), mp_.map_voxel_num_(0) - 1), 0);
  id1(1) = max(min(id(1), mp_.map_voxel_num_(1) - 1), 0);
  id1(2) = max(min(id(2), mp_.map_voxel_num_(2) - 1), 0);
  id = id1;
}

inline bool GridMap::isUnknown(const Eigen::Vector3i& id) {
  Eigen::Vector3i id1 = id;
  boundIndex(id1);
  return md_.occupancy_buffer_[toAddress(id1)] < mp_.clamp_min_log_ - 1e-3;
}

inline bool GridMap::isUnknown(const Eigen::Vector3d& pos) {
  Eigen::Vector3i idc;
  posToIndex(pos, idc);
  return isUnknown(idc);
}

inline bool GridMap::isKnownFree(const Eigen::Vector3i& id) {
  Eigen::Vector3i id1 = id;
  boundIndex(id1);
  int adr = toAddress(id1);

  // return md_.occupancy_buffer_[adr] >= mp_.clamp_min_log_ &&
  //     md_.occupancy_buffer_[adr] < mp_.min_occupancy_log_;
  return md_.occupancy_buffer_[adr] >= mp_.clamp_min_log_ && md_.occupancy_buffer_inflate_[adr] == 0;
}

inline bool GridMap::isKnownOccupied(const Eigen::Vector3i& id) {
  Eigen::Vector3i id1 = id;
  boundIndex(id1);
  int adr = toAddress(id1);

  return md_.occupancy_buffer_inflate_[adr] == 1;
}

inline void GridMap::setOccupied(Eigen::Vector3d pos) {
  if (!isInMap(pos)) return;

  Eigen::Vector3i id;
  posToIndex(pos, id);

  md_.occupancy_buffer_inflate_[id(0) * mp_.map_voxel_num_(1) * mp_.map_voxel_num_(2) +
                                id(1) * mp_.map_voxel_num_(2) + id(2)] = 1;
}

inline void GridMap::setOccupancy(Eigen::Vector3d pos, double occ) {
  if (occ != 1 && occ != 0) {
    cout << "occ value error!" << endl;
    return;
  }

  if (!isInMap(pos)) return;

  Eigen::Vector3i id;
  posToIndex(pos, id);

  md_.occupancy_buffer_[toAddress(id)] = occ;
}

inline int GridMap::getOccupancy(Eigen::Vector3d pos) {
  if (!isInMap(pos)) return -1;

  Eigen::Vector3i id;
  posToIndex(pos, id);

  return md_.occupancy_buffer_[toAddress(id)] > mp_.min_occupancy_log_ ? 1 : 0;
}

inline int GridMap::getInflateOccupancy(Eigen::Vector3d pos) {
  if (!isInMap(pos)) return -1;

  Eigen::Vector3i id;
  posToIndex(pos, id);

  return int(md_.occupancy_buffer_inflate_[toAddress(id)]);
}

inline int GridMap::getOccupancy(Eigen::Vector3i id) {
  if (id(0) < 0 || id(0) >= mp_.map_voxel_num_(0) || id(1) < 0 || id(1) >= mp_.map_voxel_num_(1) ||
      id(2) < 0 || id(2) >= mp_.map_voxel_num_(2))
    return -1;

  return md_.occupancy_buffer_[toAddress(id)] > mp_.min_occupancy_log_ ? 1 : 0;
}

inline bool GridMap::isInMap(const Eigen::Vector3d& pos) {
  if (pos(0) < mp_.map_min_boundary_(0) + 1e-4 || pos(1) < mp_.map_min_boundary_(1) + 1e-4 ||
      pos(2) < mp_.map_min_boundary_(2) + 1e-4) {
    // cout << "less than min range!" << endl;
    return false;
  }
  if (pos(0) > mp_.map_max_boundary_(0) - 1e-4 || pos(1) > mp_.map_max_boundary_(1) - 1e-4 ||
      pos(2) > mp_.map_max_boundary_(2) - 1e-4) {
    return false;
  }
  return true;
}

inline bool GridMap::isInMap(const Eigen::Vector3i& idx) {
  if (idx(0) < 0 || idx(1) < 0 || idx(2) < 0) {
    return false;
  }
  if (idx(0) > mp_.map_voxel_num_(0) - 1 || idx(1) > mp_.map_voxel_num_(1) - 1 ||
      idx(2) > mp_.map_voxel_num_(2) - 1) {
    return false;
  }
  return true;
}

inline void GridMap::posToIndex(const Eigen::Vector3d& pos, Eigen::Vector3i& id) {
  for (int i = 0; i < 3; ++i) id(i) = floor((pos(i) - mp_.map_origin_(i)) * mp_.resolution_inv_);
}

inline void GridMap::indexToPos(const Eigen::Vector3i& id, Eigen::Vector3d& pos) {
  for (int i = 0; i < 3; ++i) pos(i) = (id(i) + 0.5) * mp_.resolution_ + mp_.map_origin_(i);
}

inline void GridMap::inflatePoint(const Eigen::Vector3i& pt, int step, vector<Eigen::Vector3i>& pts) {
  int num = 0;

  /* ---------- + shape inflate ---------- */
  // for (int x = -step; x <= step; ++x)
  // {
  //   if (x == 0)
  //     continue;
  //   pts[num++] = Eigen::Vector3i(pt(0) + x, pt(1), pt(2));
  // }
  // for (int y = -step; y <= step; ++y)
  // {
  //   if (y == 0)
  //     continue;
  //   pts[num++] = Eigen::Vector3i(pt(0), pt(1) + y, pt(2));
  // }
  // for (int z = -1; z <= 1; ++z)
  // {
  //   pts[num++] = Eigen::Vector3i(pt(0), pt(1), pt(2) + z);
  // }

  /* ---------- all inflate ---------- */
  for (int x = -step; x <= step; ++x)
    for (int y = -step; y <= step; ++y)
      for (int z = -step; z <= step; ++z) {
        pts[num++] = Eigen::Vector3i(pt(0) + x, pt(1) + y, pt(2) + z);
      }
}

inline double GridMap::getResolution() { return mp_.resolution_; }

#endif