#ifndef _PLANNER_MANAGER_H_
#define _PLANNER_MANAGER_H_

#include <stdlib.h>

#include <bspline_opt/bspline_optimizer.h>
#include <bspline_opt/uniform_bspline.h>

#include <rebound_planner/DataDisp.h>

#include <plan_env/grid_map.h>

#include <plan_manage/plan_container.hpp>

#include <ros/ros.h>

#include <plan_manage/planning_visualization.h>

namespace rebound_planner {

/**
 * @brief 快速规划器管理器类
 *
 * 该类是整个规划系统的核心管理器,负责协调地图构建和路径规划的关键算法。
 * 主要功能包括:
 * - 运动学可行的轨迹重规划
 * - 基于拓扑的轨迹优化
 * - 轨迹碰撞检测
 * - 全局与局部轨迹管理
 */
class ReboundPlannerManager {
  // SECTION stable
public:
  ReboundPlannerManager();
  ~ReboundPlannerManager();

  /**
   * @brief 运动学可行的重规划接口
   * @param start_pt 起始位置
   * @param start_vel 起始速度
   * @param start_acc 起始加速度
   * @param end_pt 目标位置
   * @param end_vel 目标速度
   * @return 是否规划成功
   */
  bool kinodynamicReplan(Eigen::Vector3d start_pt, Eigen::Vector3d start_vel, Eigen::Vector3d start_acc,
                         Eigen::Vector3d end_pt, Eigen::Vector3d end_vel);

  /**
   * @brief 基于反弹策略的轨迹重规划
   * @param start_pt 起始位置
   * @param start_vel 起始速度
   * @param start_acc 起始加速度
   * @param end_pt 目标位置
   * @param end_vel 目标速度
   * @param flag_polyInit 是否使用多项式初始化
   * @param flag_randomPolyTraj 是否使用随机多项式轨迹
   * @return 是否规划成功
   */
  bool reboundReplan(Eigen::Vector3d start_pt, Eigen::Vector3d start_vel, Eigen::Vector3d start_acc,
                         Eigen::Vector3d end_pt, Eigen::Vector3d end_vel, bool flag_polyInit, bool flag_randomPolyTraj);

  /**
   * @brief 紧急停止规划
   * @param stop_pos 停止位置
   * @return 是否规划成功
   */
  bool EmergencyStop(Eigen::Vector3d stop_pos);

  /**
   * @brief 滑动窗口局部规划
   * @param local_target_pt 局部目标位置
   * @param local_target_vel 局部目标速度
   * @return 是否规划成功
   */
  bool slidingWindow(Eigen::Vector3d local_target_pt, Eigen::Vector3d local_target_vel);

  /**
   * @brief 从指定位置开始规划全局轨迹
   * @param start_pos 起始位置
   * @return 是否规划成功
   */
  bool planGlobalTraj(const Eigen::Vector3d& start_pos);

  /**
   * @brief 从指定状态开始规划全局轨迹
   * @param start_pos 起始位置
   * @param start_vel 起始速度
   * @param start_acc 起始加速度
   * @param end_pos 目标位置
   * @param end_vel 目标速度
   * @param end_acc 目标加速度
   * @return 是否规划成功
   */
  bool planGlobalTraj(const Eigen::Vector3d& start_pos, const Eigen::Vector3d& start_vel, const Eigen::Vector3d& start_acc,
                                        const Eigen::Vector3d& end_pos, const Eigen::Vector3d& end_vel, const Eigen::Vector3d& end_acc) ;

  /**
   * @brief 拓扑重规划
   * @param collide 是否发生碰撞
   * @return 是否规划成功
   */
  bool topoReplan(bool collide);

  /**
   * @brief 规划偏航角
   * @param start_yaw 起始偏航角
   */
  void planYaw(const Eigen::Vector3d& start_yaw);

  /**
   * @brief 初始化规划模块
   * @param nh ROS节点句柄
   * @param vis 可视化对象指针
   */
  void initPlanModules(ros::NodeHandle& nh, PlanningVisualization::Ptr vis = NULL);

  /**
   * @brief 检查轨迹碰撞
   * @param distance 输出碰撞距离
   * @return 是否发生碰撞
   */
  bool checkTrajCollision(double& distance);

  /**
   * @brief 使用膨胀距离检查轨迹碰撞
   * @param traj 待检查的轨迹
   * @return 是否发生碰撞
   */
  bool checkTrajCollisionInflate(UniformBspline &traj);

  // 规划参数
  PlanParameters pp_;
  // 局部轨迹数据
  LocalTrajData local_data_;
  // 全局轨迹数据
  GlobalTrajData global_data_;
  // 栅格地图指针
  GridMap::Ptr grid_map_;

private:
  /* main planning algorithms & modules */
  // 可视化对象
  PlanningVisualization::Ptr visualization_;

  // B样条优化器(反弹规划专用)
  BsplineOptimizer::Ptr bspline_optimizer_rebound_;

  // 连续失败次数计数
  int continous_failures_count_{0};

  /**
   * @brief 更新轨迹信息
   * @param position_traj 位置轨迹
   * @param time_now 当前时间
   */
  void updateTrajInfo(const UniformBspline& position_traj, const ros::Time time_now );

  // topology guided optimization - 拓扑引导优化

  /**
   * @brief 查找碰撞区间
   * @param colli_start 碰撞起点输出
   * @param colli_end 碰撞终点输出
   * @param start_pts 起点集合
   * @param end_pts 终点集合
   */
  void findCollisionRange(vector<Eigen::Vector3d>& colli_start, vector<Eigen::Vector3d>& colli_end,
                          vector<Eigen::Vector3d>& start_pts, vector<Eigen::Vector3d>& end_pts);

  /**
   * @brief 优化拓扑B样条轨迹
   * @param start_t 起始时间
   * @param duration 持续时间
   * @param guide_path 引导路径
   * @param traj_id 轨迹ID
   */
  void optimizeTopoBspline(double start_t, double duration, vector<Eigen::Vector3d> guide_path,
                           int traj_id);

  /**
   * @brief 重参数化局部轨迹(重载版本1)
   * @param start_t 起始时间
   * @param dt 时间步长(输出)
   * @param duration 持续时间(输出)
   * @return 控制点矩阵
   */
  Eigen::MatrixXd reparamLocalTraj(double start_t, double& dt, double& duration);

  /**
   * @brief 重参数化局部轨迹(重载版本2)
   * @param start_t 起始时间
   * @param duration 持续时间
   * @param seg_num 段数
   * @param dt 时间步长(输出)
   * @return 控制点矩阵
   */
  Eigen::MatrixXd reparamLocalTraj(double start_t, double duration, int seg_num, double& dt);

  /**
   * @brief 选择最佳轨迹
   * @param traj 输出的最佳轨迹
   */
  void selectBestTraj(UniformBspline& traj);

  /**
   * @brief 优化轨迹
   * @param best_traj 最佳轨迹
   * @param time_inc 时间增量
   */
  void refineTraj(UniformBspline& best_traj, double& time_inc);

  /**
   * @brief 重参数化B样条轨迹
   * @param bspline 原始B样条
   * @param start_end_derivative 起点和终点的导数
   * @param ratio 比例因子
   * @param ctrl_pts 控制点(输出)
   * @param dt 时间步长(输出)
   * @param time_inc 时间增量(输出)
   */
  void reparamBspline(UniformBspline& bspline, vector<Eigen::Vector3d>& start_end_derivative, double ratio, Eigen::MatrixXd& ctrl_pts, double& dt,
                      double& time_inc);

  /**
   * @brief 轨迹优化算法
   * @param traj 待优化的轨迹
   * @param start_end_derivative 起点和终点的导数
   * @param ratio 比例因子
   * @param ts 时间步长(输出)
   * @param optimal_control_points 最优控制点(输出)
   * @return 是否优化成功
   */
  bool refineTrajAlgo(UniformBspline& traj, vector<Eigen::Vector3d>& start_end_derivative, double ratio, double& ts, Eigen::MatrixXd& optimal_control_points);

  // !SECTION stable

  // SECTION developing

public:
  typedef unique_ptr<ReboundPlannerManager> Ptr;

  // !SECTION
};
}  // namespace rebound_planner

#endif