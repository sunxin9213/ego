#ifndef _REBO_REPLAN_FSM_H_
#define _REBO_REPLAN_FSM_H_

#include <Eigen/Eigen>
#include <algorithm>
#include <iostream>
#include <nav_msgs/Path.h>
#include <sensor_msgs/Imu.h>
#include <ros/ros.h>
#include <std_msgs/Empty.h>
#include <vector>
#include <visualization_msgs/Marker.h>

#include <bspline_opt/bspline_optimizer.h>
#include <plan_env/grid_map.h>
#include <rebound_planner/Bspline.h>
#include <rebound_planner/DataDisp.h>
#include <plan_manage/planner_manager.h>
#include <traj_utils/planning_visualization.h>

using std::vector;

namespace rebound_planner {

/**
 * @brief 反弹重规划有限状态机类
 *
 * 实现无人机轨迹规划的有限状态机,管理从目标设定到轨迹执行的完整流程。
 * 状态机包含以下状态:
 * - INIT: 初始化状态
 * - WAIT_TARGET: 等待目标状态
 * - GEN_NEW_TRAJ: 生成新轨迹状态
 * - REPLAN_TRAJ: 重规划轨迹状态
 * - EXEC_TRAJ: 执行轨迹状态
 * - EMERGENCY_STOP: 紧急停止状态
 */
class ReboReplanFSM {

private:
  /* ---------- flag ---------- */
  /**
   * @brief 有限状态机执行状态枚举
   */
  enum FSM_EXEC_STATE {
    INIT,               // 初始化
    WAIT_TARGET,        // 等待目标
    GEN_NEW_TRAJ,       // 生成新轨迹
    REPLAN_TRAJ,        // 重规划轨迹
    EXEC_TRAJ,          // 执行轨迹
    EMERGENCY_STOP      // 紧急停止
  };
  
  /**
   * @brief 目标类型枚举
   */
  enum TARGET_TYPE {
    MANUAL_TARGET = 1,    // 手动选择目标
    PRESET_TARGET = 2,    // 预设目标
    REFENCE_PATH = 3      // 参考路径
  };

  /* planning utils */
  // 规划器管理器指针
  ReboundPlannerManager::Ptr planner_manager_;
  // 可视化对象指针
  PlanningVisualization::Ptr visualization_;
  // 数据显示对象
  rebound_planner::DataDisp data_disp_;

  /* parameters */
  // 目标类型(1:手动选择, 2:预设)
  int target_type_;
  // 无重规划阈值、重规划阈值
  double no_replan_thresh_, replan_thresh_;
  // 航点数组(最多50个)
  double waypoints_[50][3];
  // 航点数量
  int waypoint_num_;
  // 规划视界距离和时间
  double planning_horizen_, planning_horizen_time_;

  /* planning data */
  // 触发标志、是否有目标、是否有里程计、是否有新目标
  bool trigger_, have_target_, have_odom_, have_new_target_;
  // 当前执行状态
  FSM_EXEC_STATE exec_state_;
  // 连续调用次数计数
  int continously_called_times_{0};

  // 里程计状态:位置、速度、加速度
  Eigen::Vector3d odom_pos_, odom_vel_, odom_acc_;
  // 里程计姿态(四元数)
  Eigen::Quaterniond odom_orient_;

  // 起始状态:初始点、起点、起始速度、起始加速度、起始偏航角
  Eigen::Vector3d init_pt_, start_pt_, start_vel_, start_acc_, start_yaw_;
  // 目标状态:目标点、目标速度
  Eigen::Vector3d end_pt_, end_vel_;
  // 局部目标状态:局部目标点、局部目标速度
  Eigen::Vector3d local_target_pt_, local_target_vel_;
  // 当前航点索引
  int current_wp_;

  // 是否正在从紧急状态中恢复
  bool flag_escape_emergency_;

  /* ROS utils */
  // ROS节点句柄
  ros::NodeHandle node_;
  // 定时器:执行定时器、安全检查定时器
  ros::Timer exec_timer_, safety_timer_;
  // 订阅者:航点订阅、里程计订阅
  ros::Subscriber waypoint_sub_, odom_sub_;
  // 发布者:重规划发布、新轨迹发布、B样条发布、数据显示发布
  ros::Publisher replan_pub_, new_pub_, bspline_pub_, data_disp_pub_;

  /* helper functions */
  /**
   * @brief 调用反弹重规划
   * @param flag_use_poly_init 是否使用多项式初始化
   * @param flag_randomPolyTraj 是否使用随机多项式轨迹
   * @return 是否重规划成功
   */
  bool callReboundReplan(bool flag_use_poly_init, bool flag_randomPolyTraj);

  /**
   * @brief 调用紧急停止
   * @param stop_pos 停止位置
   * @return 是否停止成功
   */
  bool callEmergencyStop(Eigen::Vector3d stop_pos);

  /**
   * @brief 从当前轨迹开始规划
   * @return 是否规划成功
   */
  bool planFromCurrentTraj();

  /**
   * @brief 改变状态机执行状态
   * @param new_state 新状态
   * @param pos_call 调用位置(用于调试)
   */
  void changeFSMExecState(FSM_EXEC_STATE new_state, string pos_call);

  /**
   * @brief 获取连续调用同一状态的次数
   * @return pair<连续调用次数, 当前状态>
   */
  std::pair<int, ReboReplanFSM::FSM_EXEC_STATE> timesOfConsecutiveStateCalls();

  /**
   * @brief 打印当前状态机状态
   */
  void printFSMExecState();

  /**
   * @brief 获取局部目标
   */
  void getLocalTarget();

  /* ROS functions */
  /**
   * @brief 状态机执行回调函数(定时器触发)
   * @param e 定时器事件
   */
  void execFSMCallback(const ros::TimerEvent& e);

  /**
   * @brief 碰撞检测回调函数(定时器触发)
   * @param e 定时器事件
   */
  void checkCollisionCallback(const ros::TimerEvent& e);

  /**
   * @brief 航点回调函数
   * @param msg 路径消息
   */
  void waypointCallback(const nav_msgs::PathConstPtr& msg);

  /**
   * @brief 里程计回调函数
   * @param msg 里程计消息
   */
  void odometryCallback(const nav_msgs::OdometryConstPtr& msg);

  /**
   * @brief 检查轨迹碰撞
   * @return 是否发生碰撞
   */
  bool checkCollision();

public:
  ReboReplanFSM(/* args */) {
  }
  ~ReboReplanFSM() {
  }

  /**
   * @brief 初始化状态机
   * @param nh ROS节点句柄
   */
  void init(ros::NodeHandle& nh);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

}  // namespace rebound_planner

#endif