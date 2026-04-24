#ifndef _PLAN_CONTAINER_H_
#define _PLAN_CONTAINER_H_

#include <Eigen/Eigen>
#include <vector>
#include <ros/ros.h>

#include <bspline_opt/uniform_bspline.h>
#include <poly_traj/polynomial_traj.h>

using std::vector;

namespace rebound_planner {

/**
 * @brief 全局轨迹数据结构
 *
 * 存储和管理全局与局部轨迹信息,支持轨迹的时间同步和状态查询。
 * 该类维护全局多项式轨迹和局部B样条轨迹,并提供统一的状态查询接口。
 */
class GlobalTrajData {
private:
public:
  // 全局多项式轨迹
  PolynomialTraj global_traj_;
  // 局部B样条轨迹集合(位置、速度、加速度)
  vector<UniformBspline> local_traj_;

  // 全局轨迹总时长
  double global_duration_;
  // 全局轨迹开始时间
  ros::Time global_start_time_;
  // 局部轨迹开始和结束时间
  double local_start_time_, local_end_time_;
  // 累计时间增量
  double time_increase_;
  // 上一次时间增量
  double last_time_inc_;

  GlobalTrajData(/* args */) {}

  ~GlobalTrajData() {}

  /**
   * @brief 检查局部轨迹是否到达目标
   * @return 局部轨迹是否到达终点
   */
  bool localTrajReachTarget() { return fabs(local_end_time_ - global_duration_) < 0.1; }

  /**
   * @brief 设置全局轨迹
   * @param traj 全局多项式轨迹
   * @param time 全局轨迹开始时间
   */
  void setGlobalTraj(const PolynomialTraj& traj, const ros::Time& time) {
    global_traj_ = traj;
    global_traj_.init();
    global_duration_ = global_traj_.getTimeSum();
    global_start_time_ = time;

    local_traj_.clear();
    local_start_time_ = -1;
    local_end_time_ = -1;
    time_increase_ = 0.0;
    last_time_inc_ = 0.0;
  }

  /**
   * @brief 设置局部轨迹
   * @param traj 局部B样条轨迹
   * @param local_ts 局部轨迹开始时间
   * @param local_te 局部轨迹结束时间
   * @param time_inc 时间增量
   */
  void setLocalTraj(UniformBspline traj, double local_ts, double local_te, double time_inc) {
    local_traj_.resize(3);
    local_traj_[0] = traj;
    local_traj_[1] = local_traj_[0].getDerivative();
    local_traj_[2] = local_traj_[1].getDerivative();

    local_start_time_ = local_ts;
    local_end_time_ = local_end_time_;
    global_duration_ += time_inc;
    time_increase_ += time_inc;
    last_time_inc_ = time_inc;
  }

  /**
   * @brief 获取指定时间的位置
   * @param t 时间(相对于全局轨迹开始时间)
   * @return 位置向量
   */
  Eigen::Vector3d getPosition(double t) {
    if (t >= -1e-3 && t <= local_start_time_) {
      return global_traj_.evaluate(t - time_increase_ + last_time_inc_);
    } else if (t >= local_end_time_ && t <= global_duration_ + 1e-3) {
      return global_traj_.evaluate(t - time_increase_);
    } else {
      double tm, tmp;
      local_traj_[0].getTimeSpan(tm, tmp);
      return local_traj_[0].evaluateDeBoor(tm + t - local_start_time_);
    }
  }

  /**
   * @brief 获取指定时间的速度
   * @param t 时间(相对于全局轨迹开始时间)
   * @return 速度向量
   */
  Eigen::Vector3d getVelocity(double t) {
    if (t >= -1e-3 && t <= local_start_time_) {
      return global_traj_.evaluateVel(t);
    } else if (t >= local_end_time_ && t <= global_duration_ + 1e-3) {
      return global_traj_.evaluateVel(t - time_increase_);
    } else {
      double tm, tmp;
      local_traj_[0].getTimeSpan(tm, tmp);
      return local_traj_[1].evaluateDeBoor(tm + t - local_start_time_);
    }
  }

  /**
   * @brief 获取指定时间的加速度
   * @param t 时间(相对于全局轨迹开始时间)
   * @return 加速度向量
   */
  Eigen::Vector3d getAcceleration(double t) {
    if (t >= -1e-3 && t <= local_start_time_) {
      return global_traj_.evaluateAcc(t);
    } else if (t >= local_end_time_ && t <= global_duration_ + 1e-3) {
      return global_traj_.evaluateAcc(t - time_increase_);
    } else {
      double tm, tmp;
      local_traj_[0].getTimeSpan(tm, tmp);
      return local_traj_[2].evaluateDeBoor(tm + t - local_start_time_);
    }
  }

  // get Bspline paramterization data of a local trajectory within a sphere
  // start_t: start time of the trajectory
  // dist_pt: distance between the discretized points
  void getTrajByRadius(const double& start_t, const double& des_radius, const double& dist_pt,
                       vector<Eigen::Vector3d>& point_set, vector<Eigen::Vector3d>& start_end_derivative,
                       double& dt, double& seg_duration) {
    double seg_length = 0.0;  // length of the truncated segment
    double seg_time = 0.0;    // duration of the truncated segment
    double radius = 0.0;      // distance to the first point of the segment

    double delt = 0.2;
    Eigen::Vector3d first_pt = getPosition(start_t);  // first point of the segment
    Eigen::Vector3d prev_pt = first_pt;               // previous point
    Eigen::Vector3d cur_pt;                           // current point

    // go forward until the traj exceed radius or global time

    while (radius < des_radius && seg_time < global_duration_ - start_t - 1e-3) {
      seg_time += delt;
      seg_time = min(seg_time, global_duration_ - start_t);

      cur_pt = getPosition(start_t + seg_time);
      seg_length += (cur_pt - prev_pt).norm();
      prev_pt = cur_pt;
      radius = (cur_pt - first_pt).norm();
    }

    // get parameterization dt by desired density of points
    int seg_num = floor(seg_length / dist_pt);

    // get outputs

    seg_duration = seg_time;  // duration of the truncated segment
    dt = seg_time / seg_num;  // time difference between two points

    for (double tp = 0.0; tp <= seg_time + 1e-4; tp += dt) {
      cur_pt = getPosition(start_t + tp);
      point_set.push_back(cur_pt);
    }

    start_end_derivative.push_back(getVelocity(start_t));
    start_end_derivative.push_back(getVelocity(start_t + seg_time));
    start_end_derivative.push_back(getAcceleration(start_t));
    start_end_derivative.push_back(getAcceleration(start_t + seg_time));
  }

  // get Bspline paramterization data of a fixed duration local trajectory
  // start_t: start time of the trajectory
  // duration: time length of the segment
  // seg_num: discretized the segment into *seg_num* parts
  void getTrajByDuration(double start_t, double duration, int seg_num,
                         vector<Eigen::Vector3d>& point_set,
                         vector<Eigen::Vector3d>& start_end_derivative, double& dt) {
    dt = duration / seg_num;
    Eigen::Vector3d cur_pt;
    for (double tp = 0.0; tp <= duration + 1e-4; tp += dt) {
      cur_pt = getPosition(start_t + tp);
      point_set.push_back(cur_pt);
    }

    start_end_derivative.push_back(getVelocity(start_t));
    start_end_derivative.push_back(getVelocity(start_t + duration));
    start_end_derivative.push_back(getAcceleration(start_t));
    start_end_derivative.push_back(getAcceleration(start_t + duration));
  }
};

/**
 * @brief 规划参数结构体
 *
 * 存储运动规划算法所需的各种物理限制和处理时间统计。
 */
struct PlanParameters {
  /* planning algorithm parameters */
  // 物理限制: 最大速度、最大加速度、最大加加速度
  double max_vel_, max_acc_, max_jerk_;
  // B样条控制点间距
  double ctrl_pt_dist;
  // 允许超出限制的速度/加速度比例
  double feasibility_tolerance_;
  // 规划视界距离
  double planning_horizen_;

  /* processing time */
  // 各阶段处理时间
  double time_search_ = 0.0;    // 搜索时间
  double time_optimize_ = 0.0; // 优化时间
  double time_adjust_ = 0.0;   // 调整时间
};

/**
 * @brief 局部轨迹数据结构
 *
 * 存储生成的局部轨迹信息,包括轨迹ID、时间、起点状态和轨迹曲线。
 */
struct LocalTrajData {
  /* info of generated traj */

  // 轨迹ID
  int traj_id_;
  // 轨迹持续时间
  double duration_;
  // 全局时间偏移量(当局部轨迹完成切换回全局轨迹时,全局轨迹时间与世界时间不再匹配)
  double global_time_offset;
  // 轨迹开始时间
  ros::Time start_time_;
  // 轨迹起始位置
  Eigen::Vector3d start_pos_;
  // 位置、速度、加速度B样条轨迹
  UniformBspline position_traj_, velocity_traj_, acceleration_traj_;
};

}  // namespace rebound_planner

#endif