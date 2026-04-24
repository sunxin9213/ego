#ifndef _PLANNING_VISUALIZATION_H_
#define _PLANNING_VISUALIZATION_H_

#include <Eigen/Eigen>
#include <algorithm>
#include <bspline_opt/uniform_bspline.h>
#include <iostream>
#include <poly_traj/polynomial_traj.h>
#include <ros/ros.h>
#include <vector>
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>
#include <stdlib.h>

using std::vector;

namespace rebound_planner {

/**
 * @brief 规划可视化类
 *
 * 负责在ROS RViz中可视化规划过程中的各种数据,包括:
 * - 目标点显示
 * - 全局路径显示
 * - 初始路径显示
 * - 最优轨迹显示
 * - A*搜索路径显示
 * - 引导向量显示
 */
class PlanningVisualization {
private:

  // ROS节点句柄
  ros::NodeHandle node;

  // 各类型数据的发布者
  ros::Publisher goal_point_pub;          // 目标点发布者
  ros::Publisher global_list_pub;         // 全局路径发布者
  ros::Publisher init_list_pub;           // 初始路径发布者
  ros::Publisher optimal_list_pub;       // 最优轨迹发布者
  ros::Publisher a_star_list_pub;         // A*路径发布者
  ros::Publisher guide_vector_pub;        // 引导向量发布者
  ros::Publisher intermediate_state_pub;   // 中间状态发布者

public:
  PlanningVisualization(/* args */) {}
  ~PlanningVisualization() {}
  
  /**
   * @brief 构造函数
   * @param nh ROS节点句柄
   */
  PlanningVisualization(ros::NodeHandle& nh);

  typedef std::shared_ptr<PlanningVisualization> Ptr;

  /**
   * @brief 显示标记点列表
   * @param pub 发布者
   * @param list 坐标点列表
   * @param scale 标记大小
   * @param color RGBA颜色
   * @param id 标记ID
   */
  void displayMarkerList(ros::Publisher& pub, const vector<Eigen::Vector3d>& list, double scale,
                                                Eigen::Vector4d color, int id);

  /**
   * @brief 生成路径显示数组
   * @param array 输出MarkerArray
   * @param list 坐标点列表
   * @param scale 线条粗细
   * @param color RGBA颜色
   * @param id 标记ID
   */
  void generatePathDisplayArray(visualization_msgs::MarkerArray &array,
    const vector<Eigen::Vector3d>& list, double scale, Eigen::Vector4d color, int id);

  /**
   * @brief 生成箭头显示数组
   * @param array 输出MarkerArray
   * @param list 坐标点列表
   * @param scale 箭头大小
   * @param color RGBA颜色
   * @param id 标记ID
   */
  void generateArrowDisplayArray(visualization_msgs::MarkerArray &array,
    const vector<Eigen::Vector3d>& list, double scale, Eigen::Vector4d color, int id);

  /**
   * @brief 显示目标点
   * @param goal_point 目标点位置
   * @param color RGBA颜色
   * @param scale 显示大小
   * @param id 标记ID
   */
  void displayGoalPoint(Eigen::Vector3d goal_point, Eigen::Vector4d color, const double scale, int id);

  /**
   * @brief 显示全局路径列表
   * @param global_pts 全局路径点
   * @param scale 线条粗细
   * @param id 标记ID
   */
  void displayGlobalPathList(vector<Eigen::Vector3d> global_pts, const double scale, int id);

  /**
   * @brief 显示初始路径列表
   * @param init_pts 初始路径点
   * @param scale 线条粗细
   * @param id 标记ID
   */
  void displayInitPathList(vector<Eigen::Vector3d> init_pts, const double scale, int id);

  /**
   * @brief 显示最优轨迹
   * @param optimal_pts 最优轨迹控制点
   * @param id 标记ID
   */
  void displayOptimalList(Eigen::MatrixXd optimal_pts, int id);

  /**
   * @brief 显示A*搜索路径列表
   * @param a_star_paths A*搜索得到的多条路径
   * @param id 标记ID
   */
  void displayAStarList(std::vector<std::vector<Eigen::Vector3d>> a_star_paths, int id);

  /**
   * @brief 显示箭头列表
   * @param pub 发布者
   * @param list 坐标点列表
   * @param scale 箭头大小
   * @param color RGBA颜色
   * @param id 标记ID
   */
  void displayArrowList(ros::Publisher& pub, const vector<Eigen::Vector3d>& list, double scale, Eigen::Vector4d color, int id);

};
}  // namespace rebound_planner
#endif
