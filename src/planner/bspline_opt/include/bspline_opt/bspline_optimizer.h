#ifndef _BSPLINE_OPTIMIZER_H_
#define _BSPLINE_OPTIMIZER_H_

#include <Eigen/Eigen>
#include <path_searching/dyn_a_star.h>
#include <bspline_opt/uniform_bspline.h>
#include <plan_env/grid_map.h>
#include <ros/ros.h>

/**
 * @brief B样条轨迹优化器
 *
 * 基于梯度和弹性带(Elastic Band)优化的B样条轨迹优化器。
 * 输入: 符号距离场(SDF)和一系列控制点
 * 输出: 优化后的控制点序列
 * 点格式: N x 3 矩阵,每行是一个三维点
 */
namespace rebound_planner {

/**
 * @brief 控制点类
 *
 * 存储B样条轨迹优化所需的控制点数据,包括点集、基点和方向向量等信息。
 */
class ControlPoints
{
public:

  // 间隙/清除距离(用于碰撞避免)
  double clearance;
  // 控制点数量
  int size;
  // 控制点矩阵(3行,size列)
  Eigen::MatrixXd points;
  
  // 基点向量集合(方向向量的起点,即碰撞点)
  std::vector<std::vector<Eigen::Vector3d>> base_point;
  // 方向向量(必须归一化)
  std::vector<std::vector<Eigen::Vector3d>> direction;
  // 临时标志位
  std::vector<bool> flag_temp;

  /**
   * @brief 调整控制点集合大小
   * @param size_set 新的控制点数量
   */
  void resize( const int size_set )
  {
    size = size_set;

    base_point.clear();
    direction.clear();
    flag_temp.clear();

    points.resize(3, size_set);
    base_point.resize(size);
    direction.resize(size);
    flag_temp.resize(size);
  }
};

/**
 * @brief B样条优化器类
 *
 * 实现B样条轨迹的优化功能,结合梯度下降法和弹性带方法:
 * - 平滑度优化: 最小化加加速度(jerk)
 * - 距离优化: 保持与障碍物的安全距离
 * - 可行性优化: 满足速度和加速度约束
 * - 曲线拟合: 追踪引导路径
 */
class BsplineOptimizer {

public:

  BsplineOptimizer() {}
  ~BsplineOptimizer() {}

  /* main API */

  /**
   * @brief 设置环境地图
   * @param env 栅格地图指针
   */
  void            setEnvironment(const GridMap::Ptr& env);

  /**
   * @brief 设置优化参数
   * @param nh ROS节点句柄
   */
  void            setParam(ros::NodeHandle& nh);

  /**
   * @brief B样条轨迹优化(主接口)
   * @param points 初始控制点
   * @param ts 时间间隔
   * @param cost_function 代价函数类型
   * @param max_num_id 最大迭代次数ID
   * @param max_time_id 最大时间限制ID
   * @return 优化后的控制点
   */
  Eigen::MatrixXd BsplineOptimizeTraj(const Eigen::MatrixXd& points, const double& ts,
                                      const int& cost_function, int max_num_id, int max_time_id);

  /* helper function */

  // required inputs

  /**
   * @brief 设置控制点
   * @param points 控制点矩阵
   */
  void setControlPoints(const Eigen::MatrixXd& points);

  /**
   * @brief 设置B样条时间间隔
   * @param ts 时间间隔
   */
  void setBsplineInterval(const double& ts);

  /**
   * @brief 设置代价函数
   * @param cost_function 代价函数类型
   */
  void setCostFunction(const int& cost_function);

  /**
   * @brief 设置终止条件
   * @param max_num_id 最大迭代次数ID
   * @param max_time_id 最大时间限制ID
   */
  void setTerminateCond(const int& max_num_id, const int& max_time_id);

  // optional inputs

  /**
   * @brief 设置引导路径
   * @param guide_pt 引导点序列
   */
  void setGuidePath(const vector<Eigen::Vector3d>& guide_pt);

  /**
   * @brief 设置航点约束
   * @param waypts 航点位置
   * @param waypt_idx 航点索引(最多N-2个约束)
   */
  void setWaypoints(const vector<Eigen::Vector3d>& waypts,
                    const vector<int>&             waypt_idx);

  /**
   * @brief 执行优化
   */
  void optimize();

  /**
   * @brief 获取优化后的控制点
   * @return 控制点矩阵
   */
  Eigen::MatrixXd         getControlPoints();

  // A*搜索器指针
  AStar::Ptr a_star_;
  // 参考点序列
  std::vector<Eigen::Vector3d> ref_pts_;

  /**
   * @brief 初始化控制点
   * @param init_points 初始点
   * @param flag_first_init 是否首次初始化
   * @return 初始化后的控制点集合
   */
  std::vector<std::vector<Eigen::Vector3d>> initControlPoints(Eigen::MatrixXd& init_points, bool flag_first_init = true);

  /**
   * @brief 反弹优化轨迹(必须先调用initControlPoints)
   * @param optimal_points 最优控制点输出
   * @param ts 时间间隔
   * @return 是否优化成功
   */
  bool BsplineOptimizeTrajRebound(Eigen::MatrixXd& optimal_points, double ts);

  /**
   * @brief 精修优化轨迹
   * @param init_points 初始控制点
   * @param ts 时间间隔
   * @param optimal_points 最优控制点输出
   * @return 是否优化成功
   */
  bool BsplineOptimizeTrajRefine(const Eigen::MatrixXd& init_points, const double ts, Eigen::MatrixXd& optimal_points);

  /**
   * @brief 获取B样条阶数
   * @return 阶数
   */
  inline int getOrder(void) { return order_; }

private:
  // 栅格地图指针
  GridMap::Ptr grid_map_;

  /**
   * @brief 强制停止优化类型枚举
   */
  enum FORCE_STOP_OPTIMIZE_TYPE
  {
    DONT_STOP,          // 不停止
    STOP_FOR_REBOUND,   // 为反弹停止
    STOP_FOR_ERROR      // 为错误停止
  } force_stop_type_;

  // main input
  // B样条节点间隔
  double          bspline_interval_;
  // 轨迹终点
  Eigen::Vector3d end_pt_;
                                        //
  // 几何引导路径点(N-6个)
  vector<Eigen::Vector3d> guide_pts_;
  // 航点约束
  vector<Eigen::Vector3d> waypoints_;
  // 航点约束索引
  vector<int>             waypt_idx_;
                                        //
  // 停止标准:最大迭代次数和最大时间
  int    max_num_id_, max_time_id_;
  // 代价函数类型
  int    cost_function_;
  // 移动障碍物的全局时间
  double start_time_;

  /* optimization parameters */
  // B样条阶数
  int    order_;
  // 权重参数
  double lambda1_;                // 加加速度平滑权重
  double lambda2_;                // 距离权重
  double lambda3_;                // 可行性权重
  double lambda4_;                // 曲线拟合权重

  int a;
                                   //
  // 安全距离
  double dist0_;
  // 动态限制:最大速度和加速度
  double max_vel_, max_acc_;

  // 优化变量数量
  int                 variable_num_;
  // 求解器迭代次数
  int                 iter_num_;
  // 最佳变量值
  Eigen::VectorXd     best_variable_;
  // 最小代价
  double              min_cost_;

  // 控制点集合
  ControlPoints cps_;

  /* cost function */
  /* calculate each part of cost function with control points q as input */

  /**
   * @brief 代价函数(静态回调函数)
   * @param x 优化变量
   * @param grad 梯度输出
   * @param func_data 函数数据
   * @return 代价值
   */
  static double costFunction(const std::vector<double>& x, std::vector<double>& grad, void* func_data);

  /**
   * @brief 组合代价函数
   * @param x 优化变量
   * @param grad 梯度输出
   * @param cost 代价输出
   */
  void          combineCost(const std::vector<double>& x, vector<double>& grad, double& cost);

  // q包含所有控制点

  /**
   * @brief 计算平滑度代价
   * @param q 控制点
   * @param cost 代价输出
   * @param gradient 梯度输出
   * @param falg_use_jerk 是否使用加加速度
   */
  void calcSmoothnessCost(const Eigen::MatrixXd& q, double& cost,
                          Eigen::MatrixXd& gradient, bool falg_use_jerk = true);

  /**
   * @brief 计算可行性代价
   * @param q 控制点
   * @param cost 代价输出
   * @param gradient 梯度输出
   */
  void calcFeasibilityCost(const Eigen::MatrixXd& q, double& cost,
                           Eigen::MatrixXd& gradient);

  /**
   * @brief 计算反弹距离代价
   * @param q 控制点
   * @param cost 代价输出
   * @param gradient 梯度输出
   * @param iter_num 迭代次数
   * @param smoothness_cost 平滑度代价
   */
  void calcDistanceCostRebound(const Eigen::MatrixXd& q, double& cost, Eigen::MatrixXd& gradient, int iter_num, double smoothness_cost);

  /**
   * @brief 计算拟合代价
   * @param q 控制点
   * @param cost 代价输出
   * @param gradient 梯度输出
   */
  void calcFitnessCost(const Eigen::MatrixXd& q, double& cost, Eigen::MatrixXd& gradient);

  /**
   * @brief 检查碰撞并进行反弹
   * @return 是否需要反弹
   */
  bool check_collision_and_rebound(void);

  /**
   * @brief 反弹优化代价函数
   */
  static double costFunctionRebound(const Eigen::VectorXd& x, Eigen::VectorXd& grad, bool& force_return, void* func_data);

  /**
   * @brief 精修优化代价函数
   */
  static double costFunctionRefine(const Eigen::VectorXd& x, Eigen::VectorXd& grad, bool& force_return, void* func_data);

  /**
   * @brief 反弹优化
   * @return 是否优化成功
   */
  bool rebound_optimize();

  /**
   * @brief 精修优化
   * @return 是否优化成功
   */
  bool refine_optimize();

  /**
   * @brief 组合反弹代价
   * @param x 优化变量
   * @param grad 梯度输出
   * @param f_combine 组合代价输出
   */
  void combineCostRebound(const Eigen::VectorXd& x, Eigen::VectorXd& grad, double& f_combine);

  /**
   * @brief 组合精修代价
   * @param x 优化变量
   * @param grad 梯度输出
   * @param f_combine 组合代价输出
   */
  void combineCostRefine(const Eigen::VectorXd& x, Eigen::VectorXd& grad, double& f_combine);

  /* for benckmark evaluation only */
public:

  typedef unique_ptr<BsplineOptimizer> Ptr;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

}  // namespace rebound_planner
#endif