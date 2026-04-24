#ifndef _UNIFORM_BSPLINE_H_
#define _UNIFORM_BSPLINE_H_

#include <Eigen/Eigen>
#include <algorithm>
#include <iostream>

using namespace std;

namespace rebound_planner {

/**
 * @brief 均匀B样条曲线类
 *
 * 实现非均匀B样条(也包括均匀B样条,作为非均匀的特例)。
 * 支持不同维度的B样条曲线表示和操作。
 */
class UniformBspline {
private:
  // B样条控制点矩阵
  // 每行表示一个控制点
  // 列数决定维度,如3D空间的N个点 -> Nx3矩阵
  Eigen::MatrixXd control_points_;

  // p_: B样条阶数
  // n_: 控制点数量-1 (n+1个控制点)
  // m_: 节点向量大小 (m = n+p+1)
  int             p_, n_, m_;
  // 节点向量
  Eigen::VectorXd u_;
  // 节点间隔 \delta t
  double          interval_;

  /**
   * @brief 获取导数的控制点
   * @return 导数控制点矩阵
   */
  Eigen::MatrixXd getDerivativeControlPoints();

  // 物理限制和时间调整比例
  double limit_vel_, limit_acc_, limit_ratio_, feasibility_tolerance_;

public:
  UniformBspline() {}
  
  /**
   * @brief 构造函数
   * @param points 控制点矩阵
   * @param order B样条阶数
   * @param interval 时间间隔
   */
  UniformBspline(const Eigen::MatrixXd& points, const int& order, const double& interval);
  
  ~UniformBspline();

  /**
   * @brief 获取控制点
   * @return 控制点矩阵
   */
  Eigen::MatrixXd get_control_points(void) { return control_points_; }

  /**
   * @brief 初始化为均匀B样条
   * @param points 控制点矩阵
   * @param order B样条阶数
   * @param interval 时间间隔
   */
  void setUniformBspline(const Eigen::MatrixXd& points, const int& order, const double& interval);

  // get / set basic bspline info

  /**
   * @brief 设置节点向量
   * @param knot 节点向量
   */
  void                                   setKnot(const Eigen::VectorXd& knot);

  /**
   * @brief 获取节点向量
   * @return 节点向量
   */
  Eigen::VectorXd                        getKnot();

  /**
   * @brief 获取控制点
   * @return 控制点矩阵
   */
  Eigen::MatrixXd                        getControlPoint();

  /**
   * @brief 获取时间间隔
   * @return 时间间隔
   */
  double                                 getInterval();

  /**
   * @brief 获取时间跨度
   * @param um 起始时间(输出)
   * @param um_p 结束时间(输出)
   * @return 是否成功
   */
  bool                                   getTimeSpan(double& um, double& um_p);

  // compute position / derivative

  /**
   * @brief 使用De Boor算法计算曲线上的点
   * @param u 参数值, u \in [up, u_mp]
   * @return 曲线上的点
   */
  Eigen::VectorXd   evaluateDeBoor(const double& u);

  /**
   * @brief 使用时间t计算曲线上的点
   * @param t 时间, t \in [0, duration]
   * @return 曲线上的点
   */
  inline Eigen::VectorXd   evaluateDeBoorT(const double& t) {return evaluateDeBoor(t + u_(p_));}

  /**
   * @brief 获取B样条的导数
   * @return 导数B样条曲线
   */
  UniformBspline getDerivative();

  /**
   * @brief 将点集参数化为B样条控制点
   *
   * 对3D点集进行B样条插值,同时满足边界速度和加速度约束。
   * @param ts 时间间隔
   * @param point_set 输入点集(K+2个点,含边界速度/加速度)
   * @param start_end_derivative 起点和终点的导数(速度、加速度)
   * @param ctrl_pts 输出的控制点(K+6个)
   */
  static void parameterizeToBspline(const double& ts, const vector<Eigen::Vector3d>& point_set,
                                    const vector<Eigen::Vector3d>& start_end_derivative,
                                    Eigen::MatrixXd&               ctrl_pts);

  /* check feasibility, adjust time */

  /**
   * @brief 设置物理限制
   * @param vel 最大速度
   * @param acc 最大加速度
   * @param tolerance 可行性容差
   */
  void   setPhysicalLimits(const double& vel, const double& acc, const double& tolerance);

  /**
   * @brief 检查可行性
   * @param ratio 超出限制的比例(输出)
   * @param show 是否显示详情
   * @return 是否可行
   */
  bool   checkFeasibility(double& ratio, bool show = false);

  /**
   * @brief 延长时间
   * @param ratio 延长比例
   */
  void   lengthenTime(const double& ratio);

  /* for performance evaluation */

  /**
   * @brief 获取轨迹总时间
   * @return 总时间
   */
  double getTimeSum();

  /**
   * @brief 获取轨迹长度
   * @param res 分辨率
   * @return 轨迹长度
   */
  double getLength(const double& res = 0.01);

  /**
   * @brief 获取加加速度(jerk)
   * @return 加加速度值
   */
  double getJerk();

  /**
   * @brief 获取平均和最大速度
   * @param mean_v 平均速度(输出)
   * @param max_v 最大速度(输出)
   */
  void   getMeanAndMaxVel(double& mean_v, double& max_v);

  /**
   * @brief 获取平均和最大加速度
   * @param mean_a 平均加速度(输出)
   * @param max_a 最大加速度(输出)
   */
  void   getMeanAndMaxAcc(double& mean_a, double& max_a);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
}  // namespace rebound_planner
#endif