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
 * 【B样条曲线简介】
 * B样条(B-spline)是贝塞尔曲线(Bezier)的一般化形式,具有以下优点:
 * - 局部支撑性: 移动一个控制点只影响曲线的局部形状
 * - 连续性: 可以通过控制点间接控制曲线的光滑性(C^n连续)
 * - 仿射不变性: 仿射变换后仍是同一组控制点的B样条
 *
 * 【B样条曲线方程】
 * C(u) = Σ_{i=0}^{n} N_{i,p}(u) * P_i
 * 其中:
 * - P_i 是第i个控制点
 * - N_{i,p}(u) 是p阶B样条基函数
 * - u 是参数值,定义在节点向量上
 * - n+1 是控制点数量, p 是曲线阶数
 *
 * 本实现使用均匀B样条(节点向量均匀分布),但也兼容非均匀情况。
 * 支持3D空间中的轨迹表示,每个控制点是一个三维向量。
 */
class UniformBspline {
private:
  // B样条控制点矩阵: 每列表示一个控制点,行数表示维度(3D)
  // 例如: N个3D控制点 -> 3xN 矩阵
  Eigen::MatrixXd control_points_;

  // p_: B样条阶数(阶数=次数+1,如3阶对应2次曲线)
  //     【注】3阶B样条的加加速度(jerk)是连续的→适合无人机轨迹
  int             p_;
  // n_: 控制点数量-1 (即共有 n+1 个控制点)
  int             n_;
  // m_: 节点向量大小 (m = n + p + 1)
  int             m_;
  // 节点向量(knot vector): u_[0], u_[1], ..., u_[m]
  // 均匀B样条的节点向量是均匀递增的: u_[i+1] - u_[i] = interval
  // 前p+1个节点值相同(在起点重复),后p+1个节点值相同(在终点重复)
  // 这是Clamped B-spline(夹紧B样条)的特征,确保曲线经过首尾控制点
  Eigen::VectorXd u_;
  // 节点间隔(均匀B样条的节点间距)
  double          interval_;

  /**
   * @brief 获取导数曲线的控制点
   * 
   * B样条的一个重要性质: p阶B样条的导数是一个(p-1)阶B样条。
   * 导数控制点 Q_i 由原控制点计算得到:
   *   Q_i = p * (P_{i+1} - P_i) / (u_{i+p+1} - u_{i+1})
   * 导数曲线有 n 个控制点(比原曲线少1个)
   * 
   * @return 导数控制点矩阵
   */
  Eigen::MatrixXd getDerivativeControlPoints();

  // 物理限制参数
  double limit_vel_, limit_acc_, limit_ratio_, feasibility_tolerance_;

public:
  UniformBspline() {}
  
  /**
   * @brief 构造函数
   * @param points 控制点矩阵(维度x数量)
   * @param order B样条阶数(如3阶→三次B样条)
   * @param interval 节点时间间隔
   */
  UniformBspline(const Eigen::MatrixXd& points, const int& order, const double& interval);
  
  ~UniformBspline();

  /**
   * @brief 获取控制点矩阵
   * @return 控制点矩阵
   */
  Eigen::MatrixXd get_control_points(void) { return control_points_; }

  /**
   * @brief 初始化为均匀B样条
   * 
   * 构造一个clamped均匀B样条,节点向量分布如下:
   * - 前p+1个节点: -p*interval 到 0 (步长interval)
   * - 中间节点: 0 到 (n-p)*interval (步长interval)
   * - 后p+1个节点: 从 (n-p)*interval 继续增加到 n*interval (步长interval)
   * 这样构造的节点向量使得曲线定义域为 [0, (n-p)*interval]
   * 
   * @param points 控制点矩阵
   * @param order B样条阶数
   * @param interval 时间间隔
   */
  void setUniformBspline(const Eigen::MatrixXd& points, const int& order, const double& interval);

  /**
   * @brief 设置节点向量(用于非均匀B样条)
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
   * @return 时间间隔值
   */
  double                                 getInterval();

  /**
   * @brief 获取曲线定义域的时间范围
   * 
   * clamped B样条的有效定义域是 [u_(p), u_(m-p)]
   * 即从第p+1个节点到第m-p+1个节点
   * 
   * @param um 起始时间(输出)
   * @param um_p 结束时间(输出)
   * @return 是否成功获取
   */
  bool                                   getTimeSpan(double& um, double& um_p);

  /**
   * @brief 使用De Boor算法计算曲线上参数u对应的点
   * 
   * 【De Boor算法原理】
   * De Boor算法是B样条求值的标准算法,类似于贝塞尔曲线的De Casteljau算法。
   * 
   * 算法步骤:
   * 1. 确定u所在的节点区间 [u_k, u_{k+1})
   * 2. 取与该区间相关的p+1个控制点: P_{k-p}, ..., P_k
   * 3. 对r=1,2,...,p层递归,每层在两两控制点之间做线性插值:
   *    d_i^{(r)} = (1-α_i^{(r)}) * d_{i-1}^{(r-1)} + α_i^{(r)} * d_i^{(r-1)}
   *    其中 α_i^{(r)} = (u - u_{i+k-p}) / (u_{i+1+k-r} - u_{i+k-p})
   * 4. 最终 d_p^{(p)} 即为曲线上的点
   * 
   * @param u 参数值, u ∈ [u_(p), u_(m-p)]
   * @return 曲线上的点(维度与控制点相同)
   */
  Eigen::VectorXd   evaluateDeBoor(const double& u);

  /**
   * @brief 使用时间t(从0开始)计算曲线上的点
   * 
   * 时间t从0开始,通过 t + u_(p) 映射到参数u的定义域
   * 
   * @param t 时间, t ∈ [0, duration]
   * @return 曲线上的点
   */
  inline Eigen::VectorXd   evaluateDeBoorT(const double& t) {return evaluateDeBoor(t + u_(p_));}

  /**
   * @brief 获取B样条的导数B样条
   * 
   * B样条求导: p阶B样条求导后得到(p-1)阶B样条。
   * 可以对轨迹求多次导获得速度、加速度、加加速度曲线。
   * 例如: position.getDerivative() → velocity
   *       velocity.getDerivative() → acceleration
   * 
   * @return 导数B样条曲线
   */
  UniformBspline getDerivative();

  /**
   * @brief 将离散点集参数化为B样条控制点
   * 
   * 【核心功能:离散点→连续B样条轨迹】
   * 
   * 输入一系列离散路径点(如A*搜索得到的路径),通过求解线性方程组
   * 得到对应的B样条控制点,从而将离散路径转化为光滑的B样条曲线。
   * 
   * 具体方法:
   * 1. 每个B样条曲线段上的点由三个相邻控制点加权平均得到,
   *    权重为(1,4,1)/6(三次均匀B样条在节点中点的求值公式)
   * 2. 同时附加起点和终点的速度、加速度约束(可通过边界导数方程表达)
   * 3. 构造线性系统 Ax = b 求解,其中A是约束矩阵,x是待求控制点,b是已知点位置/导数
   * 
   * 约束方程的维度:
   * - K个离散路径点 → K个位置约束方程
   * - 4个边界导数约束(起点速度、终点速度、起点加速度、终点加速度)
   * - 共K+4个方程,求解K+2个控制点(超定系统,最小二乘求解)
   * 
   * @param ts 时间步长(节点间隔)
   * @param point_set 离散路径点序列(共K个点)
   * @param start_end_derivative 边界导数(4个向量: v0, vf, a0, af)
   * @param ctrl_pts 输出的控制点矩阵(K+2个控制点)
   */
  static void parameterizeToBspline(const double& ts, const vector<Eigen::Vector3d>& point_set,
                                    const vector<Eigen::Vector3d>& start_end_derivative,
                                    Eigen::MatrixXd&               ctrl_pts);

  /**
   * @brief 设置物理限制
   * 
   * 设置无人机(或机器人)的速度、加速度限制,用于后续的可行性检查
   * 
   * @param vel 最大速度
   * @param acc 最大加速度
   * @param tolerance 可行性容差(允许超出限制的比例)
   */
  void   setPhysicalLimits(const double& vel, const double& acc, const double& tolerance);

  /**
   * @brief 检查轨迹的可行性(速度和加速度是否在限制范围内)
   * 
   * 利用B样条凸包性质: B样条曲线上的速度和加速度一定在
   * 其对应控制点构成的速度/加速度控制多边形的凸包内。
   * 
   * 因此只需检查速度控制点和加速度控制点是否在限制范围内即可,
   * 无需密集采样整条轨迹。
   * 
   * 如果不可行,输出超出最大限制的比例ratio。
   * 
   * @param ratio 超出最大限制的比例(输出)
   * @param show 是否打印详细信息
   * @return true 可行; false 不可行
   */
  bool   checkFeasibility(double& ratio, bool show = false);

  /**
   * @brief 延长轨迹时间(降低速度/加速度以满足可行性)
   * 
   * 当速度和加速度超出限制时,通过延长轨迹的总时间来降低
   * 速度/加速度幅值。延长时间等价于拉伸节点向量,同时保持
   * 控制点位置不变,从而降低控制点间的变化率(即速度/加速度)。
   * 
   * @param ratio 时间延长比例
   */
  void   lengthenTime(const double& ratio);

  /* 下面对轨迹进行性能评估 */

  /**
   * @brief 获取轨迹的总时间(定义域长度)
   * @return 总时间值
   */
  double getTimeSum();

  /**
   * @brief 获取轨迹的弧长
   * @param res 采样分辨率
   * @return 轨迹长度
   */
  double getLength(const double& res = 0.01);

  /**
   * @brief 获取轨迹的加加速度积分值(平滑度指标)
   * 
   * 加加速度(jerk)是位置的三阶导数。
   * 加加速度的平方积分 ∫(jerk)²dt 是衡量轨迹平滑度的重要指标,
   * 最小化该值可以得到光滑的轨迹。
   * 
   * @return 加加速度积分值
   */
  double getJerk();

  /**
   * @brief 获取平均速度和最大速度
   * @param mean_v 平均速度(输出)
   * @param max_v 最大速度(输出)
   */
  void   getMeanAndMaxVel(double& mean_v, double& max_v);

  /**
   * @brief 获取平均加速度和最大加速度
   * @param mean_a 平均加速度(输出)
   * @param max_a 最大加速度(输出)
   */
  void   getMeanAndMaxAcc(double& mean_a, double& max_a);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
}  // namespace rebound_planner
#endif