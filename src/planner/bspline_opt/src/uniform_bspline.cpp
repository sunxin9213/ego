#include "bspline_opt/uniform_bspline.h"
#include <ros/ros.h>

namespace rebound_planner {

/**
 * @brief 构造函数: 用控制点、阶数和时间间隔初始化均匀B样条
 * 
 * @param points 控制点矩阵 (维度 × 控制点数量, 如3×N)
 * @param order B样条阶数(p值), 3阶=三次B样条
 * @param interval 节点向量中相邻节点的时间间隔
 */
UniformBspline::UniformBspline(const Eigen::MatrixXd& points, const int& order,
                                     const double& interval) {
  setUniformBspline(points, order, interval);
}

UniformBspline::~UniformBspline() {}

/**
 * @brief 初始化/重置为均匀B样条
 * 
 * 【节点向量构造细节】
 * 均匀clamped B样条的节点向量分为三部分:
 * 
 * 1. 起始重复段: u[0]~u[p] = (-p*interval, ..., 0)
 *    前p+1个节点从-p*interval以步长interval递增到0
 *    这部分使曲线从第一个控制点开始(clamped性质)
 *    
 * 2. 中间均匀段: u[p+1]~u[m-p] 以步长interval均匀递增
 *    这是曲线的实际定义域 [u[p], u[m-p]]
 *    对应的时间范围是 [0, (n-p)*interval]
 *    
 * 3. 结束重复段: u[m-p+1]~u[m] 继续以步长interval递增
 *    后p+1个节点使曲线在最后一个控制点结束(clamped性质)
 * 
 * 【示例】n=5个控制点(控制点0~4), p=3阶, interval=0.1
 *   m = n + p + 1 = 9
 *   节点向量: [-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
 *                           ↑定义域[0, 0.3]         ↑
 *   曲线定义域: u[3]~u[6] = [0, 0.3], 即时间范围0~0.3秒
 * 
 * @param points 控制点矩阵
 * @param order B样条阶数
 * @param interval 节点时间间隔
 */
void UniformBspline::setUniformBspline(const Eigen::MatrixXd& points, const int& order,
                                          const double& interval) {
  control_points_ = points;
  p_              = order;
  interval_       = interval;

  n_ = points.cols() - 1;  // 控制点数量-1
  m_ = n_ + p_ + 1;         // 节点向量大小

  u_ = Eigen::VectorXd::Zero(m_ + 1);
  for (int i = 0; i <= m_; ++i) {

    if (i <= p_) {
      // 起始段: 从 -p*interval 到 0
      u_(i) = double(-p_ + i) * interval_;
    } else if (i > p_ && i <= m_ - p_) {
      // 中间均匀段: 从 0 到 (n-p)*interval
      u_(i) = u_(i - 1) + interval_;
    } else if (i > m_ - p_) {
      // 结束段: 继续均匀递增  
      u_(i) = u_(i - 1) + interval_;
    }
  }
}

void UniformBspline::setKnot(const Eigen::VectorXd& knot) { this->u_ = knot; }

Eigen::VectorXd UniformBspline::getKnot() { return this->u_; }

/**
 * @brief 获取曲线定义域的时间范围
 * 
 * 对于clamped B样条,有效定义域是 [u_(p_), u_(m_-p_)]
 * 这个范围的物理含义:
 * - u_(p_) 通常等于0(起点时间)
 * - u_(m_-p_) = (n_-p_)*interval_(总持续时间)
 * 
 * @param um 输出: 起始时间(通常为0)
 * @param um_p 输出: 结束时间(轨迹总时长)
 * @return true 成功获取; false 索引越界
 */
bool UniformBspline::getTimeSpan(double& um, double& um_p) {
  if ( p_ > u_.rows() || m_ - p_ > u_.rows() ) return false;

  um   = u_(p_);
  um_p = u_(m_ - p_);
  
  return true;
}

Eigen::MatrixXd UniformBspline::getControlPoint() { return control_points_; }


/**
 * @brief De Boor算法: B样条曲线求值的核心算法
 * 
 * 【算法背景】
 * B样条曲线 C(u) 由控制点 P_i 和基函数 N_{i,p}(u) 定义:
 *   C(u) = Σ_{i=0}^{n} N_{i,p}(u) * P_i,  u ∈ [u_p, u_{m-p}]
 * 
 * De Boor算法是B样条求值的标准数值方法,其核心思想是:
 * 将B样条的求值转化为一系列线性插值,逐步消去中间控制点,
 * 最终得到曲线上的点。
 * 
 * 【算法步骤详解】
 * 设参数 u 位于区间 [u_k, u_{k+1}) 内:
 * 
 * 第0层(d[0]~d[p]): 从控制点矩阵中提取p+1个相关控制点
 *   d_i^(0) = P_{k-p+i},  i=0,1,...,p
 * 
 * 第r层(r=1,2,...,p):
 *   在相邻的d之间做线性插值,每层减少一个点:
 *   d_i^(r) = (1-α_i^(r)) * d_{i-1}^(r-1) + α_i^(r) * d_i^(r-1)
 *   其中插值系数:
 *   α_i^(r) = (u - u_{i+k-p}) / (u_{i+1+k-r} - u_{i+k-p})
 *         ↑当前参数到左节点距离 / 两节点间距
 * 
 * 最终结果:
 *   C(u) = d_p^(p), 即最后一层剩下的唯一一个点
 * 
 * 【可视化理解】
 * 对一个3阶B样条(p=3):
 *   第0层: [P₀] [P₁] [P₂] [P₃]  ← 4个相关控制点
 *          ↘↙   ↘↙   ↘↙
 *   第1层:  [Q₀] [Q₁] [Q₂]       ← 3个线性插值点
 *           ↘↙   ↘↙
 *   第2层:   [R₀] [R₁]            ← 2个线性插值点
 *            ↘↙
 *   第3层:    [S₀]                ← 曲线上的点
 * 
 * @param u 参数值,范围在 [u_(p_), u_(m_-p_)] 之间
 * @return 曲线在参数u处的点(三维向量)
 */
Eigen::VectorXd UniformBspline::evaluateDeBoor(const double& u) {

  // 将u钳制在有效定义域内
  double ub = min(max(u_(p_), u), u_(m_ - p_));

  // ---- 第1步: 确定参数u所在的节点区间 ----
  // 找到k使得 u_k ≤ u < u_{k+1}
  // 从k=p_开始搜索(因为定义域从u_(p_)开始)
  int k = p_;
  while (true) {
    if (u_(k + 1) >= ub) break;
    ++k;
  }

  // ---- 第2步: 提取相关控制点(第0层) ----
  // 在区间[u_k, u_{k+1})上,曲线由p+1个控制点决定:
  // P_{k-p}, P_{k-p+1}, ..., P_k
  vector<Eigen::VectorXd> d;
  for (int i = 0; i <= p_; ++i) {
    d.push_back(control_points_.col(k - p_ + i));
  }

  // ---- 第3步: 逐层插值(第1层到第p层) ----
  // 外层循环r控制层数,内层循环i控制当前层的计算
  for (int r = 1; r <= p_; ++r) {
    for (int i = p_; i >= r; --i) {
      // 计算插值系数alpha:
      // 分子: u - u_{i+k-p}  (当前参数到左端点的距离)
      // 分母: u_{i+1+k-r} - u_{i+k-p} (节点跨度)
      double alpha = (ub - u_[i + k - p_]) / (u_[i + 1 + k - r] - u_[i + k - p_]);
      // 线性插值: 新旧点的凸组合
      d[i] = (1 - alpha) * d[i - 1] + alpha * d[i];
    }
  }

  // 返回经过p层插值后的结果
  return d[p_];
}

// Eigen::VectorXd UniformBspline::evaluateDeBoorT(const double& t) {
//   return evaluateDeBoor(t + u_(p_));
// }

/**
 * @brief 计算B样条导数曲线的控制点
 * 
 * 【B样条求导公式】
 * B样条的一个重要性质: p阶B样条的导数是一个(p-1)阶B样条。
 * 其控制点 Q_i 由原控制点通过差分运算得到:
 * 
 *   Q_i = p * (P_{i+1} - P_i) / (u_{i+p+1} - u_{i+1})
 * 
 * 其中:
 * - P_i 是原B样条的控制点
 * - P_{i+1} - P_i 是控制点间的差分(速度方向)
 * - u_{i+p+1} - u_{i+1} 是对应的节点跨度(时间差)
 * - p 的乘法因子来自B样条基函数求导
 * 
 * 【物理意义】
 * 对于无人机轨迹(position),导数控制点具有明确的物理含义:
 * - 位置B样条 → 速度控制点(Q_i代表速度控制多边形)
 * - 速度B样条 → 加速度控制点
 * - 加速度B样条 → 加加速度控制点
 * 
 * 【注意】导数控制点数量比原控制点少1个
 *   原: n+1个控制点 → 导数: n个控制点
 * 
 * @return 导数控制点矩阵 (维度 × (n-1))
 */
Eigen::MatrixXd UniformBspline::getDerivativeControlPoints() {
  // 导数控制点数量比原控制点少1个
  Eigen::MatrixXd ctp(control_points_.rows(), control_points_.cols() - 1);
  for (int i = 0; i < ctp.cols(); ++i) {
    // Qi = p * (P_{i+1} - P_i) / (u_{i+p+1} - u_{i+1})
    ctp.col(i) =
        p_ * (control_points_.col(i + 1) - control_points_.col(i)) / (u_(i + p_ + 1) - u_(i + 1));
  }
  return ctp;
}

/**
 * @brief 获取B样条的导数曲线(如: 位置→速度, 速度→加速度)
 * 
 * 【使用示例】
 *   UniformBspline pos_traj = ...;        // 位置轨迹
 *   UniformBspline vel_traj = pos_traj.getDerivative();   // 速度轨迹
 *   UniformBspline acc_traj = vel_traj.getDerivative();   // 加速度轨迹
 *   UniformBspline jerk_traj = acc_traj.getDerivative();  // 加加速度轨迹
 * 
 * 在Ego算法中,常使用:
 *   - 三次B样条(p=3)表示位置
 *   - 求导一次得二次B样条(p=2)表示速度
 *   - 再求导得一次B样条(p=1)表示加速度
 *   - 三次求导得零次B样条(p=0)表示加加速度
 * 
 * 【注意】由于clamped B样条的节点向量首尾重复,
 * 求导时需要移除首尾各一个节点以保证导数B样条的正确性。
 * 
 * @return 导数B样条曲线对象
 */
UniformBspline UniformBspline::getDerivative() {
  Eigen::MatrixXd   ctp = getDerivativeControlPoints();
  UniformBspline derivative(ctp, p_ - 1, interval_);

  /* 裁剪首尾节点: 导数B样条的节点向量比原向量少2个 */
  Eigen::VectorXd knot(u_.rows() - 2);
  knot = u_.segment(1, u_.rows() - 2);
  derivative.setKnot(knot);

  return derivative;
}

double UniformBspline::getInterval() { return interval_; }

/**
 * @brief 设置物理速度/加速度限制
 * 
 * 这些限制用于可行性检查,确保生成的轨迹满足无人机的动力学约束。
 * feasibility_tolerance_ 允许有微小超额(容差),避免过于严格的约束导致优化失败。
 * 
 * @param vel 最大速度限制 (m/s)
 * @param acc 最大加速度限制 (m/s²)
 * @param tolerance 容差比例 (如0.1表示允许超出10%)
 */
void UniformBspline::setPhysicalLimits(const double& vel, const double& acc, const double& tolerance) {
  limit_vel_   = vel;
  limit_acc_   = acc;
  limit_ratio_ = 1.1;
  feasibility_tolerance_ = tolerance;
}

/**
 * @brief 检查B样条轨迹的可行性(速度和加速度是否在限制范围内)
 * 
 * 【理论基础: B样条凸包性质】
 * B样条的一个重要性质是凸包性: B样条曲线完全位于其控制点构成的凸包内。
 * 这意味着:
 * - 速度曲线上的任意点,都在"速度控制点"构成的多边形凸包内
 * - 加速度曲线同理
 * 
 * 因此,我们只需要检查"速度控制点"和"加速度控制点"是否在限制范围内,
 * 就足以保证整条轨迹上的速度和加速度都满足约束。
 * 这大大简化了可行性检查的计算复杂度(O(n) vs O(n/Δt))。
 * 
 * 【检查方法】
 * 1. 速度检查: 遍历每个速度控制点,检查各轴分量是否超过速度限制
 *    速度控制点 = p * (P_{i+1} - P_i) / (u_{i+p+1} - u_{i+1})
 *    
 * 2. 加速度检查: 遍历每个加速度控制点,检查各轴分量是否超过加速度限制
 *    加速度控制点 = p*(p-1)*((P_{i+2}-P_{i+1})/(u_{i+p+2}-u_{i+2}) -
 *                           (P_{i+1}-P_i)/(u_{i+p+1}-u_{i+1})) / (u_{i+p+1}-u_{i+2})
 * 
 * 【输出】
 * - ratio: 超出最大限制的比例,用于后续时间延长(见lengthenTime)
 *   ratio = max( max_vel/limit_vel_, sqrt(max_acc/limit_acc_) )
 *   开方是因为加速度和时间是二次关系,延长t倍时间→加速度降为1/t²
 * 
 * @param ratio 输出: 超出比例,用于时间延长
 * @param show 是否打印详细信息
 * @return true 可行; false 不可行(存在超限)
 */
bool UniformBspline::checkFeasibility(double& ratio, bool show) {
  bool fea = true;

  Eigen::MatrixXd P         = control_points_;
  int             dimension = control_points_.rows();

  /* ---- 速度可行性检查 ---- */
  double max_vel = -1.0;
  double enlarged_vel_lim = limit_vel_ * (1.0 + feasibility_tolerance_) + 1e-4;
  for (int i = 0; i < P.cols() - 1; ++i) {
    // 计算速度控制点 Qi = p * (P_{i+1} - P_i) / (u_{i+p+1} - u_{i+1})
    Eigen::VectorXd vel = p_ * (P.col(i + 1) - P.col(i)) / (u_(i + p_ + 1) - u_(i + 1));

    // 检查各轴分量(无人机通常xyz三轴独立限速)
    if (fabs(vel(0)) > enlarged_vel_lim || fabs(vel(1)) > enlarged_vel_lim ||
        fabs(vel(2)) > enlarged_vel_lim) {

      if (show) cout << "[Check]: Infeasible vel " << i << " :" << vel.transpose() << endl;
      fea = false;

      // 记录最大速度(用于计算ratio)
      for (int j = 0; j < dimension; ++j) {
        max_vel = max(max_vel, fabs(vel(j)));
      }
    }
  }

  /* ---- 加速度可行性检查 ---- */
  double max_acc = -1.0;
  double enlarged_acc_lim = limit_acc_ * (1.0 + feasibility_tolerance_) + 1e-4;
  for (int i = 0; i < P.cols() - 2; ++i) {

    // 计算加速度控制点
    // 加速度 = 速度的导数,由三个相邻位置控制点通过二阶差分得到
    Eigen::VectorXd acc = p_ * (p_ - 1) *
        ((P.col(i + 2) - P.col(i + 1)) / (u_(i + p_ + 2) - u_(i + 2)) -
         (P.col(i + 1) - P.col(i)) / (u_(i + p_ + 1) - u_(i + 1))) /
        (u_(i + p_ + 1) - u_(i + 2));

    if (fabs(acc(0)) > enlarged_acc_lim || fabs(acc(1)) > enlarged_acc_lim ||
        fabs(acc(2)) > enlarged_acc_lim) {

      if (show) cout << "[Check]: Infeasible acc " << i << " :" << acc.transpose() << endl;
      fea = false;

      for (int j = 0; j < dimension; ++j) {
        max_acc = max(max_acc, fabs(acc(j)));
      }
    }
  }

  // 计算超出比例: 速度直接取最大比值,加速度取平方根
  // 因为时间缩放: 速度∝1/t, 加速度∝1/t², 所以加速度需要开方才能与速度统一
  ratio = max(max_vel / limit_vel_, sqrt(fabs(max_acc) / limit_acc_));

  return fea;
}

/**
 * @brief 延长轨迹时间以降低速度/加速度(使轨迹满足可行性)
 * 
 * 【原理说明】
 * 当轨迹的速度或加速度超出物理限制时,最简单有效的方法就是
 * 延长轨迹的执行时间。时间延长后,相同的位移在更长时间内完成,
 * 速度和加速度自然降低。
 * 
 * 定量的关系(假设均匀缩放):
 * - 时间延长为原来的 ratio 倍
 * - 速度变为原来的 1/ratio
 * - 加速度变为原来的 1/ratio²
 * 
 * 【实现方法】
 * 保持控制点位置不变,仅拉伸节点向量:
 * 1. 选取节点向量中间段(避开首尾重复节点,从索引5到倒数第5)
 * 2. 计算需要增加的总时间 delta_t = (ratio-1) * (u_end - u_start)
 * 3. 在中间段内,每个节点递增增加 t_inc 时间
 * 4. 在结束段(最后一个中间节点之后),统一增加 delta_t
 * 
 * 这样做的原因是:
 * - 保持起始段节点不变(确保起点位置/速度不变)
 * - 均匀拉伸中间段(相当于时间轴的线性拉伸)
 * 
 * @param ratio 时间延长倍率(>1)
 */
void UniformBspline::lengthenTime(const double& ratio) {
  // 避开首尾重复节点区域,选择中间段进行拉伸
  int num1 = 5;                          // 起始索引(避开前p+1个重复节点)
  int num2 = getKnot().rows() - 1 - 5;   // 结束索引(避开后p+1个重复节点)

  double delta_t = (ratio - 1.0) * (u_(num2) - u_(num1));  // 需要增加的总时间
  double t_inc   = delta_t / double(num2 - num1);           // 每个节点的增量
  // 中间段: 均匀递增增加
  for (int i = num1 + 1; i <= num2; ++i) u_(i) += double(i - num1) * t_inc;
  // 结束段: 统一增加相同时间量
  for (int i = num2 + 1; i < u_.rows(); ++i) u_(i) += delta_t;
}

// void UniformBspline::recomputeInit() {}

/**
 * @brief 将离散路径点参数化为B样条控制点(核心接口: 离散→连续)
 * 
 * 【功能概述】
 * 这是Ego算法中将A*等搜索得到的离散路径转化为连续光滑B样条轨迹的关键函数。
 * 输入一系列路径点(如A*路径)和边界导数条件,通过求解线性方程组得到B样条控制点。
 * 
 * 【数学原理: 三次均匀B样条离散化】
 * 对于三次均匀B样条(p=3,均匀节点),在节点中点处的求值公式为:
 *   C(u_k + interval/2) = (P_k + 4*P_{k+1} + P_{k+2}) / 6
 * 即: 相邻三个控制点的加权平均,权重比为1:4:1
 * 
 * 反过来,如果已知路径点 Q_k (希望曲线经过的点),可以构造方程:
 *   (P_k + 4*P_{k+1} + P_{k+2}) / 6 = Q_k,  k=0,...,K-1
 * 其中K是路径点数量,K+2是控制点数量。
 * 
 * 【约束系统构造】
 * 完整约束包括(K+4个方程,求解K+2个控制点):
 * 
 * 1. 位置约束(K个): 每个路径点对应一个方程
 *    (P_i + 4*P_{i+1} + P_{i+2}) / 6 = Q_i
 *    
 * 2. 起点速度约束(1个): 
 *    (P_1 - P_0) / (2*ts) = v_start
 *    → -P_0 + P_1 = 2*ts*v_start
 *    
 * 3. 终点速度约束(1个):
 *    (P_{K+1} - P_K) / (2*ts) = v_end  
 *    → -P_{K-1} + P_K = 2*ts*v_end
 *    
 * 4. 起点加速度约束(1个):
 *    (P_0 - 2*P_1 + P_2) / ts² = a_start
 *    
 * 5. 终点加速度约束(1个):
 *    (P_{K-1} - 2*P_K + P_{K+1}) / ts² = a_end
 * 
 * 写成矩阵形式:
 *   A * P = b
 * 其中A是(K+4)×(K+2)矩阵,P是(K+2)个控制点,b是已知条件
 * 超定系统,使用QR分解求最小二乘解。
 * 
 * @param ts 时间步长(节点间隔)
 * @param point_set 离散路径点序列
 * @param start_end_derivative 边界导数[v0, vf, a0, af]
 * @param ctrl_pts 输出: B样条控制点矩阵
 */
void UniformBspline::parameterizeToBspline(const double& ts, const vector<Eigen::Vector3d>& point_set,
                                              const vector<Eigen::Vector3d>& start_end_derivative,
                                              Eigen::MatrixXd&               ctrl_pts) {
  if (ts <= 0) {
    cout << "[B-spline]:time step error." << endl;
    return;
  }

  if (point_set.size() <= 3) {
    cout << "[B-spline]:point set have only " << point_set.size() << " points." << endl;
    return;
  }

  if (start_end_derivative.size() != 4) {
    cout << "[B-spline]:derivatives error." << endl;
  }

  int K = point_set.size();  // 路径点数量

  // ---- 构造约束矩阵A ----
  // 各行的权重模板
  Eigen::Vector3d prow(3), vrow(3), arow(3);
  prow << 1, 4, 1;                      // 位置: [1,4,1] /6
  vrow << -1, 0, 1;                     // 速度: [-1,0,1] /(2*ts)
  arow << 1, -2, 1;                     // 加速度: [1,-2,1] /ts²

  // A矩阵: (K+4)行 × (K+2)列
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(K + 4, K + 2);

  // 位置约束行(前K行): 每行对应一个路径点的位置方程
  for (int i = 0; i < K; ++i) 
    A.block(i, i, 1, 3) = (1 / 6.0) * prow.transpose();

  // 速度约束行: 起点速度(第K行)和终点速度(第K+1行)
  A.block(K, 0, 1, 3)         = (1 / 2.0 / ts) * vrow.transpose();
  A.block(K + 1, K - 1, 1, 3) = (1 / 2.0 / ts) * vrow.transpose();

  // 加速度约束行: 起点加速度(第K+2行)和终点加速度(第K+3行)
  A.block(K + 2, 0, 1, 3)     = (1 / ts / ts) * arow.transpose();
  A.block(K + 3, K - 1, 1, 3) = (1 / ts / ts) * arow.transpose();

  // ---- 构造右侧向量b(x/y/z三轴独立求解) ----
  Eigen::VectorXd bx(K + 4), by(K + 4), bz(K + 4);
  // 位置分量
  for (int i = 0; i < K; ++i) {
    bx(i) = point_set[i](0);
    by(i) = point_set[i](1);
    bz(i) = point_set[i](2);
  }
  // 边界导数分量
  for (int i = 0; i < 4; ++i) {
    bx(K + i) = start_end_derivative[i](0);
    by(K + i) = start_end_derivative[i](1);
    bz(K + i) = start_end_derivative[i](2);
  }

  // ---- 求解线性方程组 Ax = b (QR分解,最小二乘) ----
  Eigen::VectorXd px = A.colPivHouseholderQr().solve(bx);
  Eigen::VectorXd py = A.colPivHouseholderQr().solve(by);
  Eigen::VectorXd pz = A.colPivHouseholderQr().solve(bz);

  // ---- 合并为控制点矩阵 ----
  ctrl_pts.resize(3, K + 2);
  ctrl_pts.row(0) = px.transpose();
  ctrl_pts.row(1) = py.transpose();
  ctrl_pts.row(2) = pz.transpose();
}

/**
 * @brief 获取轨迹总持续时间
 * 
 * 对于均匀B样条,总时长 = u_(m_-p_) - u_(p_) = (n_-p_) * interval_
 * 这是B样条曲线定义域的总长度,即无人机飞完整条轨迹所需的时间。
 * 
 * @return 总时间长度(秒),失败返回-1.0
 */
double UniformBspline::getTimeSum() {
  double tm, tmp;
  if ( getTimeSpan(tm, tmp) )
    return tmp - tm;
  else
    return -1.0;
  
}

/**
 * @brief 计算轨迹的近似弧长
 * 
 * 【方法】
 * 对轨迹曲线进行密集采样,累加相邻采样点之间的欧氏距离。
 * 采样分辨率控制精度,默认0.01秒间隔。
 * 这是一种数值积分近似,精度随分辨率提高而提高。
 * 
 * @param res 采样时间分辨率(秒)
 * @return 轨迹总长度(米)
 */
double UniformBspline::getLength(const double& res) {
  double          length = 0.0;
  double          dur    = getTimeSum();
  Eigen::VectorXd p_l    = evaluateDeBoorT(0.0), p_n;
  for (double t = res; t <= dur + 1e-4; t += res) {
    p_n = evaluateDeBoorT(t);
    length += (p_n - p_l).norm();  // 累加相邻采样点的距离
    p_l = p_n;
  }
  return length;
}

/**
 * @brief 计算轨迹的加加速度积分(平滑度评价指标)
 * 
 * 【加加速度(Jerk)】
 * 加加速度是位置的三阶导数,即加速度的变化率。
 * 对于无人机轨迹,加加速度积分 ∫(jerk(t))²dt 是衡量轨迹平滑度的标准指标:
 * - 较小的jerk积分 → 轨迹平滑,加速度变化缓慢
 * - 较大的jerk积分 → 轨迹抖动,加速度变化剧烈
 * 
 * 【计算方法】
 * 1. 对位置B样条求三次导得到加加速度B样条(0阶)
 * 2. 加加速度B样条的控制点就是各时间段内的jerk值
 * 3. jerk积分 = Σ (time_span_i) * (jerk_i)²
 *    其中 time_span_i 是每个控制点对应的时间区间长度
 * 
 * 【注意】这里计算的其实是jerk²的积分,不是纯jerk积分
 * 在优化中通常使用jerk²积分作为代价函数,因为:
 * - 正负jerk都产生代价
 * - 数学上可导,便于梯度优化
 * 
 * @return 加加速度平方积分值
 */
double UniformBspline::getJerk() {
  // 求导三次: 位置→速度→加速度→加加速度(jerk)
  UniformBspline jerk_traj = getDerivative().getDerivative().getDerivative();

  Eigen::VectorXd times     = jerk_traj.getKnot();
  Eigen::MatrixXd ctrl_pts  = jerk_traj.getControlPoint();
  int             dimension = ctrl_pts.rows();

  double jerk = 0.0;
  for (int i = 0; i < ctrl_pts.cols(); ++i) {
    for (int j = 0; j < dimension; ++j) {
      // jerk²积分 = Σ Δt_i * jerk_i²
      jerk += (times(i + 1) - times(i)) * ctrl_pts(j, i) * ctrl_pts(j, i);
    }
  }

  return jerk;
}

/**
 * @brief 获取轨迹的平均速度和最大速度(用于性能评估)
 * 
 * 对速度B样条进行密集采样,统计速度的均值和最大值。
 * 采样间隔0.01秒,结果可用于评估轨迹是否满足速度约束。
 * 
 * @param mean_v 输出: 平均速度 (m/s)
 * @param max_v 输出: 最大速度 (m/s)
 */
void UniformBspline::getMeanAndMaxVel(double& mean_v, double& max_v) {
  UniformBspline vel = getDerivative();  // 获取速度B样条
  double            tm, tmp;
  vel.getTimeSpan(tm, tmp);

  double max_vel = -1.0, mean_vel = 0.0;
  int    num = 0;
  for (double t = tm; t <= tmp; t += 0.01) {  // 10ms采样间隔
    Eigen::VectorXd vxd = vel.evaluateDeBoor(t);
    double          vn  = vxd.norm();  // 速度幅值

    mean_vel += vn;
    ++num;
    if (vn > max_vel) {
      max_vel = vn;
    }
  }

  mean_vel = mean_vel / double(num);
  mean_v   = mean_vel;
  max_v    = max_vel;
}

/**
 * @brief 获取轨迹的平均加速度和最大加速度(用于性能评估)
 * 
 * 对加速度B样条进行密集采样,统计加速度的均值和最大值。
 * 可用于评估轨迹的动力学可行性。
 * 
 * @param mean_a 输出: 平均加速度 (m/s²)
 * @param max_a 输出: 最大加速度 (m/s²)
 */
void UniformBspline::getMeanAndMaxAcc(double& mean_a, double& max_a) {
  UniformBspline acc = getDerivative().getDerivative();  // 获取加速度B样条
  double            tm, tmp;
  acc.getTimeSpan(tm, tmp);

  double max_acc = -1.0, mean_acc = 0.0;
  int    num = 0;
  for (double t = tm; t <= tmp; t += 0.01) {
    Eigen::VectorXd axd = acc.evaluateDeBoor(t);
    double          an  = axd.norm();

    mean_acc += an;
    ++num;
    if (an > max_acc) {
      max_acc = an;
    }
  }

  mean_acc = mean_acc / double(num);
  mean_a   = mean_acc;
  max_a    = max_acc;
}
}  // namespace rebound_planner