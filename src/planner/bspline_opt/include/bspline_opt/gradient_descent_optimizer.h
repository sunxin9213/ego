/**
 * @file gradient_descent_optimizer.h
 * @brief 梯度下降优化器头文件
 *
 * 【算法说明】
 * 实现了一个带有Armijo线搜索的梯度下降优化器。
 * 用于B样条轨迹优化中的无约束非线性优化问题。
 * 
 * 核心特点:
 * 1. 使用BB(Barzilai-Borwein)步长自动计算初始步长
 * 2. Armijo条件确保充分下降
 * 3. 支持通过force_return标志提前终止(用于反弹机制)
 * 4. 支持最小梯度模长收敛判定
 */

#ifndef _GRADIENT_DESCENT_OPT_H_
#define _GRADIENT_DESCENT_OPT_H_

#include <iostream>
#include <vector>
#include <Eigen/Eigen>

using namespace std;

/**
 * @brief 梯度下降优化器类
 *
 * 实现梯度下降优化算法,用于B样条轨迹优化。
 * 使用BB方法估计步长,结合Armijo线搜索保证收敛。
 */
class GradientDescentOptimizer
{

public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  /**
   * @brief 目标函数指针类型
   * @param x 优化变量
   * @param grad 输出的梯度
   * @param force_return 是否强制提前返回(用于反弹触发)
   * @param data 用户数据指针
   * @return 代价值
   */
  typedef double (*objfunDef)(const Eigen::VectorXd &x, Eigen::VectorXd &grad, bool& force_return, void* data);
  
  /**
   * @brief 优化结果枚举
   */
  enum RESULT
  {
    FIND_MIN,              // 找到最小值(梯度满足收敛条件)
    FAILED,                // 优化失败
    RETURN_BY_ORDER,       // 被强制返回(检测到新障碍物触发反弹)
    REACH_MAX_ITERATION    // 达到最大迭代次数
  };
    
  /**
   * @brief 构造函数
   * @param v_num 优化变量维度
   * @param objf 目标函数指针
   * @param f_data 目标函数所需的数据
   */
  GradientDescentOptimizer(int v_num, objfunDef objf, void* f_data)
  {
    variable_num_ = v_num;
    objfun_ = objf;
    f_data_ = f_data;
  };

  /* 参数设置接口 */
  void set_maxiter(int limit){ iter_limit_ = limit; }         // 设置最大迭代次数
  void set_maxeval(int limit){ invoke_limit_ = limit; }       // 设置最大函数调用次数
  void set_xtol_rel(double xtol_rel){ xtol_rel_ = xtol_rel; } // 设置相对变量容差
  void set_xtol_abs(double xtol_abs){ xtol_abs_ = xtol_abs; } // 设置绝对变量容差
  void set_min_grad(double min_grad){ min_grad_ = min_grad; } // 设置最小梯度模长(收敛条件)

  /**
   * @brief 执行优化
   * @param x_init_optimal 初始值/最优值(输入初值,输出最优值)
   * @param opt_f 输出的最优代价值
   * @return 优化结果状态
   */
  RESULT optimize(Eigen::VectorXd &x_init_optimal, double &opt_f);

  
private:

  int variable_num_{0};      // 变量维度
  int iter_limit_{1e10};     // 最大迭代次数
  int invoke_limit_{1e10};   // 最大函数调用次数
  double xtol_rel_;          // 相对变量变化容差
  double xtol_abs_;          // 绝对变量变化容差
  double min_grad_;          // 最小梯度模长(用于收敛判定)
  double time_limit_;        // 时间限制
  void* f_data_;             // 用户数据指针
  objfunDef objfun_;         // 目标函数指针

};


#endif