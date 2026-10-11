/**
 * @file gradient_descent_optimizer.cpp
 * @brief 梯度下降优化器实现
 *
 * 【算法详解】
 * 
 * 实现了具有BB(Barzilai-Borwein)步长和Armijo线搜索的梯度下降算法。
 * 
 * BB步长:
 *   α_BB = s·y / y·y
 *   其中 s = x_{k+1} - x_k, y = ∇f_{k+1} - ∇f_k
 *   BB步长利用最近两次迭代的信息估计Hessian的逆,比固定步长收敛更快
 * 
 * Armijo线搜索条件(充分下降条件):
 *   f(x_k - α∇f_k) ≤ f(x_k) - c·α·||∇f_k||²
 *   其中 c = 1e-4(标准值)
 *   如果不满足条件,则α *= 0.5(回溯)
 * 
 * 交替更新策略:
 *   奇偶迭代交替更新x_k和x_{k+1},避免数据拷贝开销
 */

#include <bspline_opt/gradient_descent_optimizer.h>

#define RESET   "\033[0m"
#define RED     "\033[31m"

/**
 * @brief 执行梯度下降优化
 * 
 * 【算法流程】
 * 1. 初始化: 计算初始梯度,使用限制步长做第一次迭代
 * 2. 主循环:
 *    a. 计算BB步长 α = s·y / y·y
 *    b. 沿负梯度方向更新: x_{new} = x_{old} - α·∇f(x_{old})
 *    c. Armijo线搜索: 如果不满足充分下降条件,回溯减小步长(α *= 0.5)
 *    d. 检查收敛: 如果梯度模长 < min_grad_, 返回FIND_MIN
 *    e. 检查force_return: 如果目标函数设置force_return, 返回RETURN_BY_ORDER
 * 3. 达到最大迭代次数返回REACH_MAX_ITERATION
 * 
 * 奇偶交替优化:
 *   偶数迭代: x_k不变,更新x_{kp1}
 *   奇数迭代: x_{kp1}不变,更新x_k
 *   这样可以避免不必要的内存拷贝
 * 
 * @param x_init_optimal 初始变量值,优化后保存最优值
 * @param opt_f 输出的最优代价值
 * @return 优化结果
 */
GradientDescentOptimizer::RESULT 
GradientDescentOptimizer::optimize(Eigen::VectorXd &x_init_optimal, double &opt_f)
{
    /* 参数有效性检查 */
    if (min_grad_ < 1e-10) 
        { cout << RED << "min_grad_ is invalid:" << min_grad_ << RESET << endl; return FAILED;}
    if (iter_limit_ <= 2)  
        { cout << RED  << "iter_limit_ is invalid:" << iter_limit_ << RESET << endl; return FAILED;}

    void* f_data = f_data_;
    int iter = 2;                    // 已迭代次数计数(已做了2次初始评估)
    int invoke_count = 2;            // 函数调用计数
    bool force_return;               // 强制返回标志
    Eigen::VectorXd x_k( x_init_optimal ), x_kp1( x_init_optimal.rows() );  // 当前和下一步变量
    double cost_km1, cost_k, cost_kp1, cost_min;
    Eigen::VectorXd grad_k( x_init_optimal.rows() ), grad_kp1( x_init_optimal.rows() );

    /*** 第一步:计算初始点的代价和梯度 ***/
    cost_k = objfun_(x_k, grad_k, force_return, f_data);
    if ( force_return ) return RETURN_BY_ORDER;
    cost_min = cost_k;

    /*** 第二步:计算第一次迭代的步长(限制最大移动距离) ***/
    double max_grad = max( abs(grad_k.maxCoeff()), abs(grad_k.minCoeff()) );
    constexpr double MAX_MOVEMENT_AT_FIRST_ITERATION = 0.1; // 首次迭代最大移动距离(米)
    double alpha0 = max_grad < MAX_MOVEMENT_AT_FIRST_ITERATION ?
        1.0 : (MAX_MOVEMENT_AT_FIRST_ITERATION/max_grad);
    
    /*** 第三步:执行第一次迭代 ***/
    x_kp1 = x_k - alpha0 * grad_k;
    cost_kp1 = objfun_(x_kp1, grad_kp1, force_return, f_data);
    if ( force_return ) return RETURN_BY_ORDER;
    if ( cost_min > cost_kp1 ) cost_min = cost_kp1;

    /*** 主迭代循环开始 ***/
    while( ++iter <= iter_limit_ && invoke_count <= invoke_limit_ )
    {
        /*** 计算BB步长: α = s·y / y·y ***/
        Eigen::VectorXd s = x_kp1 - x_k;       // 变量变化
        Eigen::VectorXd y = grad_kp1 - grad_k;  // 梯度变化
        double alpha = s.dot(y) / y.dot(y);      // BB步长
        if ( isnan(alpha) || isinf(alpha) )
            { cout << RED << "step size invalid! alpha=" << alpha << RESET << endl; return FAILED;}

        /**
         * 奇偶交替更新策略:
         * 偶数迭代: 固定x_k, 更新x_kp1
         * 奇数迭代: 固定x_kp1, 更新x_k
         * 这样可以减少变量交换的开销
         */
        if (iter % 2)  // 奇数迭代: 固定x_kp1, 更新x_k
        {
            do
            {
                x_k = x_kp1 - alpha * grad_kp1;
                cost_k = objfun_(x_k, grad_k, force_return, f_data);
                invoke_count++;
                if ( force_return ) return RETURN_BY_ORDER;
                alpha *= 0.5;  // Armijo回溯: 步长减半
            } while ( cost_k > cost_kp1 - 1e-4*alpha*grad_kp1.transpose()*grad_kp1 ); // Armijo条件

            /* 检查梯度是否收敛 */
            if ( grad_k.norm() < min_grad_ )
            {
                opt_f = cost_k;
                return FIND_MIN;
            }
        }
        else  // 偶数迭代: 固定x_k, 更新x_kp1
        {
            do
            {
                x_kp1 = x_k - alpha * grad_k;
                cost_kp1 = objfun_(x_kp1, grad_kp1, force_return, f_data);
                invoke_count++;
                if ( force_return ) return RETURN_BY_ORDER;
                alpha *= 0.5;  // Armijo回溯: 步长减半
            } while ( cost_kp1 > cost_k - 1e-4*alpha*grad_k.transpose()*grad_k ); // Armijo条件
            
            /* 检查梯度是否收敛 */
            if ( grad_kp1.norm() < min_grad_ )
            {
                opt_f = cost_kp1;
                return FIND_MIN;
            }
        }

    }

    /*** 达到最大迭代次数,返回当前最优值 ***/
    opt_f = iter_limit_ % 2 ? cost_k : cost_kp1;
    return REACH_MAX_ITERATION;
}