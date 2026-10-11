/**
 * @file bspline_optimizer.cpp
 * @brief B样条轨迹优化器实现
 *
 * 【Ego Planner核心优化模块】
 * 
 * 该文件实现了基于"弹性带(Elastic Band)"和梯度下降的B样条轨迹优化算法。
 * 核心思想:
 *   1. 使用B样条表示轨迹,通过优化控制点来获得平滑、安全、动力学可行的轨迹
 *   2. "Rebound(反弹)"机制:当检测到碰撞时,通过A*搜索生成引导路径,
 *      然后优化控制点"弹回"到安全区域
 *   3. 优化过程综合考虑:平滑度、碰撞距离、动力学可行性、路径拟合
 *
 * 【优化变量】
 *   优化变量是B样条内部的(N-2*order)个控制点(3维位置),
 *   固定首尾各order个控制点以保证边界条件不变。
 *
 * 【代价函数组成】
 *   J_total = λ1*J_smooth + λ2*J_distance + λ3*J_feasibility + λ4*J_fitness
 *   - J_smooth: 平滑度代价,最小化加加速度(jerk)
 *   - J_distance: 距离代价,保持与障碍物的安全距离
 *   - J_feasibility: 可行性代价,满足速度和加速度约束
 *   - J_fitness: 拟合代价,追踪参考路径
 */

#include "bspline_opt/bspline_optimizer.h"
#include "bspline_opt/gradient_descent_optimizer.h"
// using namespace std;

namespace rebound_planner {

void BsplineOptimizer::setParam(ros::NodeHandle& nh) {
  /*** 从ROS参数服务器读取优化权重和约束参数 ***/

  /* 平滑度权重:控制轨迹的平滑程度,越大越平滑 */
  nh.param("optimization/lambda_smooth", lambda1_, -1.0);
  /* 碰撞代价权重:控制避障的强度,越大离障碍物越远 */
  nh.param("optimization/lambda_collision", lambda2_, -1.0);
  /* 可行性权重:控制速度和加速度约束的满足程度 */
  nh.param("optimization/lambda_feasibility", lambda3_, -1.0);
  /* 拟合权重:控制轨迹对参考路径的追踪程度 */
  nh.param("optimization/lambda_fitness", lambda4_, -1.0);

  /* 安全距离:控制点与障碍物之间期望保持的最小距离 */
  nh.param("optimization/dist0", dist0_, -1.0);
  /* 最大速度限制(m/s) */
  nh.param("optimization/max_vel", max_vel_, -1.0);
  /* 最大加速度限制(m/s²) */
  nh.param("optimization/max_acc", max_acc_, -1.0);
  /* B样条阶数(默认3阶,即三次B样条,保证加加速度连续) */
  nh.param("optimization/order", order_, 3);
}

void BsplineOptimizer::setEnvironment(const GridMap::Ptr& env) {
  this->grid_map_ = env;
}

void BsplineOptimizer::setControlPoints(const Eigen::MatrixXd& points) {
  cps_.points = points;
}

void BsplineOptimizer::setBsplineInterval(const double& ts) { bspline_interval_ = ts; }

/**
 * @brief 初始化控制点(弹性带/反弹优化的核心预处理步骤)
 * 
 * 【算法流程】
 * 步骤1: 碰撞检测与分段
 *   - 沿初始轨迹逐点检查是否被障碍物占据
 *   - 将轨迹分为多个"碰撞段"(segment),每段由碰撞边界in_id/out_id标记
 * 
 * 步骤2: A*搜索
 *   - 对每个碰撞段,使用A*算法搜索一条绕过障碍物的无碰撞路径
 *   - A*路径将作为后续优化的引导
 * 
 * 步骤3: 计算边界约束
 *   - 确定每段中可以被优化的控制点索引范围(bounds)
 *   - 段与段之间保持一定的重叠区域,确保连续性
 * 
 * 步骤4: 生成基点和方向向量(核心弹性带机制)
 *   - 对每个碰撞段中的控制点,计算其在A*路径上的投影点(交点)
 *   - 生成从控制点指向A*路径的方向向量
 *   - 这些方向向量定义了控制点应移动的方向和距离
 *   - 相邻控制点的方向向量会向内/外传播,形成"弹性带"效果
 * 
 * @param init_points 初始控制点矩阵(3行,N列)
 * @param flag_first_init 是否首次初始化(首次需重置控制点数据)
 * @return A*搜索路径集合(用于后续反弹优化)
 */
std::vector<std::vector<Eigen::Vector3d>> BsplineOptimizer::initControlPoints(Eigen::MatrixXd& init_points, bool flag_first_init /*= true*/)
{
  if ( flag_first_init )
  {
    /*** 首次初始化:设置安全距离,调整控制点大小,保存初始控制点 ***/
    cps_.clearance = dist0_;
    cps_.resize( init_points.cols() );
    cps_.points = init_points;
  }

  /*** 步骤1:根据障碍物分割初始轨迹 ***/
  constexpr int ENOUGH_INTERVAL = 2;  // 连续无障碍点的阈值数量
  double step_size = grid_map_->getResolution() / ( (init_points.col(0) - init_points.rightCols(1)).norm() / (init_points.cols()-1) ) / 2;
  int in_id, out_id;  // 碰撞段的入口和出口控制点索引
  vector<std::pair<int,int>> segment_ids;  // 存储所有碰撞段(入口,出口)对
  int same_occ_state_times = ENOUGH_INTERVAL + 1;
  bool occ, last_occ = false;
  bool flag_got_start = false, flag_got_end = false, flag_got_end_maybe = false;
  /* 只检查前2/3段的碰撞(靠近起点的部分更重要) */
  int i_end = (int)init_points.cols() - order_ - ((int)init_points.cols()-2*order_)/3;
  for ( int i=order_; i<=i_end; ++i )
  {
    /**
     * 沿相邻控制点连线方向以栅格分辨率为步长采样,
     * 检查是否有障碍物占据。
     * 使用线性插值确保不错过两个控制点之间的障碍物。
     */
    for ( double a=1.0; a>=0.0; a-=step_size )
    {
      occ = grid_map_->getInflateOccupancy(a * init_points.col(i-1) + (1-a) * init_points.col(i));

      /*** 状态转换检测:从自由空间进入障碍物 ***/
      if ( occ && ! last_occ)
      {
        if ( same_occ_state_times > ENOUGH_INTERVAL || i == order_ )
        {
          in_id = i-1;  // 记录碰撞段入口
          flag_got_start = true;
        }
        same_occ_state_times = 0;
        flag_got_end_maybe = false;
      }
      /*** 状态转换检测:从障碍物回到自由空间 ***/
      else if( !occ && last_occ )
      {
        out_id = i;  // 记录碰撞段出口
        flag_got_end_maybe = true;
        same_occ_state_times = 0;
      }
      else
      {
        ++ same_occ_state_times;  // 连续相同状态的计数
      }

      /*** 确认已经离开障碍物(连续ENOUGH_INTERVAL个点无障碍) ***/
      if ( flag_got_end_maybe && ( same_occ_state_times > ENOUGH_INTERVAL || ( i == (int)init_points.cols()-order_ ) ) )
      {
        flag_got_end_maybe = false;
        flag_got_end = true;
      }
      
      last_occ = occ;
      
      /*** 找到一个完整的碰撞段,记录其入口和出口索引 ***/
      if ( flag_got_start && flag_got_end )
      {
        flag_got_start = false;
        flag_got_end = false;
        segment_ids.push_back( std::pair<int,int>(in_id, out_id) );
      }
    }
  }


  /*** 步骤2:对每个碰撞段运行A*搜索,生成绕过障碍物的引导路径 ***/
  vector<vector<Eigen::Vector3d>> a_star_pathes;
  for ( size_t i=0; i<segment_ids.size(); ++i )
  {
    Eigen::Vector3d in( init_points.col(segment_ids[i].first) ), out( init_points.col(segment_ids[i].second) );
    /* 使用0.1m的搜索步长在3D栅格地图中进行A*搜索 */
    if ( a_star_->AstarSearch( 0.1 , in, out) )
    {
      a_star_pathes.push_back( a_star_->getPath() );
    }
    else
    {
      ROS_ERROR("a star error, force return!");
      return a_star_pathes;
    }
  } 

  /*** 步骤3:计算每个碰撞段的控制点边界约束 ***/
  int id_low_bound, id_up_bound;
  vector<std::pair<int,int>> bounds( segment_ids.size() );
  for (size_t i=0; i<segment_ids.size(); i++)
  {
    /**
     * 边界计算原则:
     * - 每个段的边界不能超出相邻段的中点
     * - 确保不同段之间不会重叠冲突
     * - 保留首尾各order个控制点不受优化影响(保持边界条件)
     */
    if ( i == 0 ) // 第一个段
    {
      id_low_bound = order_;
      if ( segment_ids.size() > 1 )
      {
        id_up_bound =  (int)(((segment_ids[0].second + segment_ids[1].first)-1.0f) / 2);
      }
      else
      {
        id_up_bound = init_points.cols() - order_ - 1;
      }
    }
    else if ( i == segment_ids.size() - 1 ) // 最后一个段
    {
      id_low_bound = (int)(((segment_ids[i].first + segment_ids[i-1].second)+1.0f) / 2);
      id_up_bound = init_points.cols() - order_ - 1;
    }
    else
    {
      id_low_bound = (int)(((segment_ids[i].first + segment_ids[i-1].second)+1.0f) / 2);
      id_up_bound = (int)(((segment_ids[i].second + segment_ids[i+1].first)-1.0f) / 2);
    }
    
    bounds[i] = std::pair<int,int>(id_low_bound, id_up_bound);
  }

  // cout << "+++++++++" << endl;
  // for ( int j=0; j<bounds.size(); ++j )
  // {
  //   cout << bounds[j].first << "  " << bounds[j].second << endl;
  // }

  /*** 调整段长度,确保每段有足够的控制点 ***/
  vector<std::pair<int,int>> final_segment_ids( segment_ids.size() );
  constexpr double MINIMUM_PERCENT = 0.0;  // 每段保证有足够点来产生足够的推力
  int minimum_points = round(init_points.cols() * MINIMUM_PERCENT), num_points;
  for (size_t i=0; i<segment_ids.size(); i++)
  {
    num_points = segment_ids[i].second - segment_ids[i].first + 1;
    if ( num_points < minimum_points )
    {
      double add_points_each_side = (int)(((minimum_points - num_points)+1.0f) / 2);

      final_segment_ids[i].first = segment_ids[i].first - add_points_each_side >= bounds[i].first ?
        segment_ids[i].first - add_points_each_side :
        bounds[i].first;

      final_segment_ids[i].second = segment_ids[i].second + add_points_each_side <= bounds[i].second ?
        segment_ids[i].second + add_points_each_side :
        bounds[i].second;
    }
    else
    {
      final_segment_ids[i].first = segment_ids[i].first;
      final_segment_ids[i].second = segment_ids[i].second;
    }
  }

  /*** 步骤4:为每个段分配基点和方向向量(弹性带机制核心) ***/
  for (size_t i=0; i<segment_ids.size(); i++)
  {
    // step 1: 初始化临时标志
    for ( int j=final_segment_ids[i].first; j <= final_segment_ids[i].second; ++j )
      cps_.flag_temp[j] = false;

    // step 2: 计算控制点在A*路径上的投影交点
    int got_intersection_id = -1;
    for ( int j=segment_ids[i].first+1; j<segment_ids[i].second; ++j )
    {
      /**
       * 核心算法:使用"二分法"在A*路径上搜索控制点的投影点
       * ctrl_pts_law 是相邻控制点的差分向量,定义了搜索方向
       * 通过点积的符号变化确定投影交点的位置
       * 找到交点后,从交点向控制点方向做射线,找到第一个安全点作为基点
       */
      Eigen::Vector3d ctrl_pts_law(cps_.points.col(j+1) - cps_.points.col(j-1)), intersection_point;
      int Astar_id = a_star_pathes[i].size() / 2, last_Astar_id;
      double val = (a_star_pathes[i][Astar_id] - cps_.points.col(j)).dot( ctrl_pts_law ), last_val = val;
      while ( Astar_id >=0 && Astar_id < (int)a_star_pathes[i].size() )
      {
        last_Astar_id = Astar_id;

        if ( val >= 0 )
          -- Astar_id;
        else
          ++ Astar_id;
        
        val = (a_star_pathes[i][Astar_id] - cps_.points.col(j)).dot( ctrl_pts_law );
        
        /* 检测到符号变化,说明找到了交点 */
        if ( val * last_val <= 0 && ( abs(val) > 0 || abs(last_val) > 0 ) )
        {
          /**
           * 使用线性插值计算精确的交点位置:
           * intersection_point = A*[id] + (A*[id] - A*[last_id]) * t
           * 其中 t = (ctrl_pts_law · (P_j - A*[id])) / (ctrl_pts_law · (A*[id] - A*[last_id]))
           */
          intersection_point = 
            a_star_pathes[i][Astar_id] + 
            ( ( a_star_pathes[i][Astar_id] - a_star_pathes[i][last_Astar_id] ) * 
              ( ctrl_pts_law.dot( cps_.points.col(j) - a_star_pathes[i][Astar_id] ) / ctrl_pts_law.dot( a_star_pathes[i][Astar_id] -  a_star_pathes[i][last_Astar_id] ) )
            );

          got_intersection_id = j;
          break;
        }
      }

      /**
       * 从交点到控制点的连线上,从交点向控制点移动,
       * 找到第一个安全(无碰撞)的位置作为"基点"。
       * 方向向量 = 从基点到交点的归一化方向。
       * 这个方向就是优化时控制点应移动的"弹性"方向。
       */
      if ( got_intersection_id >= 0 )
      {
        cps_.flag_temp[j] = true;
        double length = (intersection_point - cps_.points.col(j)).norm();
        if ( length > 1e-5 )
        {
          for ( double a=length; a>=0.0; a-=grid_map_->getResolution() )
          {
            occ =  grid_map_->getInflateOccupancy((a/length)*intersection_point + (1-a/length)*cps_.points.col(j));
      
            if ( occ || a < grid_map_->getResolution() )
            {
              if ( occ )
                a+=grid_map_->getResolution();
              cps_.base_point[j].push_back( (a/length)*intersection_point + (1-a/length)*cps_.points.col(j) );
              cps_.direction[j].push_back( (intersection_point - cps_.points.col(j)).normalized() );
              break;
            }
          }
        }
      }
    }

    /* 特殊情况:段长度太短时,控制点可能在A*路径外侧,需要特殊处理 */
    if ( segment_ids[i].second - segment_ids[i].first == 1 ) 
    {
      Eigen::Vector3d ctrl_pts_law(cps_.points.col(segment_ids[i].second) - cps_.points.col(segment_ids[i].first)), intersection_point;
      Eigen::Vector3d middle_point = (cps_.points.col(segment_ids[i].second) + cps_.points.col(segment_ids[i].first)) / 2;
      int Astar_id = a_star_pathes[i].size() / 2, last_Astar_id;
      double val = (a_star_pathes[i][Astar_id] - middle_point).dot( ctrl_pts_law ), last_val = val;
      while ( Astar_id >=0 && Astar_id < (int)a_star_pathes[i].size() )
      {
        last_Astar_id = Astar_id;

        if ( val >= 0 )
          -- Astar_id;
        else
          ++ Astar_id;
        
        val = (a_star_pathes[i][Astar_id] - middle_point).dot( ctrl_pts_law );
        
        if ( val * last_val <= 0 && ( abs(val) > 0 || abs(last_val) > 0 ) )
        {
          intersection_point = 
            a_star_pathes[i][Astar_id] + 
            ( ( a_star_pathes[i][Astar_id] - a_star_pathes[i][last_Astar_id] ) * 
              ( ctrl_pts_law.dot( middle_point - a_star_pathes[i][Astar_id] ) / ctrl_pts_law.dot( a_star_pathes[i][Astar_id] -  a_star_pathes[i][last_Astar_id] ) )
            );

          if ( (intersection_point - middle_point).norm() > 0.01 ) // 1cm阈值
          {
            cps_.flag_temp[segment_ids[i].first] = true;
            cps_.base_point[segment_ids[i].first].push_back( cps_.points.col(segment_ids[i].first) );
            cps_.direction[segment_ids[i].first].push_back( (intersection_point - middle_point).normalized() );

            got_intersection_id = segment_ids[i].first;
          }
          break;
        }
      }
    }

    // step 3: 将已计算的基点和方向传播到相邻的未标记控制点
    if ( got_intersection_id >= 0 )
    {
      for ( int j=got_intersection_id + 1; j <= final_segment_ids[i].second; ++j )
        if ( ! cps_.flag_temp[j] )
        {
          cps_.base_point[j].push_back( cps_.base_point[j-1].back() );
          cps_.direction[j].push_back( cps_.direction[j-1].back() );
        }

      for ( int j=got_intersection_id - 1; j >= final_segment_ids[i].first; --j )
        if ( ! cps_.flag_temp[j] )
        {
          cps_.base_point[j].push_back( cps_.base_point[j+1].back() );
          cps_.direction[j].push_back( cps_.direction[j+1].back() );
        }
    }
    else
    {
      // 如果找不到交点,跳过即可,影响不大 ^_^
    }

  }

  return a_star_pathes;
}


double BsplineOptimizer::costFunctionRebound(const Eigen::VectorXd& x, Eigen::VectorXd& grad, bool& force_return, void* func_data)
{
  BsplineOptimizer* opt = reinterpret_cast<BsplineOptimizer*>(func_data);

  double cost;
  opt->combineCostRebound(x, grad, cost);

  // early termination
  force_return = ( opt->force_stop_type_ == STOP_FOR_ERROR || opt->force_stop_type_ == STOP_FOR_REBOUND );

  opt->iter_num_ += 1;
  return cost;
}

double BsplineOptimizer::costFunctionRefine(const Eigen::VectorXd& x, Eigen::VectorXd& grad, bool& force_return, void* func_data)
{
  BsplineOptimizer* opt = reinterpret_cast<BsplineOptimizer*>(func_data);

  double cost;
  opt->combineCostRefine(x, grad, cost);

  force_return = false;

  opt->iter_num_ += 1;
  return cost;
}

/**
 * @brief 计算反弹距离代价(弹性带机制的核心代价函数)
 * 
 * 【算法原理】
 * 该代价函数驱动控制点远离障碍物,实现"反弹"效果。
 * 
 * 核心思想:
 *   1. 对每个控制点,计算其沿"弹性方向"到基点的距离 dist
 *   2. 期望距离为 clearance(安全距离),误差为 dist_err = clearance - dist
 *   3. 如果控制点离障碍物足够远(dist_err < 0),无代价
 *   4. 如果控制点离障碍物不够远,施加代价
 * 
 * 【代价函数设计(两段式)】
 *   为防止控制点产生过大的跳跃,使用分段代价函数:
 *   第一段(dist_err < demarcation): 三次函数 J = dist_err³
 *     - 连续且平滑(一阶导和二阶导连续)
 *     - 梯度: dJ/dx = 3*dist_err²
 *   第二段(dist_err >= demarcation): 二次函数 J = a*dist_err² + b*dist_err + c
 *     - 在demarcation处与第一段C1连续
 *     - 避免三次函数在距离过远时梯度增长过快
 * 
 * 梯度 grad = -(dJ/d(dist_err)) * direction
 * 负号表示沿方向向量(指向安全区域)移动可减小代价
 * 
 * @param q 控制点矩阵
 * @param cost 输出的代价
 * @param gradient 输出的梯度矩阵
 * @param iter_num 当前迭代次数
 * @param smoothness_cost 当前的平滑度代价(用于判断轨迹是否够平滑)
 */
void BsplineOptimizer::calcDistanceCostRebound(const Eigen::MatrixXd& q, double& cost,
                                        Eigen::MatrixXd& gradient, int iter_num, double smoothness_cost)
{
  cost = 0.0;
  int end_idx = q.cols() - order_;
  double demarcation = cps_.clearance;
  /* 两段代价函数的系数(确保C1连续性) */
  double a = 3*demarcation, b = -3*pow(demarcation, 2), c = pow(demarcation, 3);

  force_stop_type_ = DONT_STOP;
  /**
   * 当轨迹足够平滑时(每段平滑代价 < 0.1),开始检查碰撞并触发反弹。
   * 这个条件确保在轨迹质量较好时才会考虑碰撞避免,
   * 避免在轨迹还很粗糙时就过度优化碰撞。
   */
  if ( iter_num > 3 && smoothness_cost / ( cps_.size - 2*order_) < 0.1 )
  {
    check_collision_and_rebound();
  }

  /*** 计算每个控制点的距离代价和梯度 ***/
  for ( auto i=order_; i<end_idx; ++i )
  {
    for ( size_t j=0; j<cps_.direction[i].size(); ++j )
    {
      /**
       * dist: 控制点到基点在弹性方向上的投影距离
       * dist_err: 期望安全距离与实际距离的差
       * 正值表示控制点离障碍物太近,需要施加代价
       */
      double dist = (cps_.points.col(i) - cps_.base_point[i][j]).dot(cps_.direction[i][j]);
      double dist_err = cps_.clearance - dist;
      Eigen::Vector3d dist_grad = cps_.direction[i][j];
      
      if (dist_err < 0)
      {
        /* 控制点已经足够安全,无代价 */
      }
      else if ( dist_err < demarcation )
      {
        /* 第一段:三次函数代价,梯度平滑增长 */
        cost += pow(dist_err, 3);
        gradient.col(i) += -3.0*dist_err*dist_err * dist_grad;
      }
      else
      {
        /* 第二段:二次函数代价,避免梯度爆炸 */
        cost += a*dist_err*dist_err + b*dist_err + c;
        gradient.col(i) += -(2.0*a*dist_err+b) * dist_grad;
      }
    }
  }

}

/**
 * @brief 计算轨迹拟合代价(精修阶段追踪参考路径)
 * 
 * 【算法原理】
 * 该代价函数使优化后的轨迹贴近参考路径(ref_pts_)。
 * 
 * 使用"椭球距离"度量:
 *   f = |x·v|²/a² + |x×v|²/b²
 * 
 * 其中:
 *   - x: B样条曲线上的点与控制点之间的偏差
 *   - v: 参考路径的切线方向(归一化)
 *   - a: 沿切线方向的容差(较大,默认25)
 *   - b: 垂直切线方向的容差(较小,默认1)
 * 
 * 这个椭球度量的含义:
 *   - 沿路径方向(切向)允许较大的偏差(有利于平滑)
 *   - 垂直路径方向(法向)要求紧密跟踪(有利于避障)
 *   - 形成"雪茄"形的允许偏差区域
 * 
 * 使用B样条求值公式(均匀三次B样条):
 *   P(u) = (P_{i-1} + 4P_i + P_{i+1}) / 6  当 u 在节点区间中点时
 * 
 * 梯度通过链式法则传播到相邻的3个控制点。
 * 
 * @param q 控制点矩阵
 * @param cost 输出的拟合代价
 * @param gradient 输出的梯度矩阵
 */
void BsplineOptimizer::calcFitnessCost(const Eigen::MatrixXd& q, double& cost, Eigen::MatrixXd& gradient)
{

  cost = 0.0;

  int end_idx = q.cols() - order_;

  /* 椭球度量参数: a-切向容差(sqrt(25)=5m), b-法向容差(sqrt(1)=1m) */
  double a2 = 25, b2 = 1;
  for ( auto i=order_-1; i<end_idx+1; ++i )
  {
    /* B样条上的点(参数u在区间中点时的简化求值公式) */
    Eigen::Vector3d x = (q.col(i-1)+4*q.col(i)+q.col(i+1))/6.0 - ref_pts_[i-1];
    /* 参考路径的切线方向 */
    Eigen::Vector3d v = (ref_pts_[i]-ref_pts_[i-2]).normalized();

    double xdotv = x.dot(v);  // 沿切向分量
    Eigen::Vector3d xcrossv = x.cross(v);  // 垂直切向分量(叉积)

    /* 椭球距离代价 */
    double f = pow((xdotv),2)/a2 + pow(xcrossv.norm(),2)/b2;
    cost += f;

    /* 梯度计算: d(|x×v|²)/dx = 2*[v]×*x×v, 其中[v]×是v的反对称矩阵 */
    Eigen::Matrix3d m;
    m << 0,-v(2),v(1), v(2),0,-v(0), -v(1),v(0),0;
    Eigen::Vector3d df_dx = 2*xdotv/a2*v + 2/b2*m*xcrossv;

    /* 根据B样条求值公式,将梯度分配到三个相邻控制点 */
    gradient.col(i-1) += df_dx/6;
    gradient.col(i) += 4*df_dx/6;
    gradient.col(i+1) += df_dx/6;
  }

}


/**
 * @brief 计算平滑度代价(最小化加加速度或加速度)
 * 
 * 【算法原理】
 * 对于均匀三次B样条,加加速度(jerk)可以通过控制点的三阶差分近似:
 *   jerk_i = P_{i+3} - 3*P_{i+2} + 3*P_{i+1} - P_i
 * 
 * 加加速度的平方和作为平滑度代价:
 *   J_smooth = Σ ||jerk_i||²
 * 
 * 这个代价最小化意味着:
 *   - 控制点序列尽可能线性变化
 *   - 轨迹更加平滑,减少急剧转弯
 *   - 对于四旋翼飞行器,加加速度小意味着姿态变化平缓
 * 
 * 梯度计算(以jerk为例):
 *   d(||jerk_i||²)/dP_i = -2*jerk_i
 *   d(||jerk_i||²)/dP_{i+1} = 6*jerk_i
 *   d(||jerk_i||²)/dP_{i+2} = -6*jerk_i
 *   d(||jerk_i||²)/dP_{i+3} = 2*jerk_i
 *   每个控制点会收到来自最多4个jerk项的梯度贡献
 * 
 * @param q 控制点矩阵
 * @param cost 输出的平滑度代价
 * @param gradient 输出的梯度矩阵
 * @param falg_use_jerk 是否使用加加速度(否则使用加速度)
 */
void BsplineOptimizer::calcSmoothnessCost(const Eigen::MatrixXd& q, double& cost,
                                          Eigen::MatrixXd& gradient, bool falg_use_jerk/* = true*/) {
  
  cost = 0.0;

  if ( falg_use_jerk )
  {
    /*** 使用加加速度(jerk)平滑度: 三阶差分 ***/
    Eigen::Vector3d jerk, temp_j;

    for (int i = 0; i < q.cols() - 3; i++) {
      /* 计算jerk = Δ³P = P_{i+3} - 3P_{i+2} + 3P_{i+1} - P_i */
      jerk = q.col(i + 3) - 3 * q.col(i + 2) + 3 * q.col(i + 1) - q.col(i);
      cost += jerk.squaredNorm();
      temp_j = 2.0 * jerk;
      /* jerk梯度累加到4个相关控制点 */
      gradient.col(i + 0) += -temp_j;
      gradient.col(i + 1) += 3.0 * temp_j;
      gradient.col(i + 2) += -3.0 * temp_j;
      gradient.col(i + 3) += temp_j;
    }
  }
  else
  {    
    /*** 使用加速度(acceleration)平滑度: 二阶差分 ***/
    Eigen::Vector3d acc, temp_acc;

    for (int i = 0; i < q.cols() - 2; i++) {
      /* 计算acc = Δ²P = P_{i+2} - 2P_{i+1} + P_i */
      acc = q.col(i + 2) - 2 * q.col(i + 1) + q.col(i);
      cost += acc.squaredNorm();
      temp_acc = 2.0 * acc;
      /* acc梯度累加到3个相关控制点 */
      gradient.col(i + 0) += temp_acc;
      gradient.col(i + 1) += -2.0 * temp_acc;
      gradient.col(i + 2) += temp_acc;
    }
  }
  
}

/**
 * @brief 计算可行性代价(速度和加速度约束)
 * 
 * 【算法原理】
 * 确保轨迹满足动力学可行性约束:
 *   1. 速度约束: |v| ≤ max_vel (每个轴独立约束)
 *   2. 加速度约束: |a| ≤ max_acc (每个轴独立约束)
 * 
 * 对于均匀B样条,速度和加速度可以通过控制点的差分近似计算:
 *   v_i = (P_{i+1} - P_i) / Δt   (一阶差分)
 *   a_i = (P_{i+2} - 2P_{i+1} + P_i) / Δt²   (二阶差分)
 * 
 * 【分段代价函数】
 * 与距离代价类似,使用两段式惩罚函数来平滑地施加约束:
 *   当超过限制但小于限制+demarcation: 三次函数(平滑过渡)
 *   当超过限制+demarcation: 二次函数(更强惩罚)
 * 
 * @param q 控制点矩阵
 * @param cost 输出的可行性代价
 * @param gradient 输出的梯度矩阵
 */
void BsplineOptimizer::calcFeasibilityCost(const Eigen::MatrixXd& q, double& cost,
                                           Eigen::MatrixXd& gradient) {
  cost = 0.0;
  double demarcation = 1.0; // 1m/s, 1m/s/s
  double ar = 3*demarcation, br = -3*pow(demarcation, 2), cr = pow(demarcation, 3);
  double al = ar, bl = -br, cl = cr;

  /* abbreviation */
  double ts, ts_inv2;
  ts      = bspline_interval_;
  ts_inv2 = 1 / ts / ts;

  /* velocity feasibility */
  for (int i = 0; i < q.cols() - 1; i++) 
  {
    Eigen::Vector3d vi = (q.col(i + 1) - q.col(i))/ts;

    for (int j = 0; j < 3; j++) 
    {
      if ( vi(j) > max_vel_+demarcation )
      {
        double diff = vi(j) - max_vel_;
        cost += ar*diff*diff + br*diff + cr;

        double grad = (2.0*ar*diff+br)/ts;
        gradient(j, i+0) += -grad;
        gradient(j, i+1) += grad;
      }
      else if ( vi(j) > max_vel_ )
      {
        double diff = vi(j) - max_vel_;
        cost += pow( diff, 3 );

        double grad = 3*diff*diff/ts;
        gradient(j, i+0) += -grad;
        gradient(j, i+1) += grad;
      }
      else if ( vi(j) < -(max_vel_+demarcation) )
      {
        double diff = vi(j) + max_vel_;
        cost += al*diff*diff + bl*diff + cl;

        double grad = (2.0*al*diff+bl)/ts;
        gradient(j, i+0) += -grad;
        gradient(j, i+1) += grad;
      }
      else if ( vi(j) < -max_vel_ )
      {
        double diff = vi(j) + max_vel_;
        cost += -pow( diff, 3 );

        double grad = -3*diff*diff/ts;
        gradient(j, i+0) += -grad;
        gradient(j, i+1) += grad;
      }
      else
      {
        /* nothing happened */
      }
    }
  }

  /* acceleration feasibility */
  for (int i = 0; i < q.cols() - 2; i++) 
  {
    Eigen::Vector3d ai = (q.col(i + 2) - 2 * q.col(i + 1) + q.col(i))*ts_inv2;

    for (int j = 0; j < 3; j++) 
    {
      if ( ai(j) > max_acc_+demarcation )
      {
        double diff = ai(j) - max_acc_;
        cost += ar*diff*diff + br*diff + cr;

        double grad = (2.0*ar*diff+br)*ts_inv2;
        gradient(j, i+0) += grad;
        gradient(j, i+1) += -2*grad;
        gradient(j, i+2) += grad;
      }
      else if ( ai(j) > max_acc_ )
      {
        double diff = ai(j) - max_acc_;
        cost += pow( diff, 3 );

        double grad = 3*diff*diff*ts_inv2;
        gradient(j, i+0) += grad;
        gradient(j, i+1) += -2*grad;
        gradient(j, i+2) += grad;
      }
      else if ( ai(j) < -(max_acc_+demarcation) )
      {
        double diff = ai(j) + max_acc_;
        cost += al*diff*diff + bl*diff + cl;

        double grad = (2.0*al*diff+bl)*ts_inv2;
        gradient(j, i+0) += grad;
        gradient(j, i+1) += -2*grad;
        gradient(j, i+2) += grad;
      }
      else if ( ai(j) < -max_acc_ )
      {
        double diff = ai(j) + max_acc_;
        cost += -pow( diff, 3 );

        double grad = -3*diff*diff*ts_inv2;
        gradient(j, i+0) += grad;
        gradient(j, i+1) += -2*grad;
        gradient(j, i+2) += grad;
      }
      else
      {
        /* nothing happened */
      }
    }
  }

}



bool BsplineOptimizer::check_collision_and_rebound(void)
{

  int end_idx = cps_.size - order_;

  /*** Check and segment the initial trajectory according to obstacles ***/
  int in_id, out_id;
  vector<std::pair<int,int>> segment_ids;
  bool flag_new_obs_valid = false;
  int i_end = end_idx - (end_idx-order_)/3;
  for ( int i=order_-1; i<=i_end; ++i )
  {

    bool occ = grid_map_->getInflateOccupancy(cps_.points.col(i));

    /*** check if the new collision will be valid ***/
    if ( occ )
    {
      for ( size_t k=0; k<cps_.direction[i].size(); ++k )
      {
        cout.precision(2);
        if ( ( cps_.points.col(i) - cps_.base_point[i][k] ).dot(cps_.direction[i][k]) < 1 * grid_map_->getResolution() ) // current point is outside all the collision_points. 
        {
          occ = false;  // Not really takes effect, just for better hunman understanding.
          break;
        }
      }
    }

    if ( occ )
    {
      flag_new_obs_valid = true;

      int j;
      for ( j=i-1; j>=0; --j )
      {
        occ = grid_map_->getInflateOccupancy(cps_.points.col(j));
        if ( !occ )
        {
          in_id = j;
          break;
        }
      }
      if ( j < 0 ) // fail to get the obs free point
      {
        ROS_ERROR( "ERROR! the drone is in obstacle. This should not happen." );
        in_id = 0;
      }

      for ( j=i+1; j<cps_.size; ++j )
      {
        occ = grid_map_->getInflateOccupancy(cps_.points.col(j));
        
        if ( !occ )
        {
          out_id = j;
          break;
        }
      }
      if ( j >= cps_.size ) // fail to get the obs free point
      {
        ROS_WARN( "WARN! terminal point of the current trajectory is in obstacle, skip this planning." );

        force_stop_type_ = STOP_FOR_ERROR;
        return false;
      }

      i = j+1;

      segment_ids.push_back( std::pair<int,int>(in_id, out_id) );
    }
  }


  if ( flag_new_obs_valid )
  {
    vector<vector<Eigen::Vector3d>> a_star_pathes;
    for ( size_t i=0; i<segment_ids.size(); ++i )
    {
      /*** a star search ***/
      Eigen::Vector3d in( cps_.points.col(segment_ids[i].first) ), out( cps_.points.col(segment_ids[i].second) );
      if ( a_star_->AstarSearch( /*(in-out).norm()/10+0.05*/0.1, in, out) )
      {
        a_star_pathes.push_back( a_star_->getPath() );
      }
      else
      {
        ROS_ERROR("a star error");
        segment_ids.erase( segment_ids.begin() + i );
        i--;
      }

    }

    /*** Assign parameters to each segment ***/
    for (size_t i=0; i<segment_ids.size(); ++i)
    {
      // step 1
      for ( int j=segment_ids[i].first; j <= segment_ids[i].second; ++j )
        cps_.flag_temp[j] = false;

      // step 2
      int got_intersection_id = -1;
      for ( int j=segment_ids[i].first+1; j<segment_ids[i].second; ++j )
      {
        Eigen::Vector3d ctrl_pts_law(cps_.points.col(j+1) - cps_.points.col(j-1)), intersection_point;
        int Astar_id = a_star_pathes[i].size() / 2, last_Astar_id; // Let "Astar_id = id_of_the_most_far_away_Astar_point" will be better, but it needs more computation
        double val = (a_star_pathes[i][Astar_id] - cps_.points.col(j)).dot( ctrl_pts_law ), last_val = val;
        while ( Astar_id >=0 && Astar_id < (int)a_star_pathes[i].size() )
        {
          last_Astar_id = Astar_id;

          if ( val >= 0 )
            -- Astar_id;
          else
            ++ Astar_id;
          
          val = (a_star_pathes[i][Astar_id] - cps_.points.col(j)).dot( ctrl_pts_law );
          
          if ( val * last_val <= 0 && ( abs(val) > 0 || abs(last_val) > 0 ) ) // val = last_val = 0.0 is not allowed
          {
            intersection_point = 
              a_star_pathes[i][Astar_id] + 
              ( ( a_star_pathes[i][Astar_id] - a_star_pathes[i][last_Astar_id] ) * 
                ( ctrl_pts_law.dot( cps_.points.col(j) - a_star_pathes[i][Astar_id] ) / ctrl_pts_law.dot( a_star_pathes[i][Astar_id] -  a_star_pathes[i][last_Astar_id] ) ) // = t
              );

            got_intersection_id = j;
            break;
          }
        }

        if ( got_intersection_id >= 0 )
        {
          cps_.flag_temp[j] = true;
          double length = (intersection_point - cps_.points.col(j)).norm();
          if ( length > 1e-5 )
          {
            for ( double a=length; a>=0.0; a-=grid_map_->getResolution() )
            {
              bool occ = grid_map_->getInflateOccupancy((a/length)*intersection_point + (1-a/length)*cps_.points.col(j));
        
              if ( occ || a < grid_map_->getResolution() )
              {
                if ( occ )
                  a+=grid_map_->getResolution();
                cps_.base_point[j].push_back( (a/length)*intersection_point + (1-a/length)*cps_.points.col(j) );
                cps_.direction[j].push_back( (intersection_point - cps_.points.col(j)).normalized() );
                break;
              }
            }
          }
          else
          {
            got_intersection_id = -1;
          }
        }
      }

      //step 3
      if ( got_intersection_id >= 0 )
      {
        for ( int j=got_intersection_id + 1; j <= segment_ids[i].second; ++j )
          if ( ! cps_.flag_temp[j] )
          {
            cps_.base_point[j].push_back( cps_.base_point[j-1].back() );
            cps_.direction[j].push_back( cps_.direction[j-1].back() );
          }

        for ( int j=got_intersection_id - 1; j >= segment_ids[i].first; --j )
          if ( ! cps_.flag_temp[j] )
          {
            cps_.base_point[j].push_back( cps_.base_point[j+1].back() );
            cps_.direction[j].push_back( cps_.direction[j+1].back() );
          }
      }
      else
        ROS_WARN("Failed to generate direction!");
    }

    force_stop_type_ = STOP_FOR_REBOUND;
    return true;
  }

  return false;
}

bool BsplineOptimizer::BsplineOptimizeTrajRebound(Eigen::MatrixXd& optimal_points, double ts)
{
  setBsplineInterval(ts);

  bool flag_success =  rebound_optimize();

  optimal_points = cps_.points;

  return flag_success;
}


bool BsplineOptimizer::BsplineOptimizeTrajRefine(const Eigen::MatrixXd& init_points, const double ts, Eigen::MatrixXd& optimal_points)
{

  setControlPoints(init_points);
  setBsplineInterval(ts);
  
  bool flag_success =  refine_optimize();

  optimal_points = cps_.points;

  return flag_success;
}


/**
 * @brief 反弹优化主函数(核心优化循环)
 * 
 * 【算法流程】
 * 1. 初始化梯度下降求解器,设置最大迭代次数和最小梯度
 * 2. 优化变量: 内部控制点(索引order_到size-order_)
 * 3. 迭代优化,每次调用costFunctionRebound计算代价和梯度
 * 4. 优化后检查轨迹是否安全(沿轨迹采样检测碰撞)
 * 5. 如果仍有碰撞:
 *    - 增加碰撞权重(lambda2_ *= 2)
 *    - 重新初始化控制点(initControlPoints)
 *    - 重新开始优化(最多MAX_RESART_NUMS_SET=3次)
 * 6. 如果检测到新的障碍物(通过RETURN_BY_ORDER),触发反弹:
 *    - 再次调用check_collision_and_rebound更新基点和方向
 *    - 继续优化(最多20次反弹)
 * 
 * @return true 优化成功(轨迹无碰撞), false 优化失败
 */
bool BsplineOptimizer::rebound_optimize()
{
  /* ---------- 初始化求解器 ---------- */
  iter_num_ = 0;
  int start_id = order_;
  int end_id = this->cps_.size - order_;
  variable_num_ = 3 * (end_id - start_id);  // 优化变量数为内部控制点的3维坐标

  // cout << "variable_num_=" << variable_num_ << endl;

  GradientDescentOptimizer opt(variable_num_, BsplineOptimizer::costFunctionRebound, this);

  opt.set_maxeval(150);      // 最大函数调用次数
  opt.set_min_grad(1e-2);    // 最小梯度阈值(收敛判定)

  /* ---------- 初始化优化变量 ---------- */
  double final_cost;

  ros::Time t0 = ros::Time::now(), t1, t2;
  int restart_nums = 0, rebound_times = 0;;
  bool flag_force_return, flag_occ, success;
  double original_lambda2 = lambda2_;
  constexpr int MAX_RESART_NUMS_SET = 3;
  do
  {
    /* ---------- 准备 ---------- */
    min_cost_ = std::numeric_limits<double>::max();
    iter_num_ = 0;
    flag_force_return = false;
    flag_occ = false;
    success = false;

    /* 提取内部控制点作为优化变量 */
    Eigen::VectorXd q(variable_num_);
    memcpy( q.data(), cps_.points.data() + 3*order_, variable_num_*sizeof(q(0)) );

    /* ---------- 执行优化 ---------- */
    t1 = ros::Time::now();
    auto result = opt.optimize(q, final_cost);
    t2 = ros::Time::now();
    double time_ms = (t2-t1).toSec()*1000;
    double total_time_ms = (t2-t0).toSec()*1000;

    /* ---------- 优化成功,检查碰撞 ---------- */
    if ( result == GradientDescentOptimizer::FIND_MIN ||  result == GradientDescentOptimizer::REACH_MAX_ITERATION )
    {
      flag_force_return = false;

      /* 沿B样条轨迹采样,检查是否与障碍物碰撞(只检查前2/3段) */
      UniformBspline traj =  UniformBspline(cps_.points, 3, bspline_interval_);
      double tm, tmp;
      traj.getTimeSpan(tm, tmp);
      double t_step = (tmp - tm) / ( (traj.evaluateDeBoorT(tmp) - traj.evaluateDeBoorT(tm)).norm() / grid_map_->getResolution() );
      for ( double t = tm; t<tmp*2/3; t+=t_step )
      {
        flag_occ = grid_map_->getInflateOccupancy( traj.evaluateDeBoorT(t) );
        if ( flag_occ )
        {
          /* 如果前3个控制点就在障碍物中,返回失败 */
          if ( t <= bspline_interval_ )
          {
            cout << cps_.points.col(1).transpose() << "\n"  << cps_.points.col(2).transpose() << "\n"  << cps_.points.col(3).transpose() << "\n" << cps_.points.col(4).transpose() << endl;
            ROS_ERROR("First 3 control points in obstacles! return false, t=%f",t);
            return false;
          }

          break;
        }
      }

      if ( !flag_occ )
      {
        /* 优化成功且无碰撞 */
        printf("\033[32miter(+1)=%d,time(ms)=%5.3f,total_t(ms)=%5.3f,cost=%5.3f\n\033[0m", iter_num_, time_ms, total_time_ms, final_cost);
        success = true;
      }
      else
      {
        /* 仍有碰撞:增大碰撞权重,重新初始化,重新优化 */
        restart_nums++;
        initControlPoints(cps_.points, false);
        lambda2_ *= 2;

        printf("\033[32miter(+1)=%d,time(ms)=%5.3f,keep optimizing\n\033[0m", iter_num_, time_ms);
      }

    }
    else if ( result == GradientDescentOptimizer::RETURN_BY_ORDER )
    {
      /* 检测到新的障碍物,触发反弹机制 */
      flag_force_return = true;
      rebound_times ++;
      cout << "iter=" << iter_num_ << ",time(ms)=" << time_ms << ",rebound." << endl;
    }

  } while ( (flag_occ && restart_nums < MAX_RESART_NUMS_SET) || 
            (flag_force_return && force_stop_type_ == STOP_FOR_REBOUND && rebound_times <= 20) 
          );

  lambda2_ = original_lambda2;

  return success;
}


/**
 * @brief 精修优化(在最终优化阶段使用)
 * 
 * 【算法流程】
 * 与反弹优化不同,精修优化更注重轨迹的平滑度和对参考路径的拟合:
 * 1. 使用拟合代价(calcFitnessCost)替代距离代价(calcDistanceCostRebound)
 * 2. 不再检测新的碰撞和触发反弹
 * 3. 如果优化后检测到碰撞,增大拟合权重(lambda4_ *= 2)但只尝试一次
 * 
 * 适用于已经通过反弹优化获得无碰撞轨迹后的精修阶段。
 * 
 * @param init_points 初始控制点
 * @param ts B样条时间间隔
 * @param optimal_points 输出的最优控制点
 * @return true 优化成功, false 优化失败
 */
bool BsplineOptimizer::refine_optimize()
{
  /* ---------- 初始化求解器 ---------- */
  iter_num_ = 0;
  int start_id = order_;
  int end_id = this->cps_.points.cols() - order_;
  variable_num_ = 3 * (end_id - start_id);

  // cout << "variable_num_" << variable_num_ << endl;

  GradientDescentOptimizer opt(variable_num_, BsplineOptimizer::costFunctionRefine, this);

  opt.set_maxiter(100);      // 最大迭代次数
  opt.set_maxeval(100);      // 最大函数调用次数
  opt.set_min_grad(1e-2);    // 最小梯度阈值

  /* ---------- 初始化优化变量 ---------- */
  Eigen::VectorXd q(variable_num_);
  double final_cost;
  for (size_t i = start_id; i < end_id; ++i)
  {
    q.segment( 3 * (i - start_id), 3) = cps_.points.col(i);
  }

  double origin_lambda4 = lambda4_;
  bool flag_safe = true;
  int iter_count = 0;
  do
  {
    auto result = opt.optimize(q, final_cost);

    /* 检查轨迹是否安全 */
    UniformBspline traj =  UniformBspline(cps_.points, 3, bspline_interval_);
    double tm, tmp;
    traj.getTimeSpan(tm, tmp);
    double t_step = (tmp - tm) / ( (traj.evaluateDeBoorT(tmp) - traj.evaluateDeBoorT(tm)).norm() / grid_map_->getResolution() );
    for ( double t = tm; t<tmp*2/3; t+=t_step )
    {
      if ( grid_map_->getInflateOccupancy( traj.evaluateDeBoorT(t) ) )
      {
        flag_safe = false;
        break;
      }
    }

    if ( !flag_safe ) lambda4_*=2;  // 有碰撞时增大拟合权重
    
    iter_count++;
  } while ( !flag_safe && iter_count <= 0 );
  
  lambda4_ = origin_lambda4;

  return flag_safe;
}

/**
 * @brief 组合反弹代价函数(将各子代价加权求和)
 * 
 * 计算平滑度、距离(反弹)、可行性三种代价及其梯度,加权组合。
 * 这是梯度下降优化器调用的目标函数。
 * 
 * 组合代价: f = λ1*f_smooth + λ2*f_distance + λ3*f_feasibility
 * 组合梯度: grad = λ1*g_smooth + λ2*g_distance + λ3*g_feasibility
 * 
 * 注意:只返回内部控制点(索引order_到size-order_)的梯度,
 * 边界控制点固定不优化。
 * 
 * @param x 优化变量(内部控制点的3维坐标向量化)
 * @param grad 输出的梯度
 * @param f_combine 输出的组合代价值
 */
void BsplineOptimizer::combineCostRebound(const Eigen::VectorXd& x, Eigen::VectorXd& grad, double& f_combine)
{
  /* 将优化变量写回控制点矩阵 */
  memcpy( cps_.points.data() + 3*order_, x.data(), x.size()*sizeof(x(0)) );

  /* ---------- 计算各子代价和梯度 ---------- */
  double f_smoothness, f_distance, f_feasibility;

  Eigen::MatrixXd g_smoothness = Eigen::MatrixXd::Zero(3, cps_.size);
  Eigen::MatrixXd g_distance = Eigen::MatrixXd::Zero(3, cps_.size);
  Eigen::MatrixXd g_feasibility = Eigen::MatrixXd::Zero(3, cps_.size);

  calcSmoothnessCost(cps_.points, f_smoothness, g_smoothness);
  calcDistanceCostRebound(cps_.points, f_distance, g_distance, iter_num_, f_smoothness);
  calcFeasibilityCost(cps_.points, f_feasibility, g_feasibility);

  /* 加权组合代价 */
  f_combine = lambda1_ * f_smoothness + lambda2_ * f_distance + lambda3_ * f_feasibility;

  /* 加权组合梯度,只提取内部控制点的梯度 */
  Eigen::MatrixXd grad_3D = lambda1_ * g_smoothness + lambda2_ * g_distance + lambda3_ * g_feasibility;
  grad.resize(variable_num_);
  memcpy( grad.data(), grad_3D.data() + 3*order_, grad.size()*sizeof(grad(0)) );
}

/**
 * @brief 组合精修代价函数
 * 
 * 与combineCostRebound不同,精修阶段使用拟合代价替代距离代价:
 *   f = λ1*f_smooth + λ4*f_fitness + λ3*f_feasibility
 * 
 * 拟合代价使轨迹更贴近参考路径(通常是无碰撞的路径),
 * 适用于最终的精修阶段。
 * 
 * @param x 优化变量
 * @param grad 输出的梯度
 * @param f_combine 输出的组合代价值
 */
void BsplineOptimizer::combineCostRefine(const Eigen::VectorXd& x, Eigen::VectorXd& grad, double& f_combine)
{
  /* 将优化变量写回控制点矩阵 */
  memcpy( cps_.points.data() + 3*order_, x.data(), x.size()*sizeof(x(0)) );

  /* ---------- 计算各子代价和梯度 ---------- */
  double f_smoothness, f_fitness, f_feasibility;

  Eigen::MatrixXd g_smoothness = Eigen::MatrixXd::Zero(3, cps_.points.cols());
  Eigen::MatrixXd g_fitness = Eigen::MatrixXd::Zero(3, cps_.points.cols());
  Eigen::MatrixXd g_feasibility = Eigen::MatrixXd::Zero(3, cps_.points.cols());

  calcSmoothnessCost(cps_.points, f_smoothness, g_smoothness);
  calcFitnessCost(cps_.points, f_fitness, g_fitness);
  calcFeasibilityCost(cps_.points, f_feasibility, g_feasibility);

  /* ---------- 转换为NLopt格式... ---------- */
  f_combine = lambda1_ * f_smoothness + lambda4_ * f_fitness + lambda3_ * f_feasibility;

  Eigen::MatrixXd grad_3D = lambda1_ * g_smoothness + lambda4_ * g_fitness + lambda3_ * g_feasibility;
  grad.resize(variable_num_);
  memcpy( grad.data(), grad_3D.data() + 3*order_, grad.size()*sizeof(grad(0)) );
}


}  // namespace rebound_planner