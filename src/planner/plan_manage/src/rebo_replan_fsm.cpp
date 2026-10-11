
#include <plan_manage/rebo_replan_fsm.h>

namespace rebound_planner {

/**
 * @brief FSM初始化: 配置参数、初始化模块、注册回调
 * 
 * 【整体架构】
 * ReboReplanFSM是Ego算法的核心调度模块,采用有限状态机(FSM)架构
 * 管理轨迹规划和重规划的全生命周期:
 * 
 *   INIT → WAIT_TARGET → GEN_NEW_TRAJ → EXEC_TRAJ → (循环)
 *                          ↑              ↓
 *                     REPLAN_TRAJ ← EMERGENCY_STOP
 * 
 * 【各模块职责】
 * - planner_manager_: 实际的轨迹规划引擎(A*搜索+B样条优化)
 * - visualization_: 可视化工具(Rviz显示)
 * - 定时器: exec_timer_(10ms)驱动FSM状态机, safety_timer_(50ms)检测碰撞
 * 
 * @param nh ROS节点句柄,用于读取参数和创建通信接口
 */
void ReboReplanFSM::init(ros::NodeHandle& nh) {
  current_wp_  = 0;
  exec_state_  = FSM_EXEC_STATE::INIT;   // 初始状态: 等待触发
  have_target_ = false;
  have_odom_   = false;

  /* ---- FSM 配置参数 ---- */
  nh.param("fsm/flight_type", target_type_, -1);          // 飞行模式(手动/预设航点)
  nh.param("fsm/thresh_replan", replan_thresh_, -1.0);    // 触发重规划的位移阈值
  nh.param("fsm/thresh_no_replan", no_replan_thresh_, -1.0); // 接近目标时不重规划的距离阈值
  nh.param("fsm/planning_horizen", planning_horizen_, -1.0);     // 局部规划的水平距离
  nh.param("fsm/planning_horizen_time", planning_horizen_time_, -1.0); // 局部规划的水平时间

  // 预设航点(用于PRESET_TARGET模式)
  nh.param("fsm/waypoint_num", waypoint_num_, -1);
  for (int i = 0; i < waypoint_num_; i++) {
    nh.param("fsm/waypoint" + to_string(i) + "_x", waypoints_[i][0], -1.0);
    nh.param("fsm/waypoint" + to_string(i) + "_y", waypoints_[i][1], -1.0);
    nh.param("fsm/waypoint" + to_string(i) + "_z", waypoints_[i][2], -1.0);
  }

  /* ---- 初始化核心模块 ---- */
  visualization_.reset(new PlanningVisualization(nh));
  planner_manager_.reset(new ReboundPlannerManager);
  planner_manager_->initPlanModules(nh, visualization_);

  /* ---- 注册ROS定时器和通信接口 ---- */
  // 主FSM驱动定时器: 10ms执行一次状态检查和转移
  exec_timer_   = nh.createTimer(ros::Duration(0.01), &ReboReplanFSM::execFSMCallback, this);
  // 安全碰撞检测定时器: 50ms检查一次当前轨迹是否与障碍物碰撞
  safety_timer_ = nh.createTimer(ros::Duration(0.05), &ReboReplanFSM::checkCollisionCallback, this);

  // 航点订阅(来自路径生成器/用户点击)
  waypoint_sub_ = nh.subscribe("/waypoint_generator/waypoints", 1, &ReboReplanFSM::waypointCallback, this);
  // 里程计订阅(无人机位置/速度/姿态)
  odom_sub_ = nh.subscribe("/odom_world", 1, &ReboReplanFSM::odometryCallback, this);

  // 发布B样条轨迹给执行模块
  bspline_pub_ = nh.advertise<rebound_planner::Bspline>("/planning/bspline", 10);
  // 发布数据显示(用于调试)
  data_disp_pub_ = nh.advertise<rebound_planner::DataDisp>("/planning/data_display", 100);
}

/**
 * @brief 航点回调: 接收用户/路径生成器发送的目标点
 * 
 * 【处理流程】
 * 1. 根据飞行模式(MANUAL_TARGET/PRESET_TARGET)设置终点
 * 2. 调用全局规划器生成从当前位置到终点的全局轨迹
 * 3. 更新FSM状态,触发局部轨迹重规划
 * 
 * 【两种目标模式】
 * - MANUAL_TARGET: 用户在Rviz中点击的目标点(x,y,固定z=1.0)
 * - PRESET_TARGET: 从配置文件读取的预设航点序列(循环遍历)
 * 
 * @param msg 路径消息,包含目标点位置
 */
void ReboReplanFSM::waypointCallback(const nav_msgs::PathConstPtr& msg) {
  if (msg->poses[0].pose.position.z < -0.1) return;

  cout << "Triggered!" << endl;
  trigger_ = true;

  // ---- 根据飞行模式设置目标点 ----
  if (target_type_ == TARGET_TYPE::MANUAL_TARGET) {
    // 手动点击模式: 固定飞行高度1.0米
    end_pt_ << msg->poses[0].pose.position.x, msg->poses[0].pose.position.y, 1.0;

  } else if (target_type_ == TARGET_TYPE::PRESET_TARGET) {
    // 预设航点模式: 从配置文件读取航点序列
    end_pt_(0)  = waypoints_[current_wp_][0];
    end_pt_(1)  = waypoints_[current_wp_][1];
    end_pt_(2)  = waypoints_[current_wp_][2];
    current_wp_ = (current_wp_ + 1) % waypoint_num_;  // 循环到下一个航点
  }

  init_pt_ = odom_pos_;  // 记录初始位置(用于计算局部目标)

  // ---- 规划全局轨迹(从当前位置到目标点) ----
  // 使用多项式轨迹生成全局参考路径(非B样条,用于引导局部规划)
  bool success = planner_manager_->planGlobalTraj( odom_pos_, odom_vel_, Eigen::Vector3d::Zero(), end_pt_, Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero() );
  if ( success )
  {

    /*** 可视化全局轨迹 ***/
    vector<Eigen::Vector3d> gloabl_traj;
    for ( double t=0; t<planner_manager_->global_data_.global_duration_; t+=1 )
      gloabl_traj.push_back( planner_manager_->global_data_.global_traj_.evaluate(t) );

    end_vel_.setZero();
    have_target_ = true;        // 标记已接收到目标
    have_new_target_ = true;    // 标记为新目标(用于局部规划初始化)

    /*** 根据当前状态触发FSM转移 ***/
    if (exec_state_ == WAIT_TARGET)
      changeFSMExecState(GEN_NEW_TRAJ, "TRIG");  // 等待目标→生成新轨迹
    else if (exec_state_ == EXEC_TRAJ)
      changeFSMExecState(REPLAN_TRAJ, "TRIG");   // 执行中→重规划(新目标覆盖)

    // 可视化: 显示目标点和全局路径
    visualization_->displayGoalPoint(end_pt_, Eigen::Vector4d(1, 0, 0, 1), 0.3, 0);
    visualization_->displayGlobalPathList(gloabl_traj, 0.1, 0);
  }
  else
  {
    ROS_ERROR( "Unable to generate global trajectory!" );
  }
  
}

/**
 * @brief 里程计回调: 更新无人机状态估计
 * 
 * 从里程计消息中提取无人机的位置、速度、姿态信息。
 * 这些信息是FSM和规划器运行的基础:
 * - 位置: 用于确定无人机是否到达目标/偏离轨迹
 * - 速度: 用于轨迹优化的边界条件(连续性约束)
 * - 姿态: 用于偏航角控制(当前未使用,代码已注释)
 * 
 * @param msg 里程计消息(VIO/视觉惯导里程计或Mocap)
 */
void ReboReplanFSM::odometryCallback(const nav_msgs::OdometryConstPtr& msg) {
  // 提取位置信息(x, y, z)
  odom_pos_(0) = msg->pose.pose.position.x;
  odom_pos_(1) = msg->pose.pose.position.y;
  odom_pos_(2) = msg->pose.pose.position.z;

  // 提取速度信息(vx, vy, vz)
  odom_vel_(0) = msg->twist.twist.linear.x;
  odom_vel_(1) = msg->twist.twist.linear.y;
  odom_vel_(2) = msg->twist.twist.linear.z;

  //odom_acc_ = estimateAcc( msg );  // 加速度估计(当前未使用)

  // 提取姿态四元数
  odom_orient_.w() = msg->pose.pose.orientation.w;
  odom_orient_.x() = msg->pose.pose.orientation.x;
  odom_orient_.y() = msg->pose.pose.orientation.y;
  odom_orient_.z() = msg->pose.pose.orientation.z;

  have_odom_ = true;  // 标记已收到里程计数据
}

/**
 * @brief 执行FSM状态转移
 * 
 * 【功能】
 * 改变当前FSM状态,并记录状态连续调用次数。
 * 连续调用次数用于判断重规划是否陷入死循环:
 * - 如果连续处于同一个状态多次,说明可能需要更激进的策略
 * - 具体用在GEN_NEW_TRAJ状态: 首次调用使用多项式初始化,
 *   失败后启用随机初始化尝试跳出局部最优
 * 
 * @param new_state 目标状态
 * @param pos_call 调用来源标识(用于日志: "FSM"/"TRIG"/"SAFETY")
 */
void ReboReplanFSM::changeFSMExecState(FSM_EXEC_STATE new_state, string pos_call) {
  
  // 如果目标状态与当前状态相同,连续调用计数+1
  // 否则重置计数(表示状态发生了有效转移)
  if ( new_state == exec_state_ )
    continously_called_times_ ++;
  else
    continously_called_times_ = 1;

  static string state_str[7] = { "INIT", "WAIT_TARGET", "GEN_NEW_TRAJ", "REPLAN_TRAJ", "EXEC_TRAJ", "EMERGENCY_STOP" };
  int    pre_s        = int(exec_state_);
  exec_state_         = new_state;
  cout << "[" + pos_call + "]: from " + state_str[pre_s] + " to " + state_str[int(new_state)] << endl;
  
}

/**
 * @brief 获取当前状态的连续调用次数和状态值
 * 
 * 用于判断规划失败的场景:
 * 如果连续多次处于GEN_NEW_TRAJ或REPLAN_TRAJ状态,
 * 说明当前的规划策略可能失效,需要切换策略(如启用随机初始化)。
 * 
 * @return pair<连续调用次数, 当前FSM状态>
 */
std::pair<int, ReboReplanFSM::FSM_EXEC_STATE> ReboReplanFSM::timesOfConsecutiveStateCalls()
{
  return std::pair<int, FSM_EXEC_STATE>(continously_called_times_, exec_state_);
}

/**
 * @brief 打印当前FSM状态(调试用)
 * 
 * 每100次FSM回调打印一次,同时检查关键前置条件
 */
void ReboReplanFSM::printFSMExecState() {
  static string state_str[7] = { "INIT", "WAIT_TARGET", "GEN_NEW_TRAJ", "REPLAN_TRAJ", "EXEC_TRAJ", "EMERGENCY_STOP" };

  cout << "[FSM]: state: " + state_str[int(exec_state_)] << endl;
}

/**
 * @brief FSM主回调函数(10ms周期执行): Ego算法的核心调度逻辑
 * 
 * 【有限状态机架构详解】
 * 
 * 这是整个Ego算法的"大脑",以10ms的频率运行,根据当前状态执行相应操作。
 * 状态转移关系如下:
 * 
 * ┌─────────────────────────────────────────────────────┐
 * │  INIT ──(有里程计+触发)──→ WAIT_TARGET              │
 * │                           ↓ (有目标)                 │
 * │                       GEN_NEW_TRAJ ←──┐             │
 * │                           ↓           │             │
 * │                       EXEC_TRAJ ──(偏离轨迹)─→ REPLAN_TRAJ
 * │                         ↓   ↓                      │
 * │                    (到达目标) (碰撞)                 │
 * │                         ↓   ↓                      │
 * │                    WAIT_TARGET EMERGENCY_STOP       │
 * │                                         ↓          │
 * │                                    GEN_NEW_TRAJ    │
 * └─────────────────────────────────────────────────────┘
 * 
 * 【状态说明】
 * - INIT: 初始化状态,等待里程计数据
 * - WAIT_TARGET: 等待用户发送目标航点
 * - GEN_NEW_TRAJ: 生成全新的局部轨迹(使用全局路径引导)
 * - REPLAN_TRAJ: 在现有轨迹基础上重规划(保持连续性)
 * - EXEC_TRAJ: 执行轨迹,持续监控是否需要重规划
 * - EMERGENCY_STOP: 紧急停止(突然发现前方障碍物)
 * 
 * @param e ROS定时器事件
 */
void ReboReplanFSM::execFSMCallback(const ros::TimerEvent& e) {

  static int fsm_num = 0;
  fsm_num++;
  // 每100次(1秒)打印一次当前状态和前置条件状态
  if (fsm_num == 100) {
    printFSMExecState();
    if (!have_odom_) cout << "no odom." << endl;
    if (!trigger_) cout << "wait for goal." << endl;
    fsm_num = 0;
  }

  switch (exec_state_) {
    // =========================================================
    // 状态1: INIT - 初始化等待
    // 条件: 有里程计数据 && 收到用户触发指令
    // 转移: → WAIT_TARGET
    // =========================================================
    case INIT: {
      if (!have_odom_) {
        return;  // 没有里程计数据,继续等待
      }
      if (!trigger_) {
        return;  // 没有收到目标触发,继续等待
      }
      changeFSMExecState(WAIT_TARGET, "FSM");
      break;
    }

    // =========================================================
    // 状态2: WAIT_TARGET - 等待目标航点
    // 条件: waypointCallback已设置目标
    // 转移: → GEN_NEW_TRAJ
    // =========================================================
    case WAIT_TARGET: {
      if (!have_target_)
        return;  // 没有目标,继续等待
      else {
        changeFSMExecState(GEN_NEW_TRAJ, "FSM");  // 有目标→生成本地轨迹
      }
      break;
    }

    // =========================================================
    // 状态3: GEN_NEW_TRAJ - 生成全新局部轨迹
    // 
    // 【关键设计: 多级回退策略】
    // 首次进入: 使用多项式初始化(flag_random_poly_init=false)
    // 失败再进: 使用随机初始化(flag_random_poly_init=true)
    //         随机初始化通过随机抖动初始控制点跳出局部最优
    // 
    // 转移: 成功→EXEC_TRAJ; 失败→GEN_NEW_TRAJ(重试)
    // =========================================================
    case GEN_NEW_TRAJ: {
      // 设置初始状态(来自里程计)
      start_pt_  = odom_pos_;
      start_vel_ = odom_vel_;
      start_acc_.setZero();  // 假设初始加速度为0

      // 判断是否启用随机初始化
      // 如果是第一次进入本状态,使用确定性初始化
      // 如果连续多次失败,启用随机初始化尝试跳出局部最优
      bool flag_random_poly_init;
      if ( timesOfConsecutiveStateCalls().first == 1 ) 
        flag_random_poly_init = false;
      else
        flag_random_poly_init = true;

      bool success = callReboundReplan(true, flag_random_poly_init);
      if (success) {
        changeFSMExecState(EXEC_TRAJ, "FSM");  // 规划成功→执行轨迹
        flag_escape_emergency_ = true;  // 重置紧急停止标志
      } else {
        changeFSMExecState(GEN_NEW_TRAJ, "FSM");  // 规划失败→重试
      }
      break;
    }

    // =========================================================
    // 状态4: REPLAN_TRAJ - 重规划(基于当前轨迹)
    // 
    // 与GEN_NEW_TRAJ的区别:
    // - 从当前执行轨迹的末端(而非当前位置)开始规划
    // - 保持与当前轨迹的连续性(C2连续)
    // - 包含三级回退策略(见planFromCurrentTraj)
    // 
    // 转移: 成功→EXEC_TRAJ; 失败→REPLAN_TRAJ(重试)
    // =========================================================
    case REPLAN_TRAJ: {

      if (planFromCurrentTraj() )
      {
        changeFSMExecState(EXEC_TRAJ, "FSM");
      }
      else
      {
        changeFSMExecState(REPLAN_TRAJ, "FSM");
      }
      
      break;
    }

    // =========================================================
    // 状态5: EXEC_TRAJ - 执行轨迹
    // 
    // 【核心功能: 判断是否需要重规划】
    // 在执行过程中持续检查以下条件:
    // 
    // 1. 到达目标: t_cur > duration → 轨迹执行完毕
    //    → 转移到WAIT_TARGET等待下一个目标
    //    
    // 2. 接近目标: 当前位置到目标的距离 < no_replan_thresh
    //    → 继续执行,不触发重规划(避免终点附近不必要调整)
    //    
    // 3. 接近起点: 距轨迹起点的距离 < replan_thresh
    //    → 继续执行,不触发重规划(刚起飞,轨迹还新鲜)
    //    
    // 4. 需要重规划: 既不接近起点也不接近目标
    //    → 转移到REPLAN_TRAJ
    // 
    // 重规划触发阈值的设计意图:
    // - replan_thresh: 无人机飞行了足够的距离后,轨迹可能
    //   因环境变化或偏离而需要更新
    // - no_replan_thresh: 终点附近微小调整的代价可能超过收益
    // =========================================================
    case EXEC_TRAJ: {
      /* 判断是否需要重规划 */
      LocalTrajData* info     = &planner_manager_->local_data_;
      ros::Time      time_now = ros::Time::now();
      double         t_cur    = (time_now - info->start_time_).toSec();
      t_cur                   = min(info->duration_, t_cur);  // 钳制在有效范围内

      Eigen::Vector3d pos = info->position_traj_.evaluateDeBoorT(t_cur);

      // 条件1: 轨迹执行完毕(到达目标)
      if (t_cur > info->duration_ - 1e-2) {
        have_target_ = false;
        changeFSMExecState(WAIT_TARGET, "FSM");
        return;

      } 
      // 条件2: 已接近目标,不需要重规划
      else if ((end_pt_ - pos).norm() < no_replan_thresh_) {
        return;

      } 
      // 条件3: 刚离开起点不久,轨迹还新鲜
      else if ((info->start_pos_ - pos).norm() < replan_thresh_) {
        return;

      } 
      // 条件4: 条件2和3都不满足,需要进行重规划
      else {
        changeFSMExecState(REPLAN_TRAJ, "FSM");
      }
      break;
    }

    // =========================================================
    // 状态6: EMERGENCY_STOP - 紧急停止
    // 
    // 触发条件: 碰撞检测发现0.8秒内将发生碰撞
    // 
    // 行为:
    // 1. 第一次进入: 调用EmergencyStop生成急停轨迹(刹车)
    // 2. 后续进入: 等待无人机速度降到接近0
    // 3. 停止后: 转移→GEN_NEW_TRAJ重新规划
    // =========================================================
    case EMERGENCY_STOP: {

      if ( flag_escape_emergency_ ) // 防止重复调用急停(只调用一次)
      {
        callEmergencyStop( odom_pos_ );
      }
      else
      {
        // 等待无人机速度降到接近0后再重新规划
        if ( odom_vel_.norm() < 0.1 )
          changeFSMExecState(GEN_NEW_TRAJ, "FSM");
      }
      
      flag_escape_emergency_ = false;
      break;
    }
  }

  // 发布调试数据
  data_disp_.header.stamp = ros::Time::now();
  data_disp_pub_.publish(data_disp_);

}

/**
 * @brief 基于当前执行轨迹进行重规划(三级回退策略)
 * 
 * 【重规划流程】
 * 从当前轨迹的执行位置(而非当前位置)获取起始状态,
 * 确保新轨迹与当前执行轨迹在位置、速度、加速度上C2连续。
 * 
 * 【三级回退策略】
 * 第一级: 不使用多项式初始化,不使用随机初始化
 *         (直接从轨迹末端开始B样条优化)
 * 第二级: 使用多项式初始化,不使用随机初始化
 *         (重新生成B样条初始控制点)
 * 第三级: 使用多项式初始化 AND 随机初始化
 *         (随机扰动初始控制点,尝试跳出局部最优)
 * 
 * 这种逐级递进的策略在效率和鲁棒性之间取得平衡:
 * - 大多数情况下第一级就能成功(快速)
 * - 环境复杂时逐级加"剂量"(鲁棒)
 * 
 * @return true 重规划成功; false 全部三级均失败
 */
bool ReboReplanFSM::planFromCurrentTraj()
{

  LocalTrajData* info     = &planner_manager_->local_data_;
  ros::Time      time_now = ros::Time::now();
  double         t_cur    = (time_now - info->start_time_).toSec();

  // 从当前执行轨迹上获取起始状态(保证C2连续)
  start_pt_  = info->position_traj_.evaluateDeBoorT(t_cur);
  start_vel_ = info->velocity_traj_.evaluateDeBoorT(t_cur);
  start_acc_ = info->acceleration_traj_.evaluateDeBoorT(t_cur);


  // 第一级: 不使用多项式/随机初始化
  bool success = callReboundReplan(false, false);

  if ( !success) 
  {
    // 第二级: 使用多项式初始化
    success = callReboundReplan(true, false);
    if ( !success )
    {
      // 第三级: 多项式+随机初始化(最后一搏)
      success = callReboundReplan(true, true);
      if ( !success )
      {
        return false;  // 全部失败
      }
    }
  } 

  return true;
}

/**
 * @brief 碰撞检测回调(50ms周期): 安全监控
 * 
 * 【检测策略: 前2/3轨迹法】
 * 为了提高效率,对轨迹进行检查时采用以下策略:
 * - 只检查当前时间到轨迹2/3长度的部分
 * - 如果当前时间已经超过2/3点,则检查剩余全部轨迹
 * 
 * 这么做的原因是:
 * - 轨迹的后1/3可能会被后续的重规划覆盖
 * - 减少不必要的计算开销
 * 
 * 【碰撞响应策略】
 * 当检测到轨迹与障碍物碰撞时,按紧迫程度分级响应:
 * 
 * 1. 时间充裕(t-t_cur >= 0.8s): 尝试重规划
 *    - 若重规划成功→继续执行新轨迹
 *    - 若重规划失败→触发重规划状态
 * 
 * 2. 时间紧迫(t-t_cur < 0.8s): 紧急刹车
 *    - 立即生成急停轨迹(刹车到静止)
 *    - 进入EMERGENCY_STOP状态
 * 
 * @param e ROS定时器事件
 */
void ReboReplanFSM::checkCollisionCallback(const ros::TimerEvent& e) {
  LocalTrajData* info = &planner_manager_->local_data_;
  auto map = planner_manager_->grid_map_;

  // 等待目标时不检测(没有轨迹在执行)
  if ( exec_state_ == WAIT_TARGET )
    return;

  /* ---- 轨迹碰撞检测 ---- */
  constexpr double time_step = 0.01;  // 10ms采样步长
  double t_cur = ( ros::Time::now() - info->start_time_ ).toSec();
  double t_2_3 = info->duration_*2/3;  // 轨迹的2/3时间点
  for ( double t=t_cur; t<info->duration_; t+=time_step ) 
  {
    // 如果当前时间还没到2/3点,只检测到2/3处
    if ( t_cur < t_2_3 && t >= t_2_3 )
      break;

    // 检查采样点是否在障碍物(膨胀后)内
    if ( map->getInflateOccupancy( info->position_traj_.evaluateDeBoorT(t) ) )
    {
      // ---- 碰撞发生,根据时间紧迫程度响应 ----
      if ( planFromCurrentTraj() )  // 先尝试重规划
      {
        changeFSMExecState(EXEC_TRAJ, "SAFETY");
        return;
      }
      else
      {
        if ( t-t_cur < 0.8 ) // 距离碰撞时间<0.8秒→紧急刹车
        {
          ROS_ERROR("Suddenly discovered obstacles. emergency stop! time=%f",t-t_cur);
          changeFSMExecState(EMERGENCY_STOP, "SAFETY");
        }
        else // 时间还够→重规划
        {
          ROS_WARN("current traj in collision, replan.");
          changeFSMExecState(REPLAN_TRAJ, "SAFETY");
        }
        return;
      }
      break;
    } 
  }
  
}

/**
 * @brief 调用反弹重规划(Ego算法的核心接口)
 * 
 * 【功能说明】
 * 这是Ego算法中B样条轨迹优化的顶层调用接口,串联了:
 * 1. getLocalTarget(): 计算局部规划的目标点和目标速度
 * 2. reboundReplan(): 执行反弹重规划(B样条优化)
 * 3. 结果发布: 将优化后的B样条轨迹发布给执行模块
 * 
 * 【参数说明】
 * flag_use_poly_init: 是否使用多项式轨迹作为B样条初始解
 *   - true: 从多项式轨迹参数化得到初始B样条(更精确)
 *   - false: 从A*路径参数化得到初始B样条(更快)
 * 
 * flag_randomPolyTraj: 是否对初始控制点添加随机扰动
 *   - true: 添加随机扰动(帮助跳出局部最优)
 *   - false: 使用确定性初始化
 * 
 * @param flag_use_poly_init 是否使用多项式初始化
 * @param flag_randomPolyTraj 是否使用随机扰动
 * @return true 规划成功; false 规划失败
 */
bool ReboReplanFSM::callReboundReplan(bool flag_use_poly_init, bool flag_randomPolyTraj) {
  
  // 第1步: 计算局部目标(在全局路径上截取一段)
  getLocalTarget();

  // 第2步: 执行反弹重规划(B样条优化)
  bool plan_success =
      planner_manager_->reboundReplan(start_pt_, start_vel_, start_acc_, local_target_pt_, local_target_vel_, (have_new_target_ || flag_use_poly_init), flag_randomPolyTraj);
  have_new_target_ = false;

  cout << "final_plan_success=" << plan_success << endl;

  // 第3步: 发布优化结果
  if (plan_success) {

    auto info = &planner_manager_->local_data_;

    // 封装B样条轨迹消息
    rebound_planner::Bspline bspline;
    bspline.order      = 3;  // 三次B样条
    bspline.start_time = info->start_time_;
    bspline.traj_id    = info->traj_id_;

    // 填充控制点
    Eigen::MatrixXd pos_pts = info->position_traj_.getControlPoint();
    bspline.pos_pts.reserve( pos_pts.cols() );
    for (int i = 0; i < pos_pts.cols(); ++i) {
      geometry_msgs::Point pt;
      pt.x = pos_pts(0, i);
      pt.y = pos_pts(1, i);
      pt.z = pos_pts(2, i);
      bspline.pos_pts.push_back(pt);
    }

    // 填充节点向量
    Eigen::VectorXd knots = info->position_traj_.getKnot();
    bspline.knots.reserve(knots.rows()); 
    for (int i = 0; i < knots.rows(); ++i) {
      bspline.knots.push_back(knots(i));
    }

    bspline_pub_.publish(bspline);  // 发布轨迹

    visualization_->displayOptimalList( info->position_traj_.get_control_points(), 0 );

  }

  return plan_success;
}

/**
 * @brief 紧急停止: 生成刹车轨迹
 * 
 * 在检测到即将碰撞且时间紧迫时调用。
 * 生成一条从当前位置减速到静止的短轨迹。
 * 
 * @param stop_pos 当前无人机位置(刹车起点)
 * @return true 紧急停止轨迹生成成功
 */
bool ReboReplanFSM::callEmergencyStop( Eigen::Vector3d stop_pos ) {
  
  // 调用规划器生成急停轨迹
  planner_manager_->EmergencyStop( stop_pos );

  // 发布急停轨迹(同上格式)
  auto info = &planner_manager_->local_data_;

  rebound_planner::Bspline bspline;
  bspline.order      = 3;
  bspline.start_time = info->start_time_;
  bspline.traj_id    = info->traj_id_;

  Eigen::MatrixXd pos_pts = info->position_traj_.getControlPoint();
  bspline.pos_pts.reserve( pos_pts.cols() );
  for (int i = 0; i < pos_pts.cols(); ++i) {
    geometry_msgs::Point pt;
    pt.x = pos_pts(0, i);
    pt.y = pos_pts(1, i);
    pt.z = pos_pts(2, i);
    bspline.pos_pts.push_back(pt);
  }

  Eigen::VectorXd knots = info->position_traj_.getKnot();
  bspline.knots.reserve(knots.rows()); 
  for (int i = 0; i < knots.rows(); ++i) {
    bspline.knots.push_back(knots(i));
  }

  bspline_pub_.publish(bspline);

  return true;
}

/**
 * @brief 计算局部规划的目标点和目标速度
 * 
 * 【核心功能】
 * 在全局轨迹上截取一段作为局部规划的目标。
 * 因为全局轨迹可能很长(从起点到终点),但局部规划只关心
 * 前方一段距离(planning_horizen),所以需要计算截取点。
 * 
 * 【计算方法: 球面-直线相交】
 * 
 * 全局路径用直线参数化: X(t) = M + N*t
 * 其中:
 * - M = init_pt_(初始全局规划起点)
 * - N = end_pt_ - init_pt_(全局方向向量)
 * 
 * 局部规划目标是在全局路径上距离无人机planning_horizen
 * 范围内的最远点。这等价于以无人机当前位置为球心,
 * planning_horizen为半径的球面与全局路径直线的交点。
 * 
 * 球面方程: (X - X₀)² = h²   (球心X₀=start_pt_, 半径h=planning_horizen)
 * 代入直线方程 X = M + N*t:
 *   (M + N*t - X₀)² = h²
 *   → (N·N)*t² + 2*(M-X₀)·N*t + (M-X₀)² - h² = 0
 *   → a*t² + b*t + c = 0 (二次方程)
 * 
 * 解出t后代入直线方程即得目标点。
 * 
 * 特殊情况:
 * - 如果无人机到终点的距离<planning_horizen: 目标直接设为终点
 * - 如果球面与直线不相交(无人机偏离全局路径太远): 
 *   取直线上的最近点(二次方程的极值点 t = -b/(2a))
 * 
 * 【目标速度计算】
 * 目标速度方向沿全局路径方向,大小根据到终点的距离衰减:
 * - 距离终点较远: 使用最大速度
 * - 接近终点: 按比例减速(保证能在终点前停下来)
 *   减速距离: d_brake = v_max² / (2*a_max) (匀减速公式)
 * 
 */
void ReboReplanFSM::getLocalTarget()
{
  // ---- 计算局部目标点 ----
  if ( ( end_pt_ - start_pt_ ).norm() < planning_horizen_ )
  {
    // 到终点的距离小于规划水平距离: 直接以终点为目标
    local_target_pt_ = end_pt_;
  }
  else
  {
    // 求解球面-直线相交问题
    Eigen::Vector3d M = init_pt_;           // 全局路径起点
    Eigen::Vector3d N = end_pt_ - init_pt_; // 全局路径方向: X = M + N*t
    Eigen::Vector3d X0 = start_pt_;         // 球心(无人机当前位置)
    double h = planning_horizen_;           // 球半径(规划水平距离)
    
    // 二次方程系数: a*t² + b*t + c = 0
    double a = N.squaredNorm();
    double b = 2*(M-X0).dot(N);
    double c = (M-X0).squaredNorm() - h*h;
    
    if ( b*b-4*a*c > 0 )  // 有实数解(球面与直线相交)
    {
      // 取正向交点(沿路径方向)
      double t = (-b + sqrt( b*b-4*a*c )) / (2*a);
      local_target_pt_ = M + N*t;
    }
    else  // 无交点(无人机偏离路径太远)
    {
      ROS_WARN("the drone goes too far to the reference line.");
      // 取路径上距无人机最近的点
      double t = -b/(2*a);
      local_target_pt_ = M + N*t;
    }
  }

  // ---- 计算局部目标速度 ----
  // 刹车距离: d_brake = v_max² / (2*a_max)
  double d_brake = (planner_manager_->pp_.max_vel_*planner_manager_->pp_.max_vel_)/(2*planner_manager_->pp_.max_acc_);
  
  if ( ( end_pt_ - local_target_pt_ ).norm() < d_brake )
  {
    // 在刹车距离内: 按比例减速(保证终点速度为零)
    local_target_vel_ = (end_pt_ - init_pt_).normalized() * 
        planner_manager_->pp_.max_vel_ * 
        (( end_pt_ - local_target_pt_ ).norm() / d_brake);
  }
  else
  {
    // 距离终点较远: 用最大速度巡航
    local_target_vel_ = (end_pt_ - init_pt_).normalized() * planner_manager_->pp_.max_vel_; 
  }
  
}

}  // namespace rebound_planner