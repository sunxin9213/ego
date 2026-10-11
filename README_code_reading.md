# EGO Planner 代码阅读指南

> 本文档梳理 EGO（反弹重规划 / Rebound Re-planning）无人机局部运动规划器的代码阅读逻辑，帮助按依赖关系和运行流程快速定位核心模块。

## 一、项目定位

这是一个基于 ROS + Eigen 的**无人机局部运动规划器**（本地代码中的 `rebound_planner` 包），属于 EGO（Elastic Goal-based Online planner）一类算法。核心思路：

- 用**全局多项式轨迹**做参考基准（粗引导）——见 `poly_traj`
- 用**局部三次均匀 B 样条**表示实际执行轨迹——见 `bspline_opt`
- 用**弹性带 + 梯度下降**优化控制点，兼顾平滑、避障、动力学可行——见 `bspline_optimizer`
- 用**有限状态机（FSM）**驱动"生成/执行/重规划/急停"的完整生命周期——见 `plan_manage`

> 注：本代码为 EGO 的独立重实现（类似 `Fast-Planner` 的 EGO/Re-Opt 思路），只保留了运动学 B 样条优化路径；`dyn_a_star` 仅用于为碰撞段生成引导路径，非动力学搜索。

---

## 二、整体架构（模块与依赖关系）

代码主体在 `src/planner/` 下，建议按"底层→顶层"方向阅读——先理解数据结构，再理解算法，最后理解调度。

```
planner/
├── plan_env/      # 环境感知层：栅格地图 GridMap（体素、膨胀、射线建图、SDF）
├── poly_traj/     # 轨迹表示层：多项式轨迹（最小Snap闭式解，生成全局参考）
├── path_searching/# 搜索层：dyn_a_star（A*，为碰撞段生成无碰撞引导路径）
├── bspline_opt/   # 轨迹表示+优化层：UniformBspline（B样条）+ BsplineOptimizer（弹性带优化）
├── plan_manage/   # 调度与控制层：FSM + PlannerManager + 主节点 + 轨迹服务器
└── traj_utils/    # 可视化工具（PlanningVisualization）
```

依赖方向（建议阅读顺序）：

```
plan_env (GridMap)
   └──> poly_traj (minSnapTraj)
   └──> path_searching (AStar)
             └──> bspline_opt (UniformBspline, BsplineOptimizer)
                        └──> plan_manage (PlannerManager → ReboReplanFSM → 节点)
```

---

## 三、推荐阅读顺序

### Stage 0：入口与启动 — [rebound_planner_node.cpp](file:///f:/ego/src/planner/plan_manage/src/rebound_planner_node.cpp)

```cpp
int main() {
  ReboReplanFSM rebo_replan;
  rebo_replan.init(nh);   // 一切从这里开始
  ros::spin();
}
```
只做一件事：实例化并初始化 `ReboReplanFSM`，然后进入 ROS 事件循环。整个系统全部由 `FSM` 内部的定时器与回调驱动。

---

### Stage 1：环境感知层（先理解"世界长什么样"）

#### [grid_map.h](file:///f:/ego/src/planner/plan_env/include/plan_env/grid_map.h)
- `MappingParameters`：地图原点/尺寸/分辨率/膨胀距离等参数结构。
- `GridMap` 类：体素地图、`getInflateOccupancy()`（膨胀后占用查询，避障核心）、射线建图、局部地图更新。

#### [raycast.h](file:///f:/ego/src/planner/plan_env/include/plan_env/raycast.h)
- 射线光栅化，用于将深度相机点云投影到占据栅格。

> 阅读重点：理解"占用查询"接口 `getInflateOccupancy(pts)`，它是后面 A*、碰撞检测、优化距离代价都要调用的底层函数。

---

### Stage 2：全局轨迹生成（"往哪走"的参考）

#### [polynomial_traj.cpp](file:///f:/ego/src/planner/poly_traj/src/polynomial_traj.cpp)
- `minSnapTraj()`：**最小加加速度（Snap）轨迹的闭式解法**，生成全局参考多项式轨迹。
- 8 步流程：构造导数向量 D → 映射矩阵 A → 选择矩阵 C → 最小 Snap 矩阵 Q → 求 R → 分解 → 解自由变量 → 还原系数。
- 每段用 5 次多项式（6 个系数恰好满足位置/速度/加速度边界）。代码中已对每一步做了详细数学注释。

> 阅读重点：这段是"全局线路"的来源，局部规划会沿着它截取目标点。它只被 `PlanGlobalTraj` 调用。

---

### Stage 3：搜索层（避障引导线）

#### [dyn_a_star.cpp](file:///f:/ego/src/planner/path_searching/src/dyn_a_star.cpp)
- `AStar` 类：三维栅格 A* 搜索。
- 在 B 样条优化中充当"弹性带锚点"——为每个碰撞段搜索一条绕过障碍物的无碰撞路径，作为控制点被"弹"过去的方向（见 Stage 4 步骤 4）。

> 阅读重点：理解 `AstarSearch()/getPath()` 的输出，它供 `BsplineOptimizer::initControlPoints()` 使用。

---

### Stage 4：B 样条与优化层（核心算法）

#### [uniform_bspline.h](file:///f:/ego/src/planner/bspline_opt/include/bspline_opt/uniform_bspline.h) / [.cpp](file:///f:/ego/src/planner/bspline_opt/src/uniform_bspline.cpp)
- 均匀 B 样条曲线类：控制点矩阵、节点向量、基函数求值 `evaluateDeBoorT()/evaluate()`。
- 关键性质：局部支撑、`getDerivativeControlPoints()` 求导（速度/加速度轨迹）、`setPhysicalLimits()` 与 `checkFeasibility()`（速度/加速度约束检查）、`reparamBspline` 时间重分配。

#### [bspline_optimizer.h](file:///f:/ego/src/planner/bspline_opt/include/bspline_opt/bspline_optimizer.h) / [.cpp](file:///f:/ego/src/planner/bspline_opt/src/bspline_optimizer.cpp)
这是**核心优化器**。核心接口：
- `initControlPoints()`：**弹性带/反弹预处理**，分 4 步：
  1. 沿初始控制点逐点检测、切分碰撞段（in_id/out_id）；
  2. 对每段跑 A* 得到无碰撞引导路径；
  3. 计算各段可优化控制点的索引边界（bounds，保护首尾 order 个控制点不变）；
  4. 计算"基点 + 方向向量"（控制点指向 A* 路径的交点方向），形成弹性带。
- `BsplineOptimizeTrajRebound()` / `BsplineOptimizeTrajRefine()`：优化入口。

**代价函数**（文件头部有详细注释）：
```
J = λ1·J_smooth(平滑/jerk) + λ2·J_distance(距离/避障)
    + λ3·J_feasibility(可行性/速度加速度) + λ4·J_fitness(拟合/追踪参考)
```

#### [gradient_descent_optimizer.cpp](file:///f:/ego/src/planner/bspline_opt/src/gradient_descent_optimizer.cpp) / [lbfgs.hpp](file:///f:/ego/src/planner/bspline_opt/include/bspline_opt/lbfgs.hpp)
- 梯度下降求解器（含 L-BFGS 实现），实际执行控制点梯度下降优化。

> 阅读重点：把"弹性带"的四个步骤和四类代价项看懂，就理解了 EGO 如何同时在平滑、避障、可行之间折中。

---

### Stage 5：调度与控制层（系统"大脑"）

#### [planner_manager.h](file:///f:/ego/src/planner/plan_manage/include/plan_manage/planner_manager.h) / [.cpp](file:///f:/ego/src/planner/plan_manage/src/planner_manager.cpp)
- `ReboundPlannerManager`：协调地图、全局轨迹、局部轨迹、优化器。
- 关键方法：
  - `planGlobalTraj()`：用 `minSnapTraj` 生成全局参考轨迹。
  - `initPlanModules()`：读参数、建 `GridMap`、`BsplineOptimizer`、`AStar`。
  - **`reboundReplan()`**：局部重规划主流程（三段式），是核心：
    1. **INIT**：构造初始控制点。根据 `flag_polyInit` 选择从多项式轨迹或上一段轨迹采样（`parameterizeToBspline`）；支持随机扰动跳局部最优。
    2. **OPTIMIZE**：`BsplineOptimizeTrajRebound`（弹性带优化）。
    3. **REFINE**：`checkFeasibility` 失败则 `reparamBspline` 重分配时间后再 `BsplineOptimizeTrajRefine`。
  - `updateTrajInfo()`：保存优化结果到 `LocalTrajData`（位置/速度/加速度轨迹）。
- 数据结构：[plan_container.hpp](file:///f:/ego/src/planner/plan_manage/include/plan_manage/plan_container.hpp) 中的 `GlobalTrajData / LocalTrajData`，负责全局轨迹与局部 B 样条轨迹的统一时间管理与状态查询。

#### [rebo_replan_fsm.h](file:///f:/ego/src/planner/plan_manage/include/plan_manage/rebo_replan_fsm.h) / [.cpp](file:///f:/ego/src/planner/plan_manage/src/rebo_replan_fsm.cpp)
- **有限状态机**，整个系统的调度核心。状态：
  ```
  INIT → WAIT_TARGET → GEN_NEW_TRAJ → EXEC_TRAJ → (循环)
                        ↑                ↓
                   REPLAN_TRAJ  ←  EMERGENCY_STOP
  ```
- 定时器驱动：
  - `exec_timer_`（10ms）→ `execFSMCallback()`：推进状态机。
  - `safety_timer_`（50ms）→ `checkCollisionCallback()`：轨迹碰撞检测，按紧迫度分级（>0.8s 重规划，<0.8s 急停）。
- 回调：`waypointCallback()`（接收目标）、`odometryCallback()`（里程计）。
- `callReboundReplan()`：串联 `getLocalTarget()` → `reboundReplan()` → 发布 B 样条。
- `getLocalTarget()`：在全局直线上按 `planning_horizen` 截取局部目标点（球面-直线求交），并计算目标速度。

#### [traj_server.cpp](file:///f:/ego/src/planner/plan_manage/src/traj_server.cpp)
- 接收 FSM 发布的 `/planning/bspline`，转换成 `PositionCommand` 发给底层控制器（P 控制），实现轨迹跟踪。

---

### Stage 6：可视化

- [PlanningVisualization](file:///f:/ego/src/planner/plan_manage/src/planning_visualization.cpp)（以及 `traj_utils`），用于在 Rviz 显示全局路径、初始路径、A* 路径、最优轨迹、目标点等，辅助调试。

---

## 四、一图读懂数据流

```
[目标航点/用户点击]  waypointCallback
        ↓ planGlobalTraj (minSnapTraj → 全局多项式轨迹)
[Rviz显示全局路径]  ←──────────────────────────────┘
        ↓ execFSM (GEN_NEW_TRAJ / REPLAN_TRAJ)
        ↓ getLocalTarget  在全局线上截取局部目标点
        ↓ reboundReplan [PlannerManager]
        ↓      ├─ INIT   : 采样构造B样条初始控制点 (多项式或上段轨迹)
        ↓      ├─ OPTIMIZE: initControlPoints(切段+A*) → BsplineOptimizeTrajRebound
        ↓      └─ REFINE : Feasibility失败→重分配时间→BsplineOptimizeTrajRefine
        ↓ updateTrajInfo  → LocalTrajData
        ↓ 发布 /planning/bspline  (B样条消息)
[traj_server] 接收 → 计算位置/速度 → 发布 PositionCommand → 底层控制器
        ↑
[odometryCallback / 碰撞检测 checkCollisionCallback] 实时反馈,决定是否重规划/急停
```

---

## 五、关键概念速查

| 概念 | 含义 | 位置 |
|------|------|------|
| 全局轨迹 | 最小 Snap 多项式，粗参考 | `poly_traj` |
| 局部轨迹 | 三次均匀 B 样条，实际执行 | `bspline_opt` |
| 弹性带/Rebound | 碰障时用 A* 引导线把控制点"弹"回安全区 | `bspline_optimizer.initControlPoints` |
| 代价加权 | 平滑/距离/可行/拟合 四类项加权 | `bspline_optimizer.cpp` 注释 |
| FSM | 生成/执行/重规划/急停 状态机 | `rebo_replan_fsm.cpp` |
| 局部目标 | 沿全局线按 `planning_horizen` 截取 | `getLocalTarget()` |
| 碰撞分级响应 | 时间充裕→重规划；紧迫→急停 | `checkCollisionCallback()` |

---

## 六、主要参数（改行为去哪找）

- 地图：`grid_map/resolution`、`grid_map/obstacles_inflation`
- 运动约束：`manager/max_vel`、`manager/max_acc`、`manager/max_jerk`、`manager/control_points_distance`
- 优化权重：`optimization/lambda_smooth`、`lambda_collision`、`lambda_feasibility`、`lambda_fitness`、`dist0`
- FSM：`fsm/thresh_replan`、`fsm/thresh_no_replan`、`fsm/planning_horizen`、`fsm/flight_type`
- 全部在 [rebo_algorithm.xml](file:///f:/ego/src/planner/plan_manage/launch/rebo_algorithm.xml) / [rebo_replan.launch](file:///f:/ego/src/planner/plan_manage/launch/rebo_replan.launch) 中集中配置。

---

## 七、模拟环境（可选辅助）

真正的仿真闭环还包括 `uav_simulator`（`so3_control` 控制器、`so3_quadrotor_simulator` 无人机模型、`pointcloud_render_node` 深度渲染、`multi_map_server`、`waypoint_generator` 等）。但**理解规划算法逻辑只需关注 `src/planner/`**，模拟器只是数据来源与执行对象。