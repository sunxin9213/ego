# uav_simulator（无人机仿真环境）功能梳理

> 本文件用于梳理 `uav_simulator` 文件夹的整体功能与作用流程，**不含算法核心**（规划/避障算法在 `planner`、`ego_planner` 等目录）。

## 整体定位

`uav_simulator` 是 EGO 算法的**仿真测试环境**。它提供一套完整的"虚拟无人机世界"：**地图、感知、飞行控制器、动力学仿真、消息与可视化工具**——让算法跑在仿真里验证，而不需要真实硬件。本身不包含规划逻辑。

## 分层结构

它分为四层，从"物理世界"流向"算法"再流回"物理世界"：

```
虚拟地图(map_generator / local_sensing)
        →  传感器/感知(local_sensing)
        →  算法接收到里程计/地图/点云
        →  so3_control 计算控制指令
        →  so3_quadrotor_simulator 物理动力学仿真
        →  发布新的 odom，循环
```

## 各模块功能

### 1. 物理仿真层

**[so3_quadrotor_simulator]** — 四旋翼动力学仿真器（模拟器的心脏）
- `src/dynamics/Quadrotor.cpp`：建立四旋翼物理模型（质量、惯性、电机时间常数、螺旋桨推力系数等），用 **boost odeint 数值积分**（RK4）求解运动微分方程，得到位置、速度、姿态 R、角速度、电机转速等 22 维内部状态。
- `src/quadrotor_simulator_so3.cpp`：ROS 节点。订阅 `cmd`（SO3 控制指令）和扰动力/力矩，按仿真频率（默认 1000Hz）接收控制量 → 反解出四个螺旋桨转速 `rpm` → 调用 `quad.step(dt)` 推进物理仿真 → 以 100Hz 发布 `odom`（里程计）和 `imu` 数据给下游算法。
- 注意：内部 `getControl` 是"逆动力学解算"（把期望力/力矩换算成电机转速），**不是控制器本身**。

### 2. 场景/地图层

**[map_generator]** — 随机障碍地图生成器
- `src/random_forest_sensing.cpp`：用随机数生成一片"森林"（随机的柱状/方块障碍物，带圆柱碰撞半径），并发布**全局静态地图点云** 和**局部地图点云**，作为感知和避障的地图真值来源。同时监听里程计，生成以无人机为中心的局部图。

### 3. 感知层（模拟传感器）

**[local_sensing]** — 模拟无人机"眼睛"
- `src/depth_render.cu`：基于 CUDA，用 GPU 从全局地图实时渲染深度图（模拟深度相机）。
- `src/pcl_render_node.cpp`、`src/pointcloud_render_node.cpp`：把渲染出的深度图转成 **PCL 点云**，并配准到世界坐标系，发布"机载传感器看到的局部地图点云"。EGO 算法的局部地图就是从这些话题拿数据。
- `src/euroc.cpp`：负责读入 EuRoC 真实数据集（用于更真实测试）。

### 4. 控制层

**[so3_control]** — 底层姿态控制器
- `src/SO3Control.cpp`：把期望位置/速度/加速度经 **PD + 前馈 + 重力补偿** 解算成推力 `force_` 和期望姿态四元数 `orientation_`（SO(3) 控制）。
- `src/so3_control_nodelet.cpp`：订阅上游算法发布的 `PositionCommand`（EGO 规划器的轨迹点），算出力/姿态后发布 `SO3Command` 给仿真器——串起了"算法 → 执行器"。

### 5. 工具/消息/可视化层（`Utils` 子包）

一组公共基础设施，为上面所有模块提供支撑：

- **[quadrotor_msgs]**：自定义 ROS 消息定义（`PositionCommand`、`SO3Command`、`Odometry`、`Gains` 等）和序列化/反序列化库，是整个系统共享的"通信协议"。
- **[pose_utils]**：姿态转换工具库。
- **[uav_utils]**：通用数学工具（坐标系转换、几何、旋转矩阵 ↔ 欧拉角）；`scripts` 下有 odom 转欧拉角、发 odom、TF 辅助等调试脚本。
- **[waypoint_generator]**：航点生成器。提供 `manual`（RViz 手动点击设置目标）和 `series`（按 yaml 参数分段生成航点，可设偏航和时间）两种模式，以 `goal` / `waypoints` 消息发布给 EGO 规划器作为任务目标。
- **[rviz_plugins]**：RViz 可视化插件（地图显示、`goal_tool` 目标点选取工具、`pose_tool` 位姿工具等），方便手动给算法下发目标、观察地图。
- **[odom_visualization]**：用无人机 3D 网格模型（hummingbird.mesh）在 RViz 里显示飞行姿态。
- **[multi_map_server]**：多分辨率/多张地图管理服务，把多张地图合成/发布成组合地图消息。

> 说明：`odom_visualization`、`pose_utils`、`multi_map_server` 等目录里带 `.so`、`.pyc`、`*.cpp~`、`ode/`（完整 boost odeint 第三方库）等，属于编译产物或第三方依赖，不是核心逻辑。

## 一句话流程总结

```
waypoint_generator 或 RViz 设定目标
        ↓ goal/waypoints
EGO 规划算法（不在本目录）
        ↓ PositionCommand（期望轨迹）
so3_control：计算期望推力与姿态
        ↓ SO3Command
so3_quadrotor_simulator：物理动力学仿真 → odom/imu
        ↓
map_generator 全局地图 + local_sensing 机载传感器点云
        → 反馈给规划算法闭环
        + RViz 插件/odom_visualization 可视化观察
```

简单说，`uav_simulator` 就是一套 **"虚拟无人机 + 虚拟世界 + 传感器 + 底层控制器" 的闭环仿真平台**，供 EGO 规划算法在上面跑通和演示，本身不含规划逻辑。