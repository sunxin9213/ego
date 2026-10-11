#ifndef _DYN_A_STAR_H_
#define _DYN_A_STAR_H_

#include <iostream>
#include <ros/ros.h>
#include <ros/console.h>
#include <Eigen/Eigen>
#include <plan_env/grid_map.h>
#include <queue>

/**
 * @brief 无穷大常量
 *
 * 用于A*算法中表示不可达或未初始化的距离值。
 * 注意: 1>>20 是整数运算,值为 1/2^20 ≈ 0，此处应为设计缺陷，
 * 实际意图可能是 1e20 或 std::numeric_limits<double>::max()
 */
constexpr double inf = 1>>20;

/**
 * @brief 前向声明GridNode结构体
 */
struct GridNode;

/**
 * @brief GridNode指针类型别名
 */
typedef GridNode* GridNodePtr;

/**
 * @brief 网格节点结构体
 *
 * 表示A*搜索算法中的一个网格节点,包含位置索引、代价值估计和节点状态信息。
 * A*算法的核心数据结构,每个节点对应3D栅格地图中的一个体素。
 * 
 * 【算法原理】
 * A*算法通过维护每个节点的两个代价值来引导搜索:
 * - gScore: 从起点到当前节点的实际移动代价(已确定)
 * - fScore: gScore + 启发式估计代价,表示从起点经当前节点到终点的预估总代价
 * 每次从开放集中选择fScore最小的节点进行扩展,最终找到最优路径。
 */
struct GridNode
{
	/**
	 * @brief 节点状态枚举
	 */
	enum enum_state {
		OPENSET=1,    // 开放集:已发现但尚未扩展的节点
		CLOSEDSET=2,  // 关闭集:已被扩展过的节点(不再处理)
		UNDEFINED=3   // 未定义:尚未被探索的节点
	};

	/**
	 * @brief 搜索轮次(用于区分多次搜索)
	 * 
	 * 由于A*搜索可能被多次调用,每次调用时轮次递增。
	 * 通过比较rounds可以判断节点是否在当前搜索轮次中被处理过,
	 * 从而避免每次搜索前都需要重置所有节点状态。
	 */
	int rounds{0};
	
	/**
	 * @brief 节点当前状态(OPENSET/CLOSEDSET/UNDEFINED)
	 */
	enum enum_state state{UNDEFINED};
	
	/**
	 * @brief 节点在网格中的三维索引坐标 (x, y, z)
	 */
	Eigen::Vector3i index;
	
	/**
	 * @brief gScore: 从起点到当前节点的实际最小代价
	 * @brief fScore: fScore = gScore + 启发式代价, 用于优先队列排序
	 */
	double gScore{inf}, fScore{inf};
	
	/**
	 * @brief 前驱节点指针(用于路径回溯)
	 * 
	 * 当搜索到达终点后,通过cameFrom指针链反向回溯即可得到
	 * 从起点到终点的完整路径。
	 */
	GridNodePtr cameFrom{NULL};
};

/**
 * @brief 节点比较器
 *
 * 用于优先队列(最小堆)中按fScore从小到大排序。
 * 这是A*算法的关键:每次取出fScore最小的节点进行扩展。
 * fScore = gScore(已走距离) + 启发式估计(到终点的估计距离)
 */
class NodeComparator {
public:
  bool operator()(GridNodePtr node1, GridNodePtr node2) {
    return node1->fScore > node2->fScore;
  }
};

/**
 * @brief 动态A*路径搜索类
 *
 * 【算法原理】
 * 实现基于3D栅格地图的标准A*路径搜索算法,用于在已知环境中
 * 寻找从起点到终点的最短无碰撞路径。
 *
 * A*算法结合了Dijkstra算法(保证最优性)和贪心最佳优先搜索(提高效率)的优点:
 * - Dijkstra算法: f(n) = g(n), 只考虑已走代价,保证找到最短路径但效率低
 * - 贪心最佳优先: f(n) = h(n), 只考虑预估代价,效率高但不保证最优
 * - A*算法: f(n) = g(n) + h(n), 兼顾最优性和效率
 *
 * 当启发式函数h(n)满足可采纳性(即h(n) ≤ 真实代价)时,A*保证找到最优路径。
 * 本实现使用对角线距离作为启发式函数,在3D网格中具有可采纳性。
 */
class AStar
{
	private:
		// 栅格地图指针,用于查询障碍物占用状态
   		GridMap::Ptr grid_map_;

		/**
		 * @brief 快速坐标转网格索引(内联函数)
		 * @param x X坐标
		 * @param y Y坐标
		 * @param z Z坐标
		 * @param id_x 输出的X索引
		 * @param id_y 输出的Y索引
		 * @param id_z 输出的Z索引
		 */
		inline void coord2gridIndexFast(const double x, const double y, const double z, int & id_x, int & id_y, int & id_z);
		
		/**
		 * @brief 对角线启发式代价
		 * 
		 * 计算3D网格中的对角线距离(也称Chebyshev距离的变体)。
		 * 这是3D空间中允许对角线移动时的最短路径距离,
		 * 比欧几里得距离更适合作为网格搜索的启发式函数。
		 * 
		 * 计算公式: 先走尽可能多的三维对角线(diag步),
		 * 然后走二维对角线(min(dy,dz)或min(dx,dz)或min(dx,dy)步),
		 * 最后走直线段(abs(dy-dz)或abs(dx-dz)或abs(dx-dy)步)
		 * 
		 * @param node1 节点1
		 * @param node2 节点2
		 * @return 对角线距离
		 */
		double getDiagHeu(GridNodePtr node1, GridNodePtr node2);

		/**
		 * @brief 曼哈顿启发式代价
		 * 
		 * 曼哈顿距离 = |dx| + |dy| + |dz|
		 * 适用于只允许上下左右前后移动(不允许对角线移动)的场景。
		 * 作为启发式函数时,它总是≤真实代价,满足可采纳性。
		 * 
		 * @param node1 节点1
		 * @param node2 节点2
		 * @return 曼哈顿距离
		 */
		double getManhHeu(GridNodePtr node1, GridNodePtr node2);

		/**
		 * @brief 欧几里得启发式代价
		 * 
		 * 欧几里得距离 = sqrt(dx² + dy² + dz²)
		 * 在网格环境中,欧几里得距离通常小于实际移动距离,
		 * 因此满足可采纳性但不一致(可能高估),会导致A*扩展更多节点。
		 * 
		 * @param node1 节点1
		 * @param node2 节点2
		 * @return 欧几里得距离
		 */
		double getEuclHeu(GridNodePtr node1, GridNodePtr node2);

		/**
		 * @brief 获取启发式代价(带tie breaker)
		 * 
		 * 【Tie-Breaker机制】
		 * 当多个节点的fScore相同时,搜索会在这些等价节点之间徘徊,
		 * 降低搜索效率。通过引入一个略大于1的tie_breaker系数
		 * (1 + 1/10000),在保持可采纳性的前提下轻微扰动启发式值,
		 * 打破平局,使搜索更倾向于向目标方向前进。
		 * 
		 * @param node1 节点1
		 * @param node2 节点2
		 * @return 估计代价(带微小扰动)
		 */
		inline double getHeu(GridNodePtr node1, GridNodePtr node2);

		/**
		 * @brief 转换并调整起点和终点坐标
		 * 
		 * 将起点和终点的世界坐标转换为网格索引坐标。
		 * 如果起点或终点恰好位于障碍物内部,则沿着起点到终点的方向
		 * 逐步移动,直到找到一个空闲位置为止。
		 * 
		 * @param start_pt 起点世界坐标
		 * @param end_pt 终点世界坐标
		 * @param start_idx 输出的起点网格索引
		 * @param end_idx 输出的终点网格索引
		 * @return 是否成功转换
		 */
		bool ConvertToIndexAndAdjustStartEndPoints(const Eigen::Vector3d start_pt, const Eigen::Vector3d end_pt, Eigen::Vector3i &start_idx, Eigen::Vector3i &end_idx);

		/**
		 * @brief 网格索引转世界坐标
		 * 
		 * 将网格索引映射回世界坐标(取体素的中心位置)
		 * 
		 * @param index 网格索引
		 * @return 世界坐标
		 */
		inline Eigen::Vector3d Index2Coord(const Eigen::Vector3i &index) const;

		/**
		 * @brief 世界坐标转网格索引
		 * 
		 * 将世界坐标映射到网格索引,同时检查是否超出地图边界。
		 * 
		 * @param pt 世界坐标
		 * @param idx 输出的网格索引
		 * @return 是否转换成功(索引在地图范围内)
		 */
		inline bool Coord2Index(const Eigen::Vector3d &pt, Eigen::Vector3i &idx) const;

		/**
		 * @brief 检查位置是否被占用(使用膨胀后的地图)
		 * 
		 * 通过查询膨胀占据地图来判断一个位置是否不可通行。
		 * 使用膨胀地图而不是原始地图,是为了给障碍物施加安全缓冲区,
		 * 确保生成的路径与障碍物保持安全距离。
		 * 
		 * @param pos 位置坐标
		 * @return true 表示被占用(不可通行)
		 */
		inline bool checkOccupancy( const Eigen::Vector3d &pos ) { return (bool)grid_map_->getInflateOccupancy(pos); }

		/**
		 * @brief 通过前驱指针链回溯路径
		 * 
		 * 从目标节点开始,沿着cameFrom指针链一直回溯到起点,
		 * 得到完整的路径节点序列(从起点到终点)。
		 * 
		 * @param current 终点节点
		 * @return 路径节点指针列表(从起点到终点)
		 */
		std::vector<GridNodePtr> retrievePath(GridNodePtr current);

		// 步长(网格分辨率)和步长倒数
		double step_size_, inv_step_size_;
		// 地图中心点(起点和终点的中点)
		Eigen::Vector3d center_;
		// 中心索引和网格池大小
		Eigen::Vector3i CENTER_IDX_, POOL_SIZE_;
		// tie breaker系数:1+1/10000,防止fScore相同时的搜索震荡
		const double tie_breaker_ = 1.0 + 1.0 / 10000;

		// A*搜索找到的网格路径(存储节点指针)
		std::vector<GridNodePtr> gridPath_;

		// 3D网格节点地图: 三维指针数组,存储所有网格节点的指针
		GridNodePtr *** GridNodeMap_;
		// 开放集(优先队列/最小堆),按fScore从小到大排列
		std::priority_queue< GridNodePtr, std::vector<GridNodePtr>, NodeComparator > openSet_;

		// 搜索轮次计数器,每次调用AstarSearch递增
		int rounds_{0};

	public:
		typedef std::shared_ptr<AStar> Ptr;

		AStar(){};
		~AStar();

		/**
		 * @brief 初始化网格地图
		 * 
		 * 分配网格节点池的内存空间,并关联占据地图。
		 * 网格池是一个以CENTER_IDX_为中心的固定大小的三维数组。
		 * 
		 * @param occ_map 占据地图指针
		 * @param pool_size 网格池在三个维度上的大小
		 */
		void initGridMap( GridMap::Ptr occ_map, const Eigen::Vector3i pool_size );

		/**
		 * @brief A*路径搜索主函数
		 * 
		 * 【算法执行流程】
		 * 1. 将起点和终点转换为网格索引,如果被占用则调整位置
		 * 2. 将起点加入开放集(优先级队列)
		 * 3. 循环执行：
		 *    a. 从开放集中取出fScore最小的节点作为当前节点
		 *    b. 如果当前节点是终点,则路径找到,返回成功
		 *    c. 将当前节点移入关闭集
		 *    d. 遍历26邻域(3x3x3邻域,排除中心)
		 *       - 跳过边界外的节点
		 *       - 跳过已在关闭集中的节点
		 *       - 跳过被障碍物占用的节点
		 *       - 计算移动代价,如果节点未探索则加入开放集
		 *       - 如果节点已在开放集中且新路径代价更小,则更新
		 * 4. 如果开放集为空仍未找到终点,则搜索失败
		 * 
		 * @param step_size 搜索步长(网格分辨率)
		 * @param start_pt 起点世界坐标
		 * @param end_pt 终点世界坐标
		 * @return true 找到路径, false 搜索失败
		 */
		bool AstarSearch( const double step_size, Eigen::Vector3d start_pt, Eigen::Vector3d end_pt);

		/**
		 * @brief 获取A*搜索得到的路径(世界坐标序列)
		 * 
		 * 将gridPath_中的节点指针转换为世界坐标序列。
		 * 
		 * @return 路径点列表(从起点到终点)
		 */
		std::vector<Eigen::Vector3d> getPath();
};

/**
 * @brief 获取带tie breaker的启发式代价
 * 
 * 使用对角线距离作为基础启发式函数,并乘以tie_breaker系数。
 * tie_breaker略大于1,在保持可采纳性的前提下使搜索更快收敛。
 */
inline double AStar::getHeu(GridNodePtr node1, GridNodePtr node2)
{
    return tie_breaker_ * getDiagHeu(node1, node2);
}

/**
 * @brief 网格索引转世界坐标
 * 
 * 公式: coord = (index - CENTER_IDX_) * step_size_ + center_
 * 即先计算相对偏移,乘以分辨率,再加中心点坐标
 */
inline Eigen::Vector3d AStar::Index2Coord(const Eigen::Vector3i &index) const
{
	return ( ( index - CENTER_IDX_ ).cast<double>() * step_size_ ) + center_;
};

/**
 * @brief 世界坐标转网格索引
 * 
 * 公式: idx = ((pt - center_) * inv_step_size_ + 0.5).cast<int>() + CENTER_IDX_
 * 加0.5是为了四舍五入取整
 * 
 * @param pt 世界坐标
 * @param idx 输出的网格索引
 * @return true 如果索引在地图范围内
 */
inline bool AStar::Coord2Index(const Eigen::Vector3d &pt, Eigen::Vector3i &idx) const
{
	idx = ((pt - center_) * inv_step_size_ + Eigen::Vector3d(0.5,0.5,0.5)).cast<int>() + CENTER_IDX_;

	if (   idx(0) < 0 || idx(0) >= POOL_SIZE_(0) 
		|| idx(1) < 0 || idx(1) >= POOL_SIZE_(1) 
		|| idx(2) < 0 || idx(2) >= POOL_SIZE_(2) ) 	
	{
		ROS_ERROR("Ran out of pool, index=%d %d %d",idx(0),idx(1),idx(2));
		return false;	
	}
	
	return true;
};

#endif