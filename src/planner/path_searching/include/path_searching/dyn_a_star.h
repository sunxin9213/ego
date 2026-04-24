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
 * 用于A*算法中表示不可达或未初始化的距离值
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
 * 表示A*搜索算法中的一个网格节点,包含位置索引、代价估计和状态信息。
 */
struct GridNode
{
	/**
	 * @brief 节点状态枚举
	 */
	enum enum_state {
		OPENSET=1,    // 在开放集中
		CLOSEDSET=2,  // 在关闭集中
		UNDEFINED=3   // 未定义
	};

	/**
	 * @brief 搜索轮次(用于区分多次搜索)
	 */
	int rounds{0};
	
	/**
	 * @brief 节点当前状态
	 */
	enum enum_state state{UNDEFINED};
	
	/**
	 * @brief 节点在网格中的索引坐标
	 */
	Eigen::Vector3i index;
	
	/**
	 * @brief gScore:从起点到当前节点的实际代价
	 * @brief fScore:从起点经当前节点到终点的估计总代价
	 */
	double gScore{inf}, fScore{inf};
	
	/**
	 * @brief 指向当前节点的前驱节点(用于路径回溯)
	 */
	GridNodePtr cameFrom{NULL};
};

/**
 * @brief 节点比较器
 *
 * 用于优先队列中按fScore从小到大排序
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
 * 实现基于3D栅格地图的A*路径搜索算法,支持动态环境下的快速重规划。
 * 该算法用于在已知环境中寻找从起点到终点的最短路径。
 */
class AStar
{
	private:
		// 栅格地图指针
   		GridMap::Ptr grid_map_;

		/**
		 * @brief 快速坐标转网格索引
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
		 * @param node1 节点1
		 * @param node2 节点2
		 * @return 对角线距离
		 */
		double getDiagHeu(GridNodePtr node1, GridNodePtr node2);

		/**
		 * @brief 曼哈顿启发式代价
		 * @param node1 节点1
		 * @param node2 节点2
		 * @return 曼哈顿距离
		 */
		double getManhHeu(GridNodePtr node1, GridNodePtr node2);

		/**
		 * @brief 欧几里得启发式代价
		 * @param node1 节点1
		 * @param node2 节点2
		 * @return 欧几里得距离
		 */
		double getEuclHeu(GridNodePtr node1, GridNodePtr node2);

		/**
		 * @brief 获取启发式代价(带tie breaker)
		 * @param node1 节点1
		 * @param node2 节点2
		 * @return 估计代价
		 */
		inline double getHeu(GridNodePtr node1, GridNodePtr node2);

		/**
		 * @brief 转换并调整起点和终点坐标
		 * @param start_pt 起点坐标
		 * @param end_pt 终点坐标
		 * @param start_idx 输出的起点索引
		 * @param end_idx 输出的终点索引
		 * @return 是否成功
		 */
		bool ConvertToIndexAndAdjustStartEndPoints(const Eigen::Vector3d start_pt, const Eigen::Vector3d end_pt, Eigen::Vector3i &start_idx, Eigen::Vector3i &end_idx);

		/**
		 * @brief 索引转坐标
		 * @param index 索引
		 * @return 坐标
		 */
		inline Eigen::Vector3d Index2Coord(const Eigen::Vector3i &index) const;

		/**
		 * @brief 坐标转索引
		 * @param pt 坐标
		 * @param idx 输出的索引
		 * @return 是否成功
		 */
		inline bool Coord2Index(const Eigen::Vector3d &pt, Eigen::Vector3i &idx) const;

		/**
		 * @brief 检查位置占用状态(使用膨胀后的地图)
		 * @param pos 位置坐标
		 * @return 是否被占用
		 */
		inline bool checkOccupancy( const Eigen::Vector3d &pos ) { return (bool)grid_map_->getInflateOccupancy(pos); }

		/**
		 * @brief 检索路径
		 * @param current 当前节点
		 * @return 路径节点列表
		 */
		std::vector<GridNodePtr> retrievePath(GridNodePtr current);

		// 步长和步长倒数
		double step_size_, inv_step_size_;
		// 地图中心点
		Eigen::Vector3d center_;
		// 中心索引和池大小
		Eigen::Vector3i CENTER_IDX_, POOL_SIZE_;
		// tie breaker系数(避免相同fScore的节点排序问题)
		const double tie_breaker_ = 1.0 + 1.0 / 10000;

		// 网格路径
		std::vector<GridNodePtr> gridPath_;

		// 3D网格节点地图
		GridNodePtr *** GridNodeMap_;
		// 开放集(优先队列)
		std::priority_queue< GridNodePtr, std::vector<GridNodePtr>, NodeComparator > openSet_;

		// 搜索轮次
		int rounds_{0};

	public:
		typedef std::shared_ptr<AStar> Ptr;

		AStar(){};
		~AStar();

		/**
		 * @brief 初始化网格地图
		 * @param occ_map 占据地图指针
		 * @param pool_size 网格池大小
		 */
		void initGridMap( GridMap::Ptr occ_map, const Eigen::Vector3i pool_size );

		/**
		 * @brief A*搜索
		 * @param step_size 搜索步长
		 * @param start_pt 起点
		 * @param end_pt 终点
		 * @return 是否找到路径
		 */
		bool AstarSearch( const double step_size, Eigen::Vector3d start_pt, Eigen::Vector3d end_pt);

		/**
		 * @brief 获取搜索路径
		 * @return 路径点列表
		 */
		std::vector<Eigen::Vector3d> getPath();
};

inline double AStar::getHeu(GridNodePtr node1, GridNodePtr node2)
{
    return tie_breaker_ * getDiagHeu(node1, node2);
}

inline Eigen::Vector3d AStar::Index2Coord(const Eigen::Vector3i &index) const
{
	return ( ( index - CENTER_IDX_ ).cast<double>() * step_size_ ) + center_;
};

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
