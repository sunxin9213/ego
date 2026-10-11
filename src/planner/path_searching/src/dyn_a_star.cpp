#include "path_searching/dyn_a_star.h"

using namespace std;
using namespace Eigen;

AStar::~AStar()
{
    /* 释放三维网格节点池中所有节点的内存 */
    for(int i = 0; i < POOL_SIZE_(0); i++)
        for(int j = 0; j < POOL_SIZE_(1); j++)
            for( int k = 0; k < POOL_SIZE_(2); k++)
                delete GridNodeMap_[i][j][k];
}

void AStar::initGridMap( GridMap::Ptr occ_map, const Eigen::Vector3i pool_size )
{   
    /**
     * 初始化网格地图
     * 
     * 分配一个三维的GridNode指针数组(网格节点池),尺寸为pool_size。
     * 网格池以CENTER_IDX_为中心对称分布,使得A*搜索可以在地图中
     * 任意位置快速定位节点,避免动态分配的开销。
     */
    POOL_SIZE_ = pool_size;
    CENTER_IDX_ = pool_size / 2;
    
    /* 动态分配三维指针数组 */
    GridNodeMap_ = new GridNodePtr** [POOL_SIZE_(0)];
    for(int i = 0; i < POOL_SIZE_(0); i++)
    {
        GridNodeMap_[i] = new GridNodePtr* [POOL_SIZE_(1)];
        for(int j = 0; j < POOL_SIZE_(1); j++)
        {
            GridNodeMap_[i][j] = new GridNodePtr [POOL_SIZE_(2)];
            for( int k = 0; k < POOL_SIZE_(2); k++)
            {
                /* 预先创建所有节点,避免运行时频繁new/delete */
                GridNodeMap_[i][j][k] = new GridNode;
            }
        }
    }

    grid_map_ = occ_map;
}

double AStar::getDiagHeu(GridNodePtr node1, GridNodePtr node2)
{   
    /**
     * 对角线启发式距离(3D网格中的最短路径距离)
     * 
     * 【算法原理】
     * 在允许三维对角线移动的网格中,从一点到另一点的最短路径可以分解为：
     * 1. 三维对角线移动: 同时在x,y,z三个方向移动,步长√3
     * 2. 二维对角线移动: 在两个方向移动,步长√2
     * 3. 一维直线移动: 在一个方向移动,步长1
     * 
     * 先走尽可能多的三维对角线,然后走二维对角线,最后走直线段,
     * 这样得到的距离是对角线移动代价的精确下界。
     */
    double dx = abs(node1->index(0) - node2->index(0));
    double dy = abs(node1->index(1) - node2->index(1));
    double dz = abs(node1->index(2) - node2->index(2));

    double h = 0.0;
    int diag = min(min(dx, dy), dz);  // 可以走三维对角线的步数
    dx -= diag;
    dy -= diag;
    dz -= diag;

    if (dx == 0) {
        /* 在yz平面上走二维对角线 + 直线 */
        h = 1.0 * sqrt(3.0) * diag + sqrt(2.0) * min(dy, dz) + 1.0 * abs(dy - dz);
    }
    if (dy == 0) {
        h = 1.0 * sqrt(3.0) * diag + sqrt(2.0) * min(dx, dz) + 1.0 * abs(dx - dz);
    }
    if (dz == 0) {
        h = 1.0 * sqrt(3.0) * diag + sqrt(2.0) * min(dx, dy) + 1.0 * abs(dx - dy);
    }
    return h;
}

double AStar::getManhHeu(GridNodePtr node1, GridNodePtr node2)
{   
    /**
     * 曼哈顿启发式距离
     * h = |dx| + |dy| + |dz|
     * 适用于仅允许四方向/六方向移动的网格
     */
    double dx = abs(node1->index(0) - node2->index(0));
    double dy = abs(node1->index(1) - node2->index(1));
    double dz = abs(node1->index(2) - node2->index(2));

    return dx + dy + dz;
}

double AStar::getEuclHeu(GridNodePtr node1, GridNodePtr node2)
{   
    /**
     * 欧几里得启发式距离
     * h = sqrt(dx² + dy² + dz²)
     * 在网格中通常低估真实代价,可采纳但效率较低
     */
    return (node2->index - node1->index).norm();
}

vector<GridNodePtr> AStar::retrievePath(GridNodePtr current)
{   
    /**
     * 回溯路径
     * 
     * 从终点节点开始,沿着cameFrom指针链一直回溯到起点。
     * 因为cameFrom是从起点指向终点的方向,所以回溯后需要反转。
     */
    vector<GridNodePtr> path;
    path.push_back(current);

    while(current->cameFrom != NULL)
    {
        current = current -> cameFrom;
        path.push_back( current );
    }

    return path;
}

bool AStar::ConvertToIndexAndAdjustStartEndPoints(Vector3d start_pt, Vector3d end_pt, Vector3i &start_idx, Vector3i &end_idx)
{
    /**
     * 坐标转换和位置调整
     * 
     * 如果起点或终点在障碍物内部,沿起点到终点的方向逐步移动,
     * 直到找到一个空闲位置。这是在实际应用中常见的处理方式:
     * 由于地图更新和定位误差,规划的起点可能恰好落在障碍物上。
     */
    if ( ! Coord2Index(start_pt, start_idx) || ! Coord2Index(end_pt, end_idx) )
        return false;
    
    if( checkOccupancy( Index2Coord(start_idx) ) )
    {
        /* 起点在障碍物内,沿(start_pt -> end_pt)方向移动直到空闲 */
        do
        {
            start_pt = (start_pt - end_pt).normalized() * step_size_ + start_pt;
            if ( ! Coord2Index(start_pt, start_idx) )
                return false;
        } while ( checkOccupancy( Index2Coord(start_idx) ) );
    }
    
    if( checkOccupancy( Index2Coord(end_idx) ) )
    {
        /* 终点在障碍物内,沿(end_pt -> start_pt)方向移动直到空闲 */
        do
        {
            end_pt = (end_pt - start_pt).normalized() * step_size_ + end_pt;
            if ( ! Coord2Index(end_pt, end_idx) )
                return false;
        } while ( checkOccupancy( Index2Coord(end_idx) ) );
    }

    return true;
}

bool AStar::AstarSearch( const double step_size, Vector3d start_pt, Vector3d end_pt)
{   
    /**
     * 【A*主搜索算法】
     * 
     * 算法步骤:
     * 1. 初始化:起点加入开放集(优先级队列)
     * 2. 主循环:从开放集中取出fScore最小的节点
     *    - 如果是终点,搜索成功
     *    - 否则扩展到26邻域(3x3x3范围,排除中心)
     * 3. 终止条件:开放集为空(搜索失败)或找到终点(搜索成功)
     * 
     * 关键优化:
     * - 使用rounds计数器避免每次搜索前重置节点状态
     * - 设置0.2秒时间限制,防止搜索时间过长
     * - 使用对角线启发式结合tie_breaker提高搜索效率
     */
    ros::Time time_1 = ros::Time::now(); 
    ++rounds_;  // 递增搜索轮次,用于区分当次搜索和其他次搜索

    step_size_ = step_size;
    inv_step_size_ = 1 / step_size;
    center_ = (start_pt + end_pt) / 2;  // 以起点和终点的中点作为地图中心

    Vector3i start_idx, end_idx;
    if ( !ConvertToIndexAndAdjustStartEndPoints(start_pt, end_pt, start_idx, end_idx) )
    {
        ROS_ERROR("Unable to handle the initial or end point, force return!");
        return false;
    }

    /* 获取起点和终点的网格节点指针 */
    GridNodePtr startPtr = GridNodeMap_[start_idx(0)][start_idx(1)][start_idx(2)];
    GridNodePtr endPtr = GridNodeMap_[end_idx(0)][end_idx(1)][end_idx(2)];

    /* 清空开放集(通过交换一个空队列来实现) */
    std::priority_queue< GridNodePtr, std::vector<GridNodePtr>, NodeComparator > empty;
    openSet_.swap(empty);

    GridNodePtr neighborPtr = NULL;
    GridNodePtr current = NULL;

    /* 初始化起点节点 */
    startPtr->index = start_idx;
    startPtr->rounds = rounds_;
    startPtr->gScore = 0;                          // 起点到自己的代价为0
    startPtr->fScore = getHeu(startPtr, endPtr);   // fScore = h(起点,终点)
    startPtr->state = GridNode::OPENSET;           // 标记为开放集
    startPtr->cameFrom = NULL;                     // 起点没有前驱
    openSet_.push( startPtr );                     // 加入优先队列

    endPtr->index = end_idx;

    double tentative_gScore;

    int num_iter = 0;
    while ( !openSet_.empty() )
    { 
        num_iter ++;
        
        /* 从开放集中取出fScore最小的节点 */
        current = openSet_.top();
        openSet_.pop();

        /* 检查是否到达终点 */
        if(current->index(0) == endPtr->index(0)
        && current->index(1) == endPtr->index(1)
        && current->index(2) == endPtr->index(2) )
        {
            /* 找到终点,回溯路径 */
            gridPath_ = retrievePath(current);
            return true;
        }
        
        /* 将当前节点移入关闭集(不再扩展) */
        current -> state = GridNode::CLOSEDSET;

        /* 遍历3x3x3的26邻域(排除自身中心点) */
        for(int dx = -1; dx <= 1; dx++)
            for(int dy = -1; dy <= 1; dy++)
                for(int dz = -1; dz <= 1; dz++)
                {
                    if(dx == 0 && dy == 0 && dz ==0)
                        continue; 

                    Vector3i neighborIdx;
                    neighborIdx(0) = (current -> index)(0) + dx;
                    neighborIdx(1) = (current -> index)(1) + dy;
                    neighborIdx(2) = (current -> index)(2) + dz;

                    /* 跳过边界外的节点(保留一层边界) */
                    if(    neighborIdx(0) < 1 || neighborIdx(0) >= POOL_SIZE_(0) - 1
                        || neighborIdx(1) < 1 || neighborIdx(1) >= POOL_SIZE_(1) - 1
                        || neighborIdx(2) < 1 || neighborIdx(2) >= POOL_SIZE_(2) - 1){
                        continue;
                    }

                    neighborPtr = GridNodeMap_[neighborIdx(0)][neighborIdx(1)][neighborIdx(2)];
                    neighborPtr->index = neighborIdx;

                    /* 检查是否在当前轮次已被探索过 */
                    bool flag_explored = neighborPtr->rounds == rounds_;

                    if( flag_explored && neighborPtr->state == GridNode::CLOSEDSET ){
                        continue; // 在关闭集中,跳过
                    }

                    neighborPtr->rounds = rounds_;

                    if( checkOccupancy( Index2Coord( neighborPtr->index ) ) ){
                        continue; // 被障碍物占据,跳过
                    }

                    /* 计算从起点到邻居节点的实际移动代价 */
                    double static_cost = sqrt(dx * dx + dy * dy + dz * dz);
                    tentative_gScore = current -> gScore + static_cost; 

                    if( ! flag_explored ){
                        /* 发现新节点,加入开放集 */
                        neighborPtr -> state     = GridNode::OPENSET;
                        neighborPtr -> cameFrom  = current;
                        neighborPtr -> gScore    = tentative_gScore;
                        neighborPtr -> fScore    = tentative_gScore + getHeu(neighborPtr, endPtr); 
                        openSet_.push( neighborPtr );
                    }
                    else if(tentative_gScore < neighborPtr-> gScore){
                        /* 节点已在开放集中且通过当前节点路径更优,更新代价和前驱 */
                        neighborPtr -> cameFrom = current;
                        neighborPtr -> gScore = tentative_gScore;
                        neighborPtr -> fScore = tentative_gScore + getHeu(neighborPtr, endPtr); 
                    }   
                }
        
        /* 时间限制:超过0.2秒则判定搜索失败 */
        ros::Time time_2 = ros::Time::now();
        if((time_2 - time_1).toSec() > 0.2){
            ROS_WARN("Failed in A star path searching !!! 0.2 seconds time limit exceeded." );
            return false;
        }
    }

    /* 开放集耗尽但未找到终点,搜索失败 */
    ros::Time time_2 = ros::Time::now();
    if((time_2 - time_1).toSec() > 0.1)
        ROS_WARN("Time consume in A star path finding is %.3fs, iter=%d", (time_2 - time_1).toSec(), num_iter );

    return false;
}

vector<Vector3d> AStar::getPath()
{   
    /**
     * 将A*搜索得到的节点指针序列转换为世界坐标序列
     * 并反转为从起点到终点的顺序
     */
    vector<Vector3d> path;

    for(auto ptr: gridPath_)
        path.push_back( Index2Coord( ptr->index ) );

    reverse(path.begin(), path.end());
    return path;
}