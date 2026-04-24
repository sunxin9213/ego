#ifndef RAYCAST_H_
#define RAYCAST_H_

#include <Eigen/Eigen>
#include <vector>

/**
 * @brief 符号函数
 * @param x 输入值
 * @return 符号值(-1, 0, 或 1)
 */
double signum(double x);

/**
 * @brief 取模运算
 * @param value 被除数
 * @param modulus 除数
 * @return 余数
 */
double mod(double value, double modulus);

/**
 * @brief 计算射线与网格边界相交的时间
 * @param s 起始位置
 * @param ds 方向增量
 * @return 相交时间
 */
double intbound(double s, double ds);

/**
 * @brief 射线投射函数(数组输出版本)
 *
 * 使用3D DDA算法计算射线经过的所有网格点。
 * @param start 射线起点
 * @param end 射线终点
 * @param min 地图最小边界
 * @param max 地图最大边界
 * @param output_points_cnt 输出的点数
 * @param output 输出数组
 */
void Raycast(const Eigen::Vector3d& start, const Eigen::Vector3d& end, const Eigen::Vector3d& min,
             const Eigen::Vector3d& max, int& output_points_cnt, Eigen::Vector3d* output);

/**
 * @brief 射线投射函数(向量输出版本)
 *
 * 使用3D DDA算法计算射线经过的所有网格点。
 * @param start 射线起点
 * @param end 射线终点
 * @param min 地图最小边界
 * @param max 地图最大边界
 * @param output 输出的点向量
 */
void Raycast(const Eigen::Vector3d& start, const Eigen::Vector3d& end, const Eigen::Vector3d& min,
             const Eigen::Vector3d& max, std::vector<Eigen::Vector3d>* output);

/**
 * @brief 射线投射器类
 *
 * 实现3D数字微分分析器(DDA)算法,用于光线追踪和射线投射。
 * 该类支持逐步遍历射线经过的所有体素,可用于地图构建中的占用率更新。
 */
class RayCaster {
private:
  /* data */
  // 射线起点和终点
  Eigen::Vector3d start_;
  Eigen::Vector3d end_;
  // 射线方向
  Eigen::Vector3d direction_;
  // 地图边界
  Eigen::Vector3d min_;
  Eigen::Vector3d max_;
  // 当前网格坐标
  int x_;
  int y_;
  int z_;
  // 终点网格坐标
  int endX_;
  int endY_;
  int endZ_;
  // 最大距离
  double maxDist_;
  // 方向增量
  double dx_;
  double dy_;
  double dz_;
  // 每步方向
  int stepX_;
  int stepY_;
  int stepZ_;
  // 到下一个边界的时间
  double tMaxX_;
  double tMaxY_;
  double tMaxZ_;
  // 到下一个网格的时间增量
  double tDeltaX_;
  double tDeltaY_;
  double tDeltaZ_;
  // 当前距离
  double dist_;

  // 步数计数
  int step_num_;

public:
  RayCaster(/* args */) {
  }
  ~RayCaster() {
  }

  /**
   * @brief 设置射线输入
   * @param start 射线起点
   * @param end 射线终点
   * @return 是否设置成功
   */
  bool setInput(const Eigen::Vector3d& start,
                const Eigen::Vector3d& end /* , const Eigen::Vector3d& min,
                const Eigen::Vector3d& max */);

  /**
   * @brief 射线步进
   * @param ray_pt 输出的射线点
   * @return 是否还有下一个点
   */
  bool step(Eigen::Vector3d& ray_pt);
};

#endif  // RAYCAST_H_