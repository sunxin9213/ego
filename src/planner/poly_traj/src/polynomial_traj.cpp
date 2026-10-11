#include <iostream>
#include <poly_traj/polynomial_traj.h>

/**
 * @brief 最小加加速度(Snap)轨迹生成 - 闭式解法
 * 
 * =====================================================
 * 【算法原理: 最小Snap轨迹优化】
 * =====================================================
 * 
 * 1. 问题定义:
 *   对于无人机轨迹规划,我们希望生成一条经过指定路径点,
 *   且使加加速度(snap,位置的四阶导数)的平方积分最小的轨迹。
 *   
 *   代价函数: J = ∫(d⁴p/dt⁴)²dt = ∫(snap)²dt
 *   最小化snap可以同时减小jerk和加速度变化,使轨迹极为平滑。
 * 
 * 2. 轨迹参数化:
 *   每段轨迹用一个5次多项式表示(6个系数):
 *   p(t) = c₅t⁵ + c₄t⁴ + c₃t³ + c₂t² + c₁t + c₀
 *   
 *   使用5次多项式是因为:
 *   - 需要指定位置、速度、加速度边界(共6个约束)
 *   - 5次多项式恰好有6个系数,刚好满足约束
 *   - 且snap(4阶导)最小化在5次多项式下是有物理意义的
 * 
 * 3. 约束条件:
 *   每段轨迹在起点和终点需要满足:
 *   - 位置约束: p(0) = p_start, p(T) = p_end
 *   - 速度约束: v(0) = v_start, p(T) = p_end  
 *   - 加速度约束: a(0) = a_start, a(T) = a_end
 *   - 中间路径点: C²连续(位置、速度、加速度连续)
 * 
 * 4. 闭式解法:
 *   传统方法直接求解带约束的QP问题,但数值不稳定。
 *   闭式解法通过变量分离将优化变量分为:
 *   - 固定变量(fixed): 位置、速度、加速度等已知约束
 *   - 自由变量(free): 中间点的未知导数
 *   
 *   通过构造选择矩阵C,将QP问题转化为只对自由变量的无约束优化,
 *   得到闭式解: Dp* = -Rpp⁻¹ * Rfpᵀ * Df
 * 
 * 5. 算法步骤:
 *   Step 1: 构造导数向量D(包含所有已知的位置/速度/加速度约束)
 *   Step 2: 构造映射矩阵A(将多项式系数映射到端点导数)
 *   Step 3: 构造选择矩阵C(分离固定/自由变量)
 *   Step 4: 构造最小Snap矩阵Q(目标函数的Hessian矩阵)
 *   Step 5: 计算R矩阵 = C*Aᵀ⁻¹*Q*A⁻¹*Cᵀ
 *   Step 6: 分解R为Rff,Rfp,Rpf,Rpp
 *   Step 7: 求解自由变量: Dp = -Rpp⁻¹*Rfpᵀ*Df
 *   Step 8: 恢复多项式系数: P = A⁻¹*Cᵀ*D
 * 
 * @param Pos 路径点矩阵 (3行 × (seg_num+1)列)
 *            每列是一个路径点,包含起点和所有中间点
 * @param start_vel 起点速度 (3D向量)
 * @param end_vel 终点速度 (3D向量)
 * @param start_acc 起点加速度 (3D向量)
 * @param end_acc 终点加速度 (3D向量)
 * @param Time 每段轨迹的时间向量 (seg_num个时间长度)
 * @return PolynomialTraj 生成的多项式轨迹对象
 */
PolynomialTraj minSnapTraj(const Eigen::MatrixXd& Pos, const Eigen::Vector3d& start_vel,
                           const Eigen::Vector3d& end_vel, const Eigen::Vector3d& start_acc,
                           const Eigen::Vector3d& end_acc, const Eigen::VectorXd& Time) {
  int seg_num = Time.size();  // 轨迹段数
  // 每段6个系数(5次多项式),xyz三轴,共3*6个系数/段
  Eigen::MatrixXd poly_coeff(seg_num, 3 * 6);
  Eigen::VectorXd Px(6 * seg_num), Py(6 * seg_num), Pz(6 * seg_num);

  int num_f, num_p;  // num_f: 固定变量数, num_p: 自由变量数
  int num_d;         // 所有段的导数量

  // 阶乘计算Lambda函数
  const static auto Factorial = [](int x) {
    int fac = 1;
    for (int i = x; i > 0; i--)
      fac = fac * i;
    return fac;
  };


  // =====================================================
  // Step 1: 构造端点导数向量D
  // 每段在起点(索引偏移0)和终点(索引偏移1)有6个导数值:
  //   索引0: 位置 p(0)
  //   索引1: 位置 p(T) 
  //   索引2: 速度 v(0) [仅第一段]
  //   索引3: 速度 v(T) [仅最后一段]
  //   索引4: 加速度 a(0) [仅第一段]
  //   索引5: 加速度 a(T) [仅最后一段]
  // =====================================================
  Eigen::VectorXd Dx = Eigen::VectorXd::Zero(seg_num * 6);
  Eigen::VectorXd Dy = Eigen::VectorXd::Zero(seg_num * 6);
  Eigen::VectorXd Dz = Eigen::VectorXd::Zero(seg_num * 6);

  for (int k = 0; k < seg_num; k++) {
    // 位置约束: 每段起点和终点的位置
    Dx(k * 6)     = Pos(0, k);      // 起点x
    Dx(k * 6 + 1) = Pos(0, k + 1);  // 终点x
    Dy(k * 6)     = Pos(1, k);      // 起点y
    Dy(k * 6 + 1) = Pos(1, k + 1);  // 终点y
    Dz(k * 6)     = Pos(2, k);      // 起点z
    Dz(k * 6 + 1) = Pos(2, k + 1);  // 终点z

    if (k == 0) {
      // 第一段起点: 用户指定的速度和加速度
      Dx(k * 6 + 2) = start_vel(0);
      Dy(k * 6 + 2) = start_vel(1);
      Dz(k * 6 + 2) = start_vel(2);

      Dx(k * 6 + 4) = start_acc(0);
      Dy(k * 6 + 4) = start_acc(1);
      Dz(k * 6 + 4) = start_acc(2);
    } else if (k == seg_num - 1) {
      // 最后一段终点: 用户指定的速度和加速度
      Dx(k * 6 + 3) = end_vel(0);
      Dy(k * 6 + 3) = end_vel(1);
      Dz(k * 6 + 3) = end_vel(2);

      Dx(k * 6 + 5) = end_acc(0);
      Dy(k * 6 + 5) = end_acc(1);
      Dz(k * 6 + 5) = end_acc(2);
    }
    // 注意: 中间段的中间点速度和加速度是自由变量,暂为0
  }

  // =====================================================
  // Step 2: 构造映射矩阵A
  // 将多项式系数映射到端点导数值
  // 对于每段轨迹 p(t) = Σ c_j * t^j, 其导数为:
  //   p(t)的i阶导 = Σ (j!/(j-i)!) * c_j * t^(j-i)
  //   
  // 在t=0处: p^(i)(0) = i! * c_i
  // 在t=T处: p^(i)(T) = Σ (j!/(j-i)!) * c_j * T^(j-i)
  //   
  // 所以映射矩阵A的每一行对应一个约束,每一列对应一个系数:
  //   偶数行(2i): t=0处的i阶导 → Ab(2i, i) = i!
  //   奇数行(2i+1): t=T处的i阶导 → Ab(2i+1, j) = (j!/(j-i)!) * T^(j-i), j>=i
  // =====================================================
  Eigen::MatrixXd Ab;
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(seg_num * 6, seg_num * 6);

  for (int k = 0; k < seg_num; k++) {
    Ab = Eigen::MatrixXd::Zero(6, 6);
    for (int i = 0; i < 3; i++) {  // i=0:位置, i=1:速度, i=2:加速度
      // t=0处: p^(i)(0) = i! * c_i
      Ab(2 * i, i) = Factorial(i);
      // t=T处: p^(i)(T) = Σ (j!/(j-i)!) * c_j * T^(j-i)
      for (int j = i; j < 6; j++)
        Ab(2 * i + 1, j) = Factorial(j) / Factorial(j - i) * pow(Time(k), j - i);
    }
    A.block(k * 6, k * 6, 6, 6) = Ab;
  }

  // =====================================================
  // Step 3: 构造选择矩阵C (变量分离的核心)
  // 
  // 【变量分离思想】
  // 在导数向量D中,有些量是已知的(固定变量),有些是未知的(自由变量):
  // - 固定变量(num_f = 2*seg_num+4):
  //   * 所有段的端点位置: 2*seg_num个(每段起点+终点)
  //   * 起点速度、终点速度: 2个
  //   * 起点加速度、终点加速度: 2个
  //   * 合计: 2*seg_num + 4
  // - 自由变量(num_p = 2*seg_num-2):
  //   * 中间点的速度和加速度(位置是连续的,但导数未知)
  //   * 每个中间点有速度和加速度2个自由量
  //   * (seg_num-1)个中间点 → 2*(seg_num-1) = 2*seg_num-2
  // 
  // 选择矩阵C的作用:
  //   D_ordered = C * D_original
  // 其中D_ordered将固定变量排列在前面,自由变量排列在后面
  // 从而将优化问题块对角化
  // =====================================================
  Eigen::MatrixXd Ct, C;

  num_f = 2 * seg_num + 4;  // 固定变量数量
  num_p = 2 * seg_num - 2;  // 自由变量数量  
  num_d = 6 * seg_num;      // 总导数量
  Ct = Eigen::MatrixXd::Zero(num_d, num_f + num_p);
  
  // --- 第一段起点(固定变量) ---
  Ct(0, 0) = 1;  // 起点位置p→第0个固定变量
  Ct(2, 1) = 1;  // 起点速度v0→第1个固定变量
  Ct(4, 2) = 1;  // 起点加速度a0→第2个固定变量

  // 第一段终点也是第二段起点(连续性)
  Ct(1, 3) = 1;                              // 终点位置p→第3个固定变量
  Ct(3, 2 * seg_num + 4) = 1;               // 中间点速度→自由变量[0]
  Ct(5, 2 * seg_num + 5) = 1;               // 中间点加速度→自由变量[1]

  // --- 最后一段终点(固定变量) ---
  Ct(6 * (seg_num - 1) + 0, 2 * seg_num + 0) = 1;  // 终点位置p
  Ct(6 * (seg_num - 1) + 1, 2 * seg_num + 1) = 1;  // 终点位置p
  Ct(6 * (seg_num - 1) + 2, 4 * seg_num + 0) = 1;  // 终点速度→固定变量
  Ct(6 * (seg_num - 1) + 3, 2 * seg_num + 2) = 1;  // 终点位置p
  Ct(6 * (seg_num - 1) + 4, 4 * seg_num + 1) = 1;  // 终点加速度→固定变量
  Ct(6 * (seg_num - 1) + 5, 2 * seg_num + 3) = 1;  // 终点位置p

  // --- 中间段(j=2,...,seg_num-1): 连续性约束 ---
  for (int j = 2; j < seg_num; j++) {
    // 中间段起点位置=上一段终点位置(连续)
    Ct(6 * (j - 1) + 0, 2 + 2 * (j - 1) + 0) = 1;
    Ct(6 * (j - 1) + 1, 2 + 2 * (j - 1) + 1) = 1;
    // 中间段起点速度(自由变量)
    Ct(6 * (j - 1) + 2, 2 * seg_num + 4 + 2 * (j - 2) + 0) = 1;
    Ct(6 * (j - 1) + 3, 2 * seg_num + 4 + 2 * (j - 1) + 0) = 1;
    // 中间段起点加速度(自由变量)
    Ct(6 * (j - 1) + 4, 2 * seg_num + 4 + 2 * (j - 2) + 1) = 1;
    Ct(6 * (j - 1) + 5, 2 * seg_num + 4 + 2 * (j - 1) + 1) = 1;
  }

  C = Ct.transpose();  // C是转置选择矩阵

  // 重排后的导数向量: 固定变量在前,自由变量在后
  Eigen::VectorXd Dx1 = C * Dx;
  Eigen::VectorXd Dy1 = C * Dy;
  Eigen::VectorXd Dz1 = C * Dz;

  // =====================================================
  // Step 4: 构造最小Snap矩阵Q(Hessian矩阵)
  // 
  // 代价函数: J = ∫(d⁴p/dt⁴)²dt
  // 
  // 对于5次多项式 p(t) = Σ_{i=0}^{5} c_i t^i
  //   snap(t) = d⁴p/dt⁴ = Σ_{i=4}^{5} i(i-1)(i-2)(i-3) c_i t^(i-4)
  //           = 24*c₄ + 120*c₅*t
  // 
  // snap²的积分:
  //   ∫(24*c₄ + 120*c₅*t)²dt = ...
  // 
  // 更一般的形式,Q矩阵的元素为:
  //   Q(i,j) = ∫(i(i-1)(i-2)(i-3))*(j(j-1)(j-2)(j-3))*t^(i+j-8)dt
  //          = i(i-1)(i-2)(i-3)*j(j-1)(j-2)(j-3)*T^(i+j-7)/(i+j-7)
  // 
  // 注意: 对于5次多项式,只有i,j≥4的项才产生非零snap
  // 所以i和j从3开始(因为代码中索引0~5对应次数0~5)
  // 但Q矩阵只需要系数i≥3(即代码索引3~5,对应4次和5次项)的交互
  // 
  // 代码中的索引对应关系:
  //   i=3 → 4次项系数c₄, snap贡献: 4*3*2*1=24
  //   i=4 → 5次项系数c₅, snap贡献: 5*4*3*2=120  
  //   i=5 → ... (对5次多项式,没有6次项)
  // 
  // 所以实际上i,j从3到4(但代码用了3到5,多算了一些0)
  // =====================================================
  Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(seg_num * 6, seg_num * 6);

  for (int k = 0; k < seg_num; k++) {
    for (int i = 3; i < 6; i++) {
      for (int j = 3; j < 6; j++) {
        // Q[i][j] = i*(i-1)*(i-2) * j*(j-1)*(j-2) * T^(i+j-5) / (i+j-5)
        // 注意: 这里i,j是代码索引对应多项式的(i+1)次项
        // 修正后的实际意义:
        //   i*(i-1)*(i-2) = (i+1)!/(i-2)! 即(i+1)阶多项式的4阶导系数
        Q(k * 6 + i, k * 6 + j) =
            i * (i - 1) * (i - 2) * j * (j - 1) * (j - 2) / (i + j - 5) * pow(Time(k), (i + j - 5));
      }
    }
  }

  // =====================================================
  // Step 5-7: 闭式求解
  // 
  // 【核心推导】
  // 原始问题: min cᵀQc, s.t. Ac = D (约束方程)
  // 
  // 变换到导数空间:
  //   代入 c = A⁻¹D:
  //   min Dᵀ(A⁻¹)ᵀQA⁻¹D, s.t. 部分D已知
  // 
  // 变量分离:
  //   令 R = C * A⁻ᵀ * Q * A⁻¹ * Cᵀ
  //   则代价 = [Dfᵀ, Dpᵀ] * [Rff, Rfp; Rpf, Rpp] * [Df; Dp]
  //   其中Df是固定变量,Dp是自由变量
  // 
  // 无约束最优条件: ∂J/∂Dp = 0
  //   → 2*Rpf*Df + 2*Rpp*Dp = 0
  //   → Dp* = -Rpp⁻¹ * Rpf * Df = -Rpp⁻¹ * Rfpᵀ * Df
  // 
  // 这就是闭式解!只需一次矩阵求逆和乘法即可得到所有自由变量。
  // =====================================================
  Eigen::MatrixXd R = C * A.transpose().inverse() * Q * A.inverse() * Ct;

  // 提取固定变量
  Eigen::VectorXd Dxf(2 * seg_num + 4), Dyf(2 * seg_num + 4), Dzf(2 * seg_num + 4);
  Dxf = Dx1.segment(0, 2 * seg_num + 4);
  Dyf = Dy1.segment(0, 2 * seg_num + 4);
  Dzf = Dz1.segment(0, 2 * seg_num + 4);

  // 分块R矩阵
  Eigen::MatrixXd Rff(2 * seg_num + 4, 2 * seg_num + 4);
  Eigen::MatrixXd Rfp(2 * seg_num + 4, 2 * seg_num - 2);
  Eigen::MatrixXd Rpf(2 * seg_num - 2, 2 * seg_num + 4);
  Eigen::MatrixXd Rpp(2 * seg_num - 2, 2 * seg_num - 2);

  Rff = R.block(0, 0, 2 * seg_num + 4, 2 * seg_num + 4);
  Rfp = R.block(0, 2 * seg_num + 4, 2 * seg_num + 4, 2 * seg_num - 2);
  Rpf = R.block(2 * seg_num + 4, 0, 2 * seg_num - 2, 2 * seg_num + 4);
  Rpp = R.block(2 * seg_num + 4, 2 * seg_num + 4, 2 * seg_num - 2, 2 * seg_num - 2);

  // --- 闭式求解自由变量 ---
  Eigen::VectorXd Dxp(2 * seg_num - 2), Dyp(2 * seg_num - 2), Dzp(2 * seg_num - 2);
  Dxp = -(Rpp.inverse() * Rfp.transpose()) * Dxf;
  Dyp = -(Rpp.inverse() * Rfp.transpose()) * Dyf;
  Dzp = -(Rpp.inverse() * Rfp.transpose()) * Dzf;

  // 回填自由变量到完整导数向量
  Dx1.segment(2 * seg_num + 4, 2 * seg_num - 2) = Dxp;
  Dy1.segment(2 * seg_num + 4, 2 * seg_num - 2) = Dyp;
  Dz1.segment(2 * seg_num + 4, 2 * seg_num - 2) = Dzp;

  // =====================================================
  // Step 8: 恢复多项式系数
  //   c = A⁻¹ * Cᵀ * D
  // 将导数空间的值映射回多项式系数空间
  // =====================================================
  Px = (A.inverse() * Ct) * Dx1;
  Py = (A.inverse() * Ct) * Dy1;
  Pz = (A.inverse() * Ct) * Dz1;

  // 整理每段轨迹的系数(6个系数×3轴)
  for (int i = 0; i < seg_num; i++) {
    poly_coeff.block(i, 0, 1, 6)   = Px.segment(i * 6, 6).transpose();
    poly_coeff.block(i, 6, 1, 6)   = Py.segment(i * 6, 6).transpose();
    poly_coeff.block(i, 12, 1, 6)  = Pz.segment(i * 6, 6).transpose();
  }

  // =====================================================
  // 构造PolynomialTraj对象
  // 
  // 【注意系数顺序】
  // 上面求解的系数是: [c₅, c₄, c₃, c₂, c₁, c₀] (从高次到低次)
  // 但PolynomialTraj的evaluate函数需要: [c₀, c₁, c₂, c₃, c₄, c₅] (从低次到高次)
  // 所以需要用reverse反转顺序
  // =====================================================
  PolynomialTraj poly_traj;
  for (int i = 0; i < poly_coeff.rows(); ++i) {
    vector<double> cx(6), cy(6), cz(6);
    for (int j = 0; j < 6; ++j) {
      cx[j] = poly_coeff(i, j), cy[j] = poly_coeff(i, j + 6), cz[j] = poly_coeff(i, j + 12);
    }
    // 反转系数顺序: [c₅,c₄,c₃,c₂,c₁,c₀] → [c₀,c₁,c₂,c₃,c₄,c₅]
    reverse(cx.begin(), cx.end());
    reverse(cy.begin(), cy.end());
    reverse(cz.begin(), cz.end());
    double ts = Time(i);
    poly_traj.addSegment(cx, cy, cz, ts);
  }

  return poly_traj;
}

/**
 * @brief 快速直线轨迹生成 (4阶/7段, jerk连续)
 * 
 * =====================================================
 * 【算法原理: 时间最优梯形速度剖面】
 * =====================================================
 * 
 * 本函数生成从起点到终点的直线轨迹,满足最大速度、加速度和加加速度约束。
 * 使用7段分段多项式构造一个jerk连续的"S曲线"速度剖面:
 * 
 *   jerk      ████      ████
 *   ↑    ████    ████    ████
 *   +────┴────┴────┴────┴────→ t
 *   ↓ ████                  ████
 *   
 *   a(t)   ┌────┐           ┌────
 *   ↑  ┌──┘    └──┐   ┌───┘
 *   +──┴─────────┴───┴────────→ t
 *   
 *   v(t)     ┌────────┐
 *   ↑  ┌─────┘        └─────
 *   +──┴────────────────────→ t
 *   
 *   7段分段:
 *   1. 加jerk阶段: jerk = +max_jerk, 加速度从0到+max_acc
 *   2. 恒加速阶段: jerk = 0, 加速度保持+max_acc, 速度线性增加
 *   3. 减jerk阶段: jerk = -max_jerk, 加速度从+max_acc降到0
 *   4. 匀速阶段: jerk = 0, 加速度=0, 速度=+max_vel
 *   5. 减jerk阶段: jerk = -max_jerk, 加速度从0到-max_acc
 *   6. 恒减速阶段: jerk = 0, 加速度保持-max_acc
 *   7. 加jerk阶段: jerk = +max_jerk, 加速度从-max_acc回到0
 * 
 * 【参数计算】
 * t1 = t3 = t5 = t7 = max_acc / max_jerk  (加减jerk时间)
 * t2 = t6 = max_vel / max_acc - t1          (恒加/减速时间)  
 * t4 = total_dist / max_vel - 2*t1 - t2     (匀速时间)
 * 
 * 每段轨迹是一个3次多项式(系数[c₀,c₁,c₂,c₃,c₄,c₅]):
 *   p(t) = c₅ + c₄*t + c₃*t²/2 + c₂*t³/6
 *   其中c₅=初始位置, c₄=初始速度, c₃=初始加速度, c₂=初始jerk
 * 
 * @param start 起点位置 (3D)
 * @param end 终点位置 (3D)
 * @param max_vel 最大速度 (标量)
 * @param max_acc 最大加速度 (标量)
 * @param max_jerk 最大加加速度 (标量)
 * @return PolynomialTraj 7段多项式轨迹
 */
PolynomialTraj fastLine4deg(Eigen::Vector3d start, Eigen::Vector3d end, double max_vel, double max_acc,
                            double max_jerk) {
  Eigen::Vector3d disp = end - start;
  double len = disp.norm();
  Eigen::Vector3d dir = disp.normalized();

  // 计算各轴缩放系数(处理各轴分量不同的情况)
  int max_id = -1;
  double max_dist = -1.0;
  for (int i = 0; i < 3; ++i) {
    if (fabs(disp(i)) > max_dist) {
      max_dist = disp(i);
      max_id = i;
    }
  }
  Eigen::Vector3d scale_vec = disp / max_dist;

  PolynomialTraj poly_traj;
  vector<double> cx(6), cy(6), cz(6), zero(6);
  for (int i = 0; i < 6; ++i)
    zero[i] = 0.0;

  Eigen::Vector3d p0 = start;
  Eigen::Vector3d v0 = Eigen::Vector3d::Zero();
  Eigen::Vector3d a0 = Eigen::Vector3d::Zero();
  Eigen::Vector3d j0 = max_jerk * scale_vec;

  // ====== 第1段: 加jerk阶段(t1时间) ======
  // 初始条件: p=p0, v=0, a=0, j=+max_jerk
  // 加速度从0增加到max_acc
  cx = cy = cz = zero;
  double t1 = max_acc / max_jerk;

  cx[5] = p0(0); cy[5] = p0(1); cz[5] = p0(2);  // 位置
  cx[4] = v0(0); cy[4] = v0(1); cz[4] = v0(2);  // 速度=0
  cx[3] = (1 / 2.0) * a0(0); cy[3] = (1 / 2.0) * a0(1); cz[3] = (1 / 2.0) * a0(2);  // 加速度=0
  cx[2] = (1 / 6.0) * j0(0); cy[2] = (1 / 6.0) * j0(1); cz[2] = (1 / 6.0) * j0(2);  // jerk

  poly_traj.addSegment(cx, cy, cz, t1);

  // 传播到第1段结束的状态
  Eigen::Vector3d p1 = p0 + v0 * t1 + (1 / 2.0) * a0 * pow(t1, 2) + (1 / 6.0) * j0 * pow(t1, 3);
  Eigen::Vector3d v1 = v0 + a0 * t1 + (1 / 2.0) * j0 * pow(t1, 2);
  Eigen::Vector3d a1 = a0 + j0 * t1;
  Eigen::Vector3d j1 = Eigen::Vector3d::Zero();  // jerk回到0

  // ====== 第2段: 恒加速阶段(t2时间) ======
  // 初始条件: p=p1, v=v1, a=+max_acc, j=0
  cx = cy = cz = zero;
  double t2 = max_vel / max_acc - t1;  // 保证速度达到max_vel

  cx[5] = p1(0); cy[5] = p1(1); cz[5] = p1(2);
  cx[4] = v1(0); cy[4] = v1(1); cz[4] = v1(2);
  cx[3] = (1 / 2.0) * a1(0); cy[3] = (1 / 2.0) * a1(1); cz[3] = (1 / 2.0) * a1(2);
  cx[2] = (1 / 6.0) * j1(0); cy[2] = (1 / 6.0) * j1(1); cz[2] = (1 / 6.0) * j1(2);

  poly_traj.addSegment(cx, cy, cz, t2);

  Eigen::Vector3d p2 = p1 + v1 * t2 + (1 / 2.0) * a1 * pow(t2, 2) + (1 / 6.0) * j1 * pow(t2, 3);
  Eigen::Vector3d v2 = v1 + a1 * t2 + (1 / 2.0) * j1 * pow(t2, 2);
  Eigen::Vector3d a2 = a1 + j1 * t2;
  Eigen::Vector3d j2 = -max_jerk * scale_vec;  // jerk反向(减jerk)

  // ====== 第3段: 减jerk阶段(t3=t1) ======
  // 初始条件: p=p2, v=v2, a=+max_acc, j=-max_jerk
  // 加速度从+max_acc线性降到0
  cx = cy = cz = zero;
  double t3 = t1;

  cx[5] = p2(0); cy[5] = p2(1); cz[5] = p2(2);
  cx[4] = v2(0); cy[4] = v2(1); cz[4] = v2(2);
  cx[3] = (1 / 2.0) * a2(0); cy[3] = (1 / 2.0) * a2(1); cz[3] = (1 / 2.0) * a2(2);
  cx[2] = (1 / 6.0) * j2(0); cy[2] = (1 / 6.0) * j2(1); cz[2] = (1 / 6.0) * j2(2);

  poly_traj.addSegment(cx, cy, cz, t3);

  Eigen::Vector3d p3 = p2 + v2 * t3 + (1 / 2.0) * a2 * pow(t3, 2) + (1 / 6.0) * j2 * pow(t3, 3);
  Eigen::Vector3d v3 = v2 + a2 * t3 + (1 / 2.0) * j2 * pow(t3, 2);
  Eigen::Vector3d a3 = a2 + j2 * t3;
  Eigen::Vector3d j3 = Eigen::Vector3d::Zero();

  // ====== 第4段: 匀速阶段(t4时间) ======
  // 初始条件: p=p3, v=max_vel, a=0, j=0
  cx = cy = cz = zero;
  double t4 = max_dist / max_vel - 2 * t1 - t2;  // 总时间减去加减速时间

  cx[5] = p3(0); cy[5] = p3(1); cz[5] = p3(2);
  cx[4] = v3(0); cy[4] = v3(1); cz[4] = v3(2);
  cx[3] = (1 / 2.0) * a3(0); cy[3] = (1 / 2.0) * a3(1); cz[3] = (1 / 2.0) * a3(2);
  cx[2] = (1 / 6.0) * j3(0); cy[2] = (1 / 6.0) * j3(1); cz[2] = (1 / 6.0) * j3(2);

  poly_traj.addSegment(cx, cy, cz, t4);

  Eigen::Vector3d p4 = p3 + v3 * t4 + (1 / 2.0) * a3 * pow(t4, 2) + (1 / 6.0) * j3 * pow(t4, 3);
  Eigen::Vector3d v4 = v3 + a3 * t4 + (1 / 2.0) * j3 * pow(t4, 2);
  Eigen::Vector3d a4 = a3 + j3 * t4;
  Eigen::Vector3d j4 = -max_jerk * scale_vec;  // 开始减速

  // ====== 第5-7段: 减速对称阶段(与第1-3段对称) ======
  // 第5段: 减jerk, 加速度从0到-max_acc
  cx = cy = cz = zero;
  double t5 = t1;

  cx[5] = p4(0); cy[5] = p4(1); cz[5] = p4(2);
  cx[4] = v4(0); cy[4] = v4(1); cz[4] = v4(2);
  cx[3] = (1 / 2.0) * a4(0); cy[3] = (1 / 2.0) * a4(1); cz[3] = (1 / 2.0) * a4(2);
  cx[2] = (1 / 6.0) * j4(0); cy[2] = (1 / 6.0) * j4(1); cz[2] = (1 / 6.0) * j4(2);

  poly_traj.addSegment(cx, cy, cz, t5);

  Eigen::Vector3d p5 = p4 + v4 * t5 + (1 / 2.0) * a4 * pow(t5, 2) + (1 / 6.0) * j4 * pow(t5, 3);
  Eigen::Vector3d v5 = v4 + a4 * t5 + (1 / 2.0) * j4 * pow(t5, 2);
  Eigen::Vector3d a5 = a4 + j4 * t5;
  Eigen::Vector3d j5 = Eigen::Vector3d::Zero();

  // 第6段: 恒减速
  cx = cy = cz = zero;
  double t6 = t2;

  cx[5] = p5(0); cy[5] = p5(1); cz[5] = p5(2);
  cx[4] = v5(0); cy[4] = v5(1); cz[4] = v5(2);
  cx[3] = (1 / 2.0) * a5(0); cy[3] = (1 / 2.0) * a5(1); cz[3] = (1 / 2.0) * a5(2);
  cx[2] = (1 / 6.0) * j5(0); cy[2] = (1 / 6.0) * j5(1); cz[2] = (1 / 6.0) * j5(2);

  poly_traj.addSegment(cx, cy, cz, t6);

  Eigen::Vector3d p6 = p5 + v5 * t6 + (1 / 2.0) * a5 * pow(t6, 2) + (1 / 6.0) * j5 * pow(t6, 3);
  Eigen::Vector3d v6 = v5 + a5 * t6 + (1 / 2.0) * j5 * pow(t6, 2);
  Eigen::Vector3d a6 = a5 + j5 * t6;
  Eigen::Vector3d j6 = max_jerk * scale_vec;  // 最后减速到0

  // 第7段: 加jerk回到零加速度和零速度
  cx = cy = cz = zero;
  double t7 = t1;

  cx[5] = p6(0); cy[5] = p6(1); cz[5] = p6(2);
  cx[4] = v6(0); cy[4] = v6(1); cz[4] = v6(2);
  cx[3] = (1 / 2.0) * a6(0); cy[3] = (1 / 2.0) * a6(1); cz[3] = (1 / 2.0) * a6(2);
  cx[2] = (1 / 6.0) * j6(0); cy[2] = (1 / 6.0) * j6(1); cz[2] = (1 / 6.0) * j6(2);

  poly_traj.addSegment(cx, cy, cz, t7);

  return poly_traj;
}

/**
 * @brief 快速直线轨迹生成 (3阶/3段, 梯形速度剖面, 无jerk约束)
 * 
 * =====================================================
 * 【算法原理: 梯形速度剖面】
 * =====================================================
 * 
 * 与fastLine4deg(7段S曲线)不同,本函数使用更简单的3段梯形速度剖面:
 * 
 *   v(t)
 *   ↑      ┌──────────┐
 *   max_vel─┘          └──
 *   +────────────────────→ t
 *   
 *   a(t)
 *   ↑  ┌────┐         ┌────
 *   +──┴────┴─────────┴────┴→ t
 *   max_acc     0     -max_acc
 *   
 *   3段分段:
 *   第1段(加速段): 恒定加速度 +max_acc, 速度从0线性增加到max_vel
 *   第2段(匀速段): 加速度=0, 保持最大速度max_vel
 *   第3段(减速段): 恒定加速度 -max_acc, 速度从max_vel线性降到0
 * 
 * 【与fastLine4deg的区别】
 * 特征          | fastLine4deg (7段) | fastLine3deg (3段)
 * --------------|-------------------|-------------------
 * 阶数          | 4阶(jerk连续)     | 3阶(加速度连续)
 * 段数          | 7                 | 3
 * jerk约束      | 有                | 无
 * 加速度变化    | 平滑(S曲线)       | 阶跃(梯形)
 * 计算开销      | 较大              | 小
 * 适用场景      | 高精度轨迹        | 快速粗略轨迹
 * 
 * 【参数计算】
 * 设总位移为len,则需要:
 *   加速段长度: len_ht = 0.5 * max_acc * t_ht²
 *     其中 t_ht = max_vel / max_acc (加速到最大速度所需时间)
 *   匀速段时间: t_m = (len - 2*len_ht) / max_vel
 *     即减去加减速段的剩余距离除以最大速度
 * 
 * 约束条件(需要保证):
 *   len >= 2 * len_ht = max_vel² / max_acc
 *   即轨迹必须有足够的距离完成加速和减速
 *   如果距离不足,需要降低最大速度或增大加速度
 * 
 * 【系数约定】
 * 多项式 p(t) = c₀ + c₁*t + c₂*t²/2 + c₃*t³/6
 * 使用6个系数的向量[c₀,c₁,c₂,c₃,c₄,c₅],但c₄=c₅=0(3次多项式)
 *   索引5 = c₀ = 初始位置 p(0)
 *   索引4 = c₁ = 初始速度 v(0)
 *   索引3 = c₂/2 = 初始加速度 a(0)/2
 *   索引2 = c₃/6 = 初始jerk/6 = 0 (3次多项式没有jerk项)
 *   索引1,0 = 0 (无更高次项)
 * 
 * @param start 起点位置 (3D向量)
 * @param end 终点位置 (3D向量)
 * @param max_vel 最大速度限制 (标量)
 * @param max_acc 最大加速度限制 (标量, 正值, 函数内部自动取负用于减速)
 * @return PolynomialTraj 3段多项式轨迹对象
 */
PolynomialTraj fastLine3deg(Eigen::Vector3d start, Eigen::Vector3d end, double max_vel, double max_acc) {
  Eigen::Vector3d disp = end - start;      // 位移向量
  double len = disp.norm();                // 总位移距离
  Eigen::Vector3d dir = disp.normalized(); // 位移方向单位向量(用于路径点计算)

  // =====================================================
  // 计算各轴缩放系数
  // 当各轴分量不同时(例如x方向位移大,y方向位移小),
  // 使用最大分量为基准进行缩放,保证各轴同时达到最大速度
  // 例如: disp = [10, 5, 0], 则 scale_vec = [1, 0.5, 0]
  // 这样x轴以全速运动,y轴以半速运动,保证各轴同步到达
  // =====================================================
  double max_dist = std::max(disp(0), disp(1));
  max_dist = std::max(max_dist, disp(2));
  Eigen::Vector3d scale_vec = disp / max_dist;

  // =====================================================
  // 计算各段时间
  // t_ht: 加/减速段持续时间
  //   梯形剖面中,速度从0到max_vel需要的时间
  //   由 v = a*t 得: t = v/a
  // =====================================================
  double t_ht = max_vel / max_acc;

  // len_ht: 加/减速段走过的距离
  //   匀加速运动中: s = 0.5 * a * t²
  //   但这里用 scaled acceleration,所以用scale_vec.norm()修正
  double len_ht = 0.5 * (max_acc * scale_vec.norm()) * t_ht * t_ht;
  
  // t_m: 匀速段时间
  //   总距离减去加减速段距离,除以最大速度
  double t_m = (len - 2 * len_ht) / (max_vel * scale_vec.norm());

  // 加减速段与匀速段的两个连接点(中间航点)
  // mpt1: 加速段结束/匀速段开始
  // mpt2: 匀速段结束/减速段开始
  Eigen::Vector3d mpt1 = start + dir * len_ht;
  Eigen::Vector3d mpt2 = end - dir * len_ht;

  // 3段轨迹的容器
  PolynomialTraj poly_traj;
  vector<double> cx(6), cy(6), cz(6), zero(6);
  for (int i = 0; i < 6; ++i)
    zero[i] = 0.0;

  // ====== 第1段: 加速段(t_ht时间) ======
  // 初始条件: 位置=start, 速度=0, 加速度=+max_acc*scale_vec
  // 多项式系数: [c₀=start, c₁=0, c₂/2=max_acc/2, ...]
  // 在t=t_ht时, 速度达到max_vel*scale_vec, 位置到达mpt1
  cx = cy = cz = zero;

  cx[5] = start(0);  // p0 = start
  cy[5] = start(1);
  cz[5] = start(2);

  cx[4] = 0.0;  // v0 = 0
  cy[4] = 0.0;
  cz[4] = 0.0;

  cx[3] = 0.5 * max_acc * scale_vec(0);  // a0/2 = max_acc/2 * scale
  cy[3] = 0.5 * max_acc * scale_vec(1);
  cz[3] = 0.5 * max_acc * scale_vec(2);

  poly_traj.addSegment(cx, cy, cz, t_ht);

  // ====== 第2段: 匀速段(t_m时间) ======
  // 初始条件: 位置=mpt1, 速度=max_vel*scale_vec, 加速度=0
  // 在t=t_m时, 位置到达mpt2
  cx = cy = cz = zero;

  cx[5] = mpt1(0);  // p0 = mpt1
  cy[5] = mpt1(1);
  cz[5] = mpt1(2);

  cx[4] = max_vel * scale_vec(0);  // v0 = max_vel * scale
  cy[4] = max_vel * scale_vec(1);
  cz[4] = max_vel * scale_vec(2);

  cx[3] = 0.0;  // a0 = 0 (匀速)
  cy[3] = 0.0;
  cz[3] = 0.0;

  poly_traj.addSegment(cx, cy, cz, t_m);

  // ====== 第3段: 减速段(t_ht时间) ======
  // 初始条件: 位置=mpt2, 速度=max_vel*scale_vec, 加速度=-max_acc*scale_vec
  // 注意加速度为负值,表示减速
  // 在t=t_ht时, 速度降到0, 位置到达end
  cx = cy = cz = zero;

  cx[5] = mpt2(0);  // p0 = mpt2
  cy[5] = mpt2(1);
  cz[5] = mpt2(2);

  cx[4] = max_vel * scale_vec(0);  // v0 = max_vel * scale
  cy[4] = max_vel * scale_vec(1);
  cz[4] = max_vel * scale_vec(2);

  cx[3] = -0.5 * max_acc * scale_vec(0);  // a0/2 = -max_acc/2 * scale (负号表示减速)
  cy[3] = -0.5 * max_acc * scale_vec(1);
  cz[3] = -0.5 * max_acc * scale_vec(2);

  poly_traj.addSegment(cx, cy, cz, t_ht);

  return poly_traj;
}

/**
 * @brief 单段多项式轨迹生成 (5次多项式, 满足位置/速度/加速度边界条件)
 * 
 * =====================================================
 * 【算法原理: 带有边界约束的多项式插值】
 * =====================================================
 * 
 * 本函数生成一段5次多项式轨迹,严格满足起点和终点的位置、速度、加速度约束。
 * 共有6个约束条件(起点3个+终点3个),因此使用5次多项式(6个系数)。
 * 
 * 【约束条件】
 *   起点(t=0):  p(0) = start_pt
 *              v(0) = start_vel
 *              a(0) = start_acc
 *   终点(t=T):  p(T) = end_pt
 *              v(T) = end_vel
 *              a(T) = end_acc
 * 
 * 【约束矩阵构造】
 *   设多项式: p(t) = c₅t⁵ + c₄t⁴ + c₃t³ + c₂t² + c₁t + c₀
 *   则:
 *   p(t)  = c₅t⁵ + c₄t⁴ + c₃t³ + c₂t² + c₁t + c₀
 *   v(t)  = 5c₅t⁴ + 4c₄t³ + 3c₃t² + 2c₂t + c₁
 *   a(t)  = 20c₅t³ + 12c₄t² + 6c₃t + 2c₂
 * 
 *   约束矩阵形式: C * c = b
 *   其中C是6×6矩阵, c是系数向量[c₀,c₁,c₂,c₃,c₄,c₅], b是边界值向量
 * 
 *   [1  0  0  0   0   0 ] [c₀]   [start_pt]
 *   [0  1  0  0   0   0 ] [c₁]   [start_vel]
 *   [0  0  2  0   0   0 ] [c₂]   [start_acc]
 *   [1  T  T² T³  T⁴  T⁵] [c₃] = [end_pt  ]
 *   [0  1  2T 3T² 4T³ 5T⁴] [c₄]   [end_vel ]
 *   [0  0  2  6T 12T² 20T³] [c₅]   [end_acc ]
 * 
 * 【求解方法】
 *   使用QR分解(colPivHouseholderQr)求解线性方程组:
 *   c = C.solve(b)
 *   QR分解比常规求逆更稳定,对接近奇异的矩阵也有较好的数值特性
 * 
 * 【应用场景】
 *   在Ego算法中,此函数用于:
 *   1. 生成单段连接轨迹(如紧急情况下的直接连接)
 *   2. 作为复杂轨迹的构建基块
 *   3. 验证轨迹的可行性(已知端点状态,插值检查中间点)
 * 
 * @param start_pt 起点位置 (3D)
 * @param start_vel 起点速度 (3D)
 * @param start_acc 起点加速度 (3D)
 * @param end_pt 终点位置 (3D)
 * @param end_vel 终点速度 (3D)
 * @param end_acc 终点加速度 (3D)
 * @param t 轨迹持续时间 (标量, >0)
 * @return PolynomialTraj 包含单段5次多项式的轨迹对象
 */
PolynomialTraj one_segment_traj_gen( const Eigen::Vector3d& start_pt, const Eigen::Vector3d& start_vel, const Eigen::Vector3d& start_acc, 
                            const Eigen::Vector3d& end_pt, const Eigen::Vector3d& end_vel, const Eigen::Vector3d& end_acc, 
                            double t )
{
  // =====================================================
  // Step 1: 构造约束矩阵C (6×6)
  // 矩阵的每一行对应一个约束条件:
  //   行0: p(0) = c₀ = start_pt
  //   行1: v(0) = c₁ = start_vel
  //   行2: a(0) = 2c₂ = start_acc
  //   行3: p(T) = c₀ + c₁T + c₂T² + c₃T³ + c₄T⁴ + c₅T⁵ = end_pt
  //   行4: v(T) = c₁ + 2c₂T + 3c₃T² + 4c₄T³ + 5c₅T⁴ = end_vel
  //   行5: a(T) = 2c₂ + 6c₃T + 12c₄T² + 20c₅T³ = end_acc
  // =====================================================
  Eigen::MatrixXd C = Eigen::MatrixXd::Zero(6,6), Crow(1,6);
  Eigen::VectorXd Bx(6), By(6), Bz(6);

  // --- 起点约束(t=0) ---
  C(0,5) = 1;              // p(0) = c₀ → 索引5对应c₀
  C(1,4) = 1;              // v(0) = c₁ → 索引4对应c₁
  C(2,3) = 2;              // a(0) = 2*c₂ → 索引3对应c₂

  // --- 终点约束(t=T) ---
  // p(T) = c₀ + c₁T + c₂T² + c₃T³ + c₄T⁴ + c₅T⁵
  Crow << pow(t,5), pow(t,4), pow(t,3), pow(t,2), t, 1;
  C.row(3) = Crow;
  
  // v(T) = c₁ + 2c₂T + 3c₃T² + 4c₄T³ + 5c₅T⁴
  Crow << 5*pow(t,4), 4*pow(t,3), 3*pow(t,2), 2*t, 1, 0;
  C.row(4) = Crow;
  
  // a(T) = 2c₂ + 6c₃T + 12c₄T² + 20c₅T³
  Crow << 20*pow(t,3), 12*pow(t,2), 6*t, 2, 0, 0;
  C.row(5) = Crow;

  // =====================================================
  // Step 2: 构造边界值向量B
  // 包含起点和终点的6个约束值
  // =====================================================
  Bx << start_pt(0), start_vel(0), start_acc(0), end_pt(0), end_vel(0), end_acc(0);
  By << start_pt(1), start_vel(1), start_acc(1), end_pt(1), end_vel(1), end_acc(1);
  Bz << start_pt(2), start_vel(2), start_acc(2), end_pt(2), end_vel(2), end_acc(2);

  // =====================================================
  // Step 3: 求解多项式系数
  // 使用列主元HouseholderQR分解求解 C*c = B
  // QR分解的优势:
  //   - 数值稳定性好,对病态矩阵有较好的鲁棒性
  //   - colPivHouseholderQr对秩亏矩阵也能给出最小范数解
  //   - 比LU分解更稳定,比SVD分解更快
  // =====================================================
  Eigen::VectorXd Cofx = C.colPivHouseholderQr().solve(Bx);
  Eigen::VectorXd Cofy = C.colPivHouseholderQr().solve(By);
  Eigen::VectorXd Cofz = C.colPivHouseholderQr().solve(Bz);

  // 将Eigen向量转为std::vector并构造轨迹
  vector<double> cx(6), cy(6), cz(6);
  for ( int i=0; i<6; i++ )
  {
    cx[i] = Cofx(i);
    cy[i] = Cofy(i);
    cz[i] = Cofz(i);
  }

  // 构造单段多项式轨迹
  PolynomialTraj poly_traj;
  poly_traj.addSegment(cx, cy, cz, t);

  return poly_traj;
}