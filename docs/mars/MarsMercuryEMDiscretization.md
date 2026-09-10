# Mars 与 Mercury 电磁离散对照

## 1. 本次统一的范围

`Z5_Mars` 的显式电磁更新现按 Mercury 的顺序组装：

1. 面上理想 MHD Rusanov 候选；
2. 从面候选构造 UCT 边电动势；
3. 电子压力边电动势；
4. 规则共享边的有向 alias 归约和奇异边组装；
5. Hall Rusanov 边电动势及其奇异边组装；
6. 固壁切向电动势条件；
7. 边场 owner/halo 同步，再用 CT curl 更新面磁通。

面工作量 `Eface_*` 始终在接收块的计算坐标基中本地构造，不跨块同步；真正的
边 1-form `E_*`、`Ehall_*` 才按拓扑方向同步。这一点对旋转块接口尤其重要。

## 2. 统一电子密度 limiter

对 Mars 的十种单价离子，先由各组分无量纲质量密度求数密度：

\[
n_s^*=\frac{\rho_s^*\rho_{ref}}{m_s n_{ref}},\qquad
n_e=\sum_s\max(n_s^*,0).
\]

与 Mercury 相同，Hall 和电子压力项统一使用平滑有效电子密度

\[
n_{e,eff}=\sqrt{n_e^2+n_{e,floor}^2},\qquad n_{e,floor}=10^{-2}.
\]

流体源项中分配给组分 $s$ 的电磁权重为

\[
w_s=\frac{n_s^*}{n_{e,eff}},\qquad
\sum_s w_s=\frac{n_e}{n_{e,eff}}.
\]

因此正常等离子体区 $w_s$ 退化为真实电荷比例；低密度区的 Hall、电子压力和
$\boldsymbol{J}\times\boldsymbol{B}$ 耦合会平滑衰减。混合离子速度仍使用真实比例
$\chi_s=n_s^*/n_e$，不允许用 density floor 改写组分比例；真空处速度置零。

## 3. 理想 UCT

每个坐标面先使用局部 Lax--Friedrichs/Rusanov 形式计算横向磁通。以法向
$\xi$、被输运分量 $\beta$ 为例：

\[
F^{\xi}_{\beta}=
\frac12\left[(u^\xi B^\beta-u^\beta B^\xi)_L+
(u^\xi B^\beta-u^\beta B^\xi)_R\right]
-\frac12 a_\xi(B^\beta_R-B^\beta_L),
\]

其中 $a_\xi=\max(|u_L^\xi|,|u_R^\xi|)$。随后把相邻四个面候选平均到边上，
并加入 Mercury 的 $C_{UCT}=0.5$ 横向迎风修正。例如

\[
\mathcal E_\xi=\overline{\mathcal E}_\xi+
C_{UCT}\left[
\frac12a_\eta\Delta_\eta B^\zeta-
\frac12a_\zeta\Delta_\zeta B^\eta\right].
\]

最终仍以离散 curl 更新面 2-form 磁通，保持 CT 的拓扑散度约束。

## 4. Hall Rusanov

单元中心 Hall 电场及 whistler 尺度为

\[
\boldsymbol E_H=\alpha\,\boldsymbol J\times\boldsymbol B,
\qquad
\alpha=\frac{C_H}{n_{e,eff}},
\qquad
\beta_H=|\alpha|\,|\boldsymbol B|.
\]

边中心项由周围四单元的 $\boldsymbol E_H$ 平均后沿边积分。以 $\xi$ 边为例，
Rusanov 修正为

\[
\mathcal E^H_\xi=
\langle\boldsymbol E_H\rangle_4\cdot\Delta\boldsymbol l_\xi
+\frac12\mu_\eta\Delta_\eta B^\zeta
-\frac12\mu_\zeta\Delta_\zeta B^\eta,
\]

\[
\mu_q=C_{wh}\frac{\beta_{H,e}}{h_q^2},\qquad
C_{wh}=0.1,
\]

其中 $\beta_{H,e}$ 取四单元最大值，$h_q$ 由相邻正交面面积除以边长得到。极点上
几何退化的切向 Hall 边按照拓扑定义显式置零；其他非正或非有限边长、面积和单元
Jacobian 会抛出包含 rank、block、索引和值的异常，不能以 `1e-30` 静默替代。

## 5. 时间步约束

按当前 Mars 设定，全局显式时间步仅为

\[
\Delta t=\min(\Delta t_{fluid},\Delta t_{fast\ MHD}).
\]

第一项逐离子计算 Euler 声学 CFL；第二项用总带电流体
$\rho_c=\sum_s\rho_s$、$p_c=\sum_s p_s$、混合离子速度和总磁场计算六个面的
快磁声谱半径。化学源项由局部点隐式更新，不进入全局 CFL。

**这里有意不施加 Hall/whistler 的 $h^2$ 显式时间步约束。** Hall Rusanov 的
$h^{-2}$ 数值耗散仍保留在空间离散中，但不参与 `ComputeTimestep()`。这符合当前
运行要求；若细网格 Hall 模式出现高频不稳定，应从解的诊断判断，而不是把极小
Hall 步长悄悄重新加入。

## 6. 仍然存在的 Mars--Mercury 差异

| 项目 | Mercury | Mars | 处理结论 |
|---|---|---|---|
| 带电组分 | H+、Na+ | 十种单价离子 | 用统一求和公式推广 |
| 外加/内部磁场 | `Bind + Badd` | 单一总磁场，可选剩磁初场 `Badd` | Mars fresh 初始化经 Edge A → CT 构造剩磁并加到总场；详见 [剩磁输入与来源](MarsCrustalDipoles.md) |
| Hall 径向 taper | 水星内部区域平滑关闭 | 无 | 不移植水星专用行星内部 taper |
| 电阻项 | 物理/人工电阻及内部电阻隐式求解 | 无 | 按 Mars 需求保持关闭 |
| Hall 时间步 | 当前 Mercury 时间步实现也未单列 Hall CFL | 仅 fluid/fast-MHD CFL | 不施加 $h^2$ 约束 |
| 单元磁场重构 | 六面最小二乘权重 | 三个正向面的局部 $3\times3$ 重构 | 尚未统一 |
| 电流离散 | face-to-edge DEC curl、Hodge 后重构 `J_cell` | 单元中心磁场的中心差分 curl | 尚未统一；Hall 中心项会受此差异影响 |
| 化学 | 两组分水星源项 | Fortran Mars 十离子局部点隐式网络 | 保留 Mars 物理 |

因此，本次已经统一的是理想 UCT、Hall Rusanov、电子密度 limiter、共享/奇异边
和固壁电动势路径；`B_cell/J_cell` 的重构算子仍是已明确记录的主要离散差异。

