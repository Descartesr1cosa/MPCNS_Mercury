# Mars 解析大气、光致电离与电离层内边界

## 1. 目的与适用范围

`Z5_Mars` 不再在运行时读取旧固定网格上的
`roeneu300km3X_EUVAC.txt`、`photoionization300km3X_EUVAC.txt` 和
`ionbound300km3X_EUVAC.txt`。现在每个 MPI 分区直接由火星中心坐标
计算中性大气与光致电离率，并在固体内边界求局地光化学平衡（PCE）。

这一实现面向全球 MHD 源项，而不是高光谱分辨率的电离层反演。其物理结构来自
MTGCM/SWIM 解析大气、EUVAC 光学深度方法和火星 MHD 常用的 100 km PCE
内边界；14 通道的顶层频率及少数组分的 100 km 锚点保留旧 EUVAC 预处理的
量级，以便与串行程序连续。文末列出了升级为完整 37 波段 EUVAC 所需的改动。

## 2. 中性原子和分子分布

令

\[
r=\|\boldsymbol{x}\|,\qquad z=(r-R_M)/1000,
\]

其中 $r,R_M$ 的单位为 m，$z$ 的单位为 km。默认
$R_M=3396\,\mathrm{km}$，模型下边界 $z_0=100\,\mathrm{km}$。

### 2.1 CO2、冷 O 和温度

采用 SWIM 比较工程中基于 Mars Thermospheric GCM 的解析拟合（密度单位
均为 $\mathrm{cm^{-3}}$）：

\[
n_{CO_2}(z)=6.04\times10^{18}e^{-z/6.98}
             +1.67\times10^{15}e^{-z/11.49},
\]

\[
n_{O,c}(z)=5.85\times10^{13}e^{-z/10.56}
           +7.02\times10^{9}e^{-z/33.97},
\]

\[
T_n(z)=196.95-64.56
\exp\!\left[-\frac12\left(\frac{z-115.7}{20.14}\right)^2\right]\quad\mathrm{K}.
\]

程序对温度施加 80 K 下限，并禁止把热层拟合外推到 100 km 以下。上述拟合及其
MTGCM 工况（$L_s=270^\circ$，火星处 $F_{10.7}=34$）可在
[SWIM/irfpy 模型说明](https://irfpy.irf.se/projects/planets/api/api_irfpy.mars.exosphere.html)
中核对。它是工程拟合而非新的经验反演，因此太阳活动、季节和地方时研究应换成
M-GITM/MTGCM 场或重新标定系数。

### 2.2 H 与热 O 日冕

H 使用冷、热两个球对称分量：

\[
n_H(r)=1.50\times10^5
 \exp\!\left[25965\left(\frac1{r_k}-\frac1{3593.5}\right)\right]
+1.90\times10^4
 \exp\!\left[10365\left(\frac1{r_k}-\frac1{3593.5}\right)\right],
\]

其中 $r_k=z+3393.5$（km）。热氧采用

\[
n_{O,h}(z)=5.23\times10^3e^{-z/626.2}
            +9.76\times10^2e^{-z/2790}
            +3.71\times10^4e^{-z/88.47}.
\]

总氧为 $n_O=n_{O,c}+n_{O,h}$。输出场的第 11、12 列分别保存冷 O 与热 O，
前十列仍严格按 `Ar CO2 CO C H N2 N NO O O2` 进入旧 Fortran 化学网络。
这些日冕表达式和系数来自
[Diéval et al. (2012)](https://link.springer.com/article/10.5047/eps.2011.08.015)；
观测表明 400 km 附近 O、H、H2 和 He 是质量密度的重要贡献者，且存在约
2–3 倍的不确定性，见
[Chaufray et al. (2009)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2008JE003130)。

### 2.3 其余组分的球对称扩散平衡

Ar、CO、C、N2、N、NO、O2 使用变温、球对称重力下的扩散平衡：

\[
n_s(r)T(r)=n_s(r_0)T(r_0)
\exp\!\left[-\frac{m_sGM_M}{k_B}
\int_{r_0}^{r}\frac{dr'}{T(r')r'^2}\right].
\]

程序用 32 段梯形积分计算公共积分。默认 100 km 锚点为：

| 组分 | Ar | CO | C | N2 | N | NO | O2 |
|---|---:|---:|---:|---:|---:|---:|---:|
| $n_s(z_0)\,[\mathrm{cm^{-3}}]$ | $1.80\times10^{10}$ | $1.93\times10^{10}$ | $3.49\times10^2$ | $2.32\times10^{10}$ | $4.50\times10^6$ | $1.04\times10^6$ | $1.03\times10^9$ |

这些锚点是与旧输入的兼容标定，不应解释为所有季节/太阳活动条件下的唯一值。
Fox 的详细火星光化学模型同样以多组分中性背景为输入，见
[Fox (2004)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2004JA010514)。

## 3. 空间光致电离率

### 3.1 当前 14 通道灰模型

对通道 $c$，母体组分为 $p(c)$，体产生率为

\[
q_c(\boldsymbol{x})=n_{p(c)}(\boldsymbol{x})\,J_{c,\infty}
\left(\frac{1.524}{d_{AU}}\right)^2e^{-\tau(\boldsymbol{x})}.
\]

14 个母体依次为

`Ar CO2 CO C H N2 N NO O O2 CO2 N2 CO2 CO2`，对应旧化学网络中的
Ar+、CO2+、CO+、C+、H+、N2+、N+、NO+、O+、O2+ 以及 CO2/N2 的四个
解离电离支路。$J_{c,\infty}$ 是旧 EUVAC 文件在光学薄区域的
$q_c/n_{p(c)}$ 回归值；它们作为代码常量保存，不再依赖旧网格文件。

光学深度直接沿指向太阳的球面射线积分，避免平板 $\sec\chi$ 在终止线发散：

\[
\tau(\boldsymbol{x})=
\int_0^{s_{top}}\sum_a\sigma_{a,eff}
n_a(\boldsymbol{x}+s\hat{\boldsymbol{s}})\,ds.
\]

当前灰模型对 Ar、CO2、CO、N2、O、O2 分别使用
$1.75,3.50,2.625,2.625,1.75,3.50\times10^{-18}\,\mathrm{cm^2}$ 的有效截面。
这些是把硬 XUV 与软 EUV 合并后的工程有效参数，不是单一波长实验截面。
积分步长从 1 km 几何增长到 100 km，顶界为 1500 km。若太阳方向射线穿过
半径 $R_M$ 的球体，则点位于几何阴影中，直接令 $q_c=0$。
PCE 求解时总产生率低于 $10^{-20}\,\mathrm{cm^{-3}\,s^{-1}}$ 的点按暗侧处理，
以避免对物理上为零、仅由浮点下溢留下的通量做病态非线性迭代。

### 3.2 与完整 EUVAC 的关系

EUVAC 原模型用 37 个波段及 $F_{10.7}$ 和其 81 日均值调节顶层光谱：

\[
F_{b,\infty}=F_{b,74113}\,[1+A_b(P-80)],\qquad
P=\frac{F_{10.7}+\overline{F}_{10.7,81}}2,
\]

完整波段形式应为

\[
q_c=n_{p(c)}\sum_b \sigma^{ion}_{p(c),c,b}F_{b,\infty}
\exp\!\left[-\sum_a\sigma^{abs}_{a,b}N_a\right].
\]

37 波段定义与截面来自
[Richards, Fennelly & Torr (1994)](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/94JA00518)。
当前灰模型保留完全相同的“母体密度 × 电离频率 × 光学衰减”结构；若后续加入
37 波段表，只需替换 `AtmosphereModel::Evaluate` 内的波段求和，无需改变并行场、
源项和边界接口。火星处 CO2/O 的典型光学薄电离频率及约 10% 次级电离贡献可参见
[Cravens et al. (2017)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1002/2017JA024582)。

## 4. 壁面离子分布

在每个壁面位置先把坐标沿径向投影到 $z_0$，随后使用同一点的中性密度和
光致电离率求十离子局地 PCE：

\[
P_s(\boldsymbol{n_i},\boldsymbol{n_n},\boldsymbol{q},T_i,T_e)
-L_s(\boldsymbol{n_i},\boldsymbol{n_n},T_i,T_e)=0,
\qquad n_e=\sum_s n_s.
\]

$P_s,L_s$ 就是 `MarsChemistry.cpp` 从 `time_advance2.for/source_species`
移植的全部反应网络。求解器在正数对数空间作阻尼固定点迭代，默认相对产生/损失
残差小于 $10^{-7}$；结果按“块号 + 壁面相邻单元索引”缓存。边界状态为

\[
\rho_s=m_sn_s,\quad \boldsymbol{u}_s=0,\quad
p_s=n_sk_BT_{wall},\quad
E_s=\frac{p_s}{\gamma-1}.
\]

H+ 延续原 Fortran 条件，单独取上游太阳风 H+ 密度的 0.3 倍。火星多组分 MHD
在约 100 km 内边界使用 PCE 是常见处理，可参见
[Ma et al. (2014)](https://doi.org/10.1002/2013ja019402)；PCE 主要适用于约
200 km 以下、输运尚未主导的区域，因此不应把此壁面表达式用作整个顶部电离层的
静态离子分布。

## 5. 参数与代码位置

运行参数在 `Z5_Mars/9999setup/4_Mars_Atmosphere.txt`：

- `mars_radius_m`：火星参考半径；
- `mars_lower_boundary_km`：PCE 下边界高度；
- `mars_orbit_au`：日心距离，光通量按 $d^{-2}$ 缩放；
- `mars_sun_x/y/z`：由火星指向太阳的单位方向（程序会自动归一化）。

实现位于 `Z5_Mars/include/MarsAtmosphereModel.h` 和
`Z5_Mars/src/MarsAtmosphereModel.cpp`；与 MPCNS 网格及壁面边界的连接位于
`Z5_Mars/src/MarsSolver.cpp`。

## 6. 已知限制

1. 中性大气当前球对称，不能表达季节、经纬度、地方时和沙尘暴变化。
2. 灰光学深度不能同时精确重现软 X 射线低层肩峰与各 EUV 通道峰高；定量电离层
   研究应升级为 37 波段 EUVAC/FISM2。
3. 几何阴影区只含数值密度底限，没有电子沉降、输运形成的夜侧电离层。
4. 次级电子电离尚未显式计算；若使用文献所述约 10% 修正，应避免与已标定的
   14 通道系数重复计算。
