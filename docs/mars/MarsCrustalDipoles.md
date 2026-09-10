# Mars 剩磁初场：Dipoles → Edge 磁矢势 → CT Face 通量

MarsZ5 新增 `Aadd_xi/eta/zeta` 与 `Badd_xi/eta/zeta`。独立输入文件
`Z5_Mars/9999setup/Dipoles` 支持任意数量的笛卡尔偶极子（最多 1000000），
无需为数千个源逐个定义 Parameters LIST。

## 文献与数据来源

1. **Langlais, B., Purucker, M. E., & Mandea, M. (2004)**,
   *Crustal magnetic field of Mars*, JGR 109, E02008,
   [doi:10.1029/2003JE002048](https://doi.org/10.1029/2003JE002048)。
   使用 MGS 观测反演等效源的位置网格上的磁化矢量；优选解
   M23/−20/14 包含 **4840 个**源，深度 **20 km**，平均间距约 **173 km**。
   这些是表示磁化地壳外场的等效源，不代表已确定的 4840 块独立磁性岩体。
   [出版社 Supporting Information](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2003JE002048)
   列出原始解 `2003JE002048_solutionM23_20_14.txt`（另有 `.14` 版本）、
   `2003JE002048_forward-code.c` 和 `2003JE002048-README.txt`。
2. **Langlais, B., Thébault, E., Houliez, A., Purucker, M. E., & Lillis, R. J. (2019)**,
   *A New Model of the Crustal Magnetic Field of Mars Using MGS and MAVEN*,
   JGR: Planets 124, 1542–1569,
   [doi:10.1029/2018JE005854](https://doi.org/10.1029/2018JE005854)。
   扩展等效源方法，结合 MGS 与 MAVEN；公开球谐解截断到 134 阶。
   [Zenodo 数据记录](https://zenodo.org/records/3876714)中的
   `Langlais2019.sh.gz` 是**球谐势系数**，不能直接当作本程序的偶极矩输入。

**本仓库的 Dipoles 是合成数值示例，不是上述任一论文的原始反演系数。**
2026-09-10 核查到了 2004 年附录文件目录，但实际下载返回 HTTP 403，
因此没有导入或验证真实 M23/−20/14 系数。不能用当前示例声称复现了火星观测剩磁。
默认 `is_Crustal_Bfield=0`。原始系数获得后，应结合作者 README 和正演代码核对
列定义、坐标基、体积因子、单位和有效半径，再转换到下述格式，并记录原文件校验和。
目前没有凭猜测实现原始文件列布局的转换器。

示例的每个位置和偶极矩由 `Z5_Mars/tools/make_demo_dipoles.py` 确定性生成：
4840 点 Fibonacci 球面；源半径 `3396 km - 20 km`；每点等效体积
`V=4πr_source²×40 km/N`；人工径向磁化
`M=2 exp(-((latitude+π/4)/0.45)²) sin(8 longitude) A/m`；`m=VM e_r`。
只有源数量和埋深参考 2004 年方法，布点方式和磁化强度/方向均为测试构造，
不等同于论文网格或数据。复现命令（从 MPCNS_Mercury 目录）：

```sh
python3 Z5_Mars/tools/make_demo_dipoles.py /tmp/Dipoles --count 4840
```

## 输入与坐标

在 `CASE/setup/3_MHD_Control.txt` 中启用：

```text
bool    is_Crustal_Bfield          1
string  Crustal_Dipoles_File       CASE/setup/Dipoles
```

路径相对于 MarsZ5 的 CASE 工作目录；也可写绝对路径。旧配置没有开关时维持关闭。
开关为 1 而路径省略时使用 `CASE/setup/Dipoles`。文件由每个 MPI rank 读取同一个
共享文件；运行期间不要修改。文件不加入 `filenames`，不经过 Parameters 解析器。
`PrepareCaseWorkdirIfNeeded` 原有的 setup 目录链接自动包含该文件；已有独立 CASE
需自行复制文件及开关。缺失、行数不符、非法值会报错并终止 MPI 作业。

```text
# 可使用空行及 # 行尾注释
MARS_DIPOLES_V1 2
# x_m       y_m z_m  mx_Am2 my_Am2 mz_Am2 core_m
0           0   0    0      0      1e15   0
1000000     0   0    0      1e14   0      0
```

每条记录严格七列。位置和 core 单位为 m；偶极矩单位 A·m²，**未吸收 μ0/(4π)**。
输入必须已经在模拟网格使用的、以火星中心为原点的右手笛卡尔坐标系中。
地理体固坐标的惯例可定义 X 指向纬度 0°/东经 0°，Y 指向东经 90°，Z 指北极。
若模拟采用太阳定向坐标，应按选定时刻的姿态对**位置及磁矩同时施加同一个旋转**。
程序不猜测季节、子太阳经纬度或 SPICE 姿态；本版本为固定姿态初场，不包含自转更新。
示例将上述体固轴直接对齐模拟轴，属于人为姿态。

若外部模型给出球坐标磁化 `M_r,M_theta,M_phi`，theta 为余纬、phi 为东经，则

```text
e_r     = (sinθ cosφ, sinθ sinφ, cosθ)
e_theta = (cosθ cosφ, cosθ sinφ, -sinθ)  # 向南
e_phi   = (-sinφ, cosφ, 0)               # 向东
position = r_source e_r
moment   = V (M_r e_r + M_theta e_theta + M_phi e_phi)
```

V 必须采用源模型自己的体积/权重，km³ 转 m³ 需乘 10⁹；不能用示例的等面积体积
替代作者源权重。若作者文件已经给偶极矩，则不能再乘体积。北向分量与 theta
分量符号相反。原始模型参考半径不能直接替换成 MarsZ5 的 `L_ref`。

## 离散方法和场含义

每个偶极子的磁矢势（SI）为

```text
A(x) = μ0/(4π) Σ m_s × (x-r_s) / (|x-r_s|² + a_s²)^(3/2).
```

输入时预计算 `r_s/L_ref`、`a_s/L_ref`、`μ0 m_s/(4π B_ref L_ref³)`，
得到无量纲 `A/(B_ref L_ref)`。物理源通常取 `a_s=0` 且位于计算域（包括几何
ghost）之外；正 core 是可选的数值软化，会改变源附近的物理模型，不是文献默认参数。
点源不能落在积分边上。非常接近源的边需加密或另行提高积分精度。

每条直边使用对称两点 Gauss 积分近似 `∫ A·dl`，返回无量纲边 1-form
（物理尺度 `B_ref L_ref²`）；退化成一个点的边积分为零。
源求和使用补偿累加，完整叠加全部源，不使用空间硬截断。
初始化复杂度为 O(边数×源数)，含 ghost；无需在每个时间步重复计算。

流程为：

1. 计算全部 Edge（含 ghost）的 Aadd 环量。
2. 用现有 orientation-aware EdgeOwner 和 Halo 同步，共享边采用同一有向值；
   不对 Aadd 施加 E 的固壁/远场边界条件。
3. 调用 **`CTOperators::CurlEdgeToFace(...,+1,true)`**，得到面 2-form
   `Badd=∮ A·dl`。公共 CT 新增可选 `include_ghost` 参数，原有调用默认仍只写内区。
4. 同步 Badd 面通量。没有 Badd 物理边界覆盖操作。
5. 检查所有活动单元六面有向通量之和；日志输出最大绝对残差以及相对全域最大
   六面绝对通量和的残差。相对值超过 `1e-11` 时终止，捕捉接口不一致。

`div(curl)=0` 来自每条单元边在六个面环量中正负抵消，因此无散性不依赖
Gauss 积分是否精确；它不代表场幅值已收敛。真正的离散散度为通量和除以单元体积。
浮点抵消留下舍入误差。MPI/旋转接口必须先统一 A 的边方向和别名值。

**MarsZ5 保留总场推进约定**：`B_xi/eta/zeta` 存总场，首次启动时
`B_total = IMF·S + Badd`；Badd 独立固定保存，Aadd 保留用于检查。
后续 CT 更新总场，`B_cell`、电流、理想 UCT、Hall、Lorentz 力和 divB 都沿用总场，
不得在这些算子里再次加 Badd。本次没有引入 Mercury 的 Bind/Badd 背景场分裂，
也没有从离散电流中扣除剩磁的离散 curl。因此解析无电流的剩磁仍可能产生
有限网格误差电流；若需要消除此类背景电流，需另行统一修改感应场分裂离散。
原有固定 IMF 远场 ghost/边电场约定保持，外边界应放到剩磁已充分衰减的位置。

重启文件继续只保存已有的总 B 面通量。续算时可以重建 Aadd/Badd 供检查，
但**绝不重新叠加到已读入总场**。关闭开关不会删除检查点中原有剩磁；更改 Dipoles
也不会替换检查点初场。应使用同一 Dipoles 和姿态续算；变更初场需 fresh run。
`Badd_*`/`Aadd_*` 在 Field 中可访问，现有默认 Tecplot 仍输出总 `B_cell`。

## 验证

```sh
cmake -S . -B build-mars -DMPCNS_BUILD_MERCURY=OFF -DCMAKE_BUILD_TYPE=Debug
cmake --build build-mars -j4
ctest --test-dir build-mars --output-on-failure
```

`MarsDipolesTest` 使用生产 CT 算子，检查单偶极子解析磁矢势/解析 curl、
量纲归一化、源线性叠加、边反向与退化边、输入错误与奇点、4840 源输入，
并在倾斜网格上验证单源及多源 `div(curl A)` 的舍入级闭合。
真实文献模型的幅值/地图验证需要取得原始数据及作者正演参考；合成测试不能替代它。

完整求解器的小网格 fresh/restart 检查（需要可运行的 MPI 环境）：

```sh
python3 Z5_Mars/tests/mars_dipoles_smoke.py build-mars/Z5_Mars/MarsZ5 --ranks 2
python3 Z5_Mars/tests/mars_dipoles_smoke.py build-mars/Z5_Mars/MarsZ5 --ranks 2 --rotated
```

脚本只在临时目录生成两个 5×5×5 单元块，使用全部 4840 源。非周期接口的
通信标记为偶数；`--rotated` 将第二块坐标轴旋转 90°，覆盖 Edge 方向与 Face
分量交换。它只做初始化（`max_Nstep=0`），继而读取自身检查点，逐字节比较
所有 rank 的 conserved/face-flux 场文件，检测重复叠加及 ghost 不一致。
`--ranks 1` 可使用同进程的两个块接口。临时直角网格使用 `corner_treat=0`，
因此这些测试不覆盖球面极点或非规则多块奇异边，也不证明长期时间推进稳定性。

本次验证结果：Debug 完整构建成功；3/3 CTest 通过；单源及 4840 源含 ghost
倾斜网格闭合测试通过（多源相对残差 `2.03e-17`）；该测试的
AddressSanitizer/UndefinedBehaviorSanitizer 检查通过；双 MPI 进程的同向和
90° 旋转接口初始化/重启测试均通过，输出 Badd 闭合残差为 0，重启场文件一致。
