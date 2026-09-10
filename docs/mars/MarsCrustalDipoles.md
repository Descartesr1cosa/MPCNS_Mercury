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
bool    Read_Badd_Dipoles          1
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

缓存未命中或禁止读取缓存时，计算流程为：

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

**MarsZ5 采用与 Mercury 相同的感应场/背景场含义**：

```text
B_xi/eta/zeta    = Bind 的面通量（CT 推进的未知量）
Badd_xi/eta/zeta = IMF·S + CT(A_dipoles)（固定背景）
B_total         = Bind + Badd
Bind_cell       = 感应场的单元中心重构
B_cell          = 总场的单元中心重构
```

上面的 Edge → CT 流程只构造偶极部分；之后将当前 INITIAL 的 IMF 按面面积向量
投影加入 Badd 一次，再同步并检查完整背景的六面闭合。关闭剩磁开关时仍保留 IMF。
重算时 Aadd 保留偶极磁矢势环量；缓存命中时 Aadd 为零，不代表已重建磁矢势，
也不包含 IMF 的磁矢势。

电流由 `Bind_cell` 的离散 curl 计算。理想 UCT 中心物理通量、Hall 的 J×B、
Lorentz 力及快磁声时间步使用总场；理想 UCT 与 Hall 的磁场跳跃耗散使用感应场。
CT 只推进 Bind，背景不随时间演化。因此非均匀指定背景本身不产生求解电流，
但会通过总场电动势驱动感应场。Mars 的三面磁场重构和单元中心电流差分仍与
Mercury 的六面重构/DEC 电流算子不同，详见 [电磁离散对照](MarsMercuryEMDiscretization.md)。
`divB` 诊断为感应场散度，`divB_total` 为总场散度；Tecplot 同时输出 `B_cell`、
`Bind_cell` 及两种散度。Farfield 中 `Bind_cell` ghost=0、`B_cell` ghost=IMF，
边电场仍为固定 `-u_inf×IMF`，与 Mercury 的分工一致；外边界应远离强剩磁区。

无论从头运行还是续算，均根据**当前 setup** 重建完整 Badd：
从头运行令 Bind=0；新版续算保留检查点 Bind 及离子状态，仅替换背景。
相同配置不会重复加入 IMF 或偶极场，更换/关闭剩磁也不会清除已经演化的感应响应。
**bin 同时保存十种离子守恒量、`B_xi/eta/zeta` 感应面通量和完整
`Badd_xi/eta/zeta` 背景面通量，共 16 个字段。** 后处理可直接逐面相加得到
该检查点时刻的总磁通，无需再读取 setup 或偶极缓存。

**保存的 Badd 仅用于后处理，不参与续算背景的构造。** 续算读入字段后，
先将 Badd 清零，再由当前 INITIAL 的 IMF 与当前偶极模型/有效缓存重建。
不从感应场中加减保存的 Badd，也不读取任何旧 IMF 参数。
`Read_Badd_Dipoles` 只控制独立偶极缓存，和 bin 中保存的完整 Badd 无关。

这是一次不兼容旧火星结果的格式更新。Mars 检查点使用 `MPCNSRST` **version 2**，
在原有 nblock/nfield 后、字段记录前写入长度前缀字符串
`Mars_Bind_Badd_IMF_Dipoles_v1`，明确 B 为感应场、Badd 为 IMF 加偶极背景。
旧 version 1 总场文件、不兼容标记及缺少任一必要字段的文件直接拒绝，
不提供旧场转换或缺失背景的兼容分支。本次更新后应从头启动计算，后续可正常续算。
后处理读取工具需识别此额外头字段；Mercury 自身仍使用原 version 1。
续算应保持网格与归一化尺度不变。

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
分量交换。它先做初始化（`max_Nstep=0`），继而读取自身检查点，逐字节比较
所有 rank 的 conserved/face-flux 场文件，检测相同模型重复叠加及 ghost 不一致。
随后将模型磁矩加倍，并在续算中关闭、再开启剩磁，检查感应场和离子守恒量保持不变、关闭后 Badd 仅含 IMF。还检查更换 IMF 时
缓存命中且 IMF 只加入一次，以及单步推进生成非零感应场后更换背景仍保留感应场。
`--ranks 1` 可使用同进程的两个块接口。临时直角网格使用 `corner_treat=0`，
因此这些测试不覆盖球面极点或非规则多块奇异边，也不证明长期时间推进稳定性。

本次验证结果：Debug 完整构建成功；4/4 CTest 通过；单源及 4840 源含 ghost
倾斜网格闭合测试通过（多源相对残差 `2.03e-17`）；该测试的
AddressSanitizer/UndefinedBehaviorSanitizer 检查通过；双 MPI 进程的同向和
90° 旋转接口初始化/重启测试均通过，输出 Badd 闭合残差为 0，重启场文件一致。


## 2000 偶极子的计算量与缓存

两点 Gauss 每条边调用两次 Potential，每次累加所有偶极子，因此约有
`2 × N_dipoles × N_edges` 次单源求值。大三维块忽略边界/ghost 时
`N_edges ≈ 3 N_cells`，2000 个源约需 `12000 N_cells` 次求值。

本机对生产 `DipoleModel::EdgeIntegral` 用 `c++ -O3 -std=c++20` 做单线程微基准：
2000 源、20000 条边，耗时约 **0.288 s**（约 6.94 万边/s）。线性外推：
10 万单元约 **4.3 s**，100 万单元约 **43 s**，1000 万单元约 **432 s**。
这些是单核外推而非完整求解器计时，不包括 ghost、网格预处理、MPI 通信；并行
分摊源求值，但分块 ghost 和通信会增加开销，不保证按核数线性加速。

缓存三个方向的 double 面通量，忽略边界/ghost 约 `24 N_cells` 字节：
百万单元约 **24 MB**。本机 24 MB 缓存微基准写入（含校验）约 **0.029 s**，
读取（含内容校验）约 **0.035 s**；这是本地缓存文件系统结果，未计当前模型/网格
指纹计算和同步，不能当成集群冷读性能。结论：**反复启动的大网格值得缓存最终
偶极面通量；仅保存 2000 条偶极矩记录不能省掉主要求值成本。** 小网格只运行一次时收益有限。

`3_MHD_Control.txt` 新开关：

| 设置 | 行为 |
|---|---|
| `is_Crustal_Bfield=0` | 当前 Badd=IMF，不读取偶极模型或缓存 |
| `is_Crustal_Bfield=1, Read_Badd_Dipoles=0` | 强制由当前 Dipoles 积分并 CT，写入/更新缓存 |
| `is_Crustal_Bfield=1, Read_Badd_Dipoles=1` | 优先读取匹配缓存，失效则全部 rank 重算并更新 |

`Read_Badd_Dipoles` 缺省为 0；示范 setup 也保持 0。启用缓存读取可直接设为 1，
第一次没有缓存会自动生成，不需要先手动运行一次开关为 0 的版本。
文件按 MPI rank 写到 `DATA/Badd_Dipoles/dipole0000.bin`、`dipole0001.bin` 等。
这是已完成 CT/同步的 **仅偶极部分的面通量（含 ghost）**，不含 IMF，也不是源参数、总场或重启文件。
IMF 总是在读缓存或计算偶极场之后加入；单独改变 IMF（归一化尺度不变）不使缓存失效。

缓存校验覆盖：Dipoles 文件原始内容、L_ref/B_ref/μ0、MPI rank/进程数、活动块布局、
全部 Edge 端点坐标、当前 rank 的 inner/parallel/physical 边界文件及边界名称文件，
以及算法版本。缓存文件头包含标识、配置指纹、元素数量和数据校验和；读取时检查
长度、有限值、尾部、内容校验和。使用 FNV-1a 64 位检测意外损坏/过期，非安全用途。
为生成当前指纹，**即使命中缓存仍须提供当前 Dipoles 文本和网格**。

任一 rank 缓存不匹配、损坏或缺失，都统一重算，避免部分读缓存、部分重算的
通信顺序冲突。写缓存先写临时文件再 rename；写失败会提示，但继续使用有效的
已计算场。二进制格式使用本机字节序，跨平台或修改算法后应强制重算；修改积分或
CT 缓存约定时需递增代码中的 revision。缓存不随 checkpoint archive 自动复制，
缺失时可重建；检查点同时保存感应场与完整背景供后处理；续算背景仍由当前 setup 与有效缓存或重算结果恢复。

缓存命中后仍做偶极面通量的有向 Face 同步和单元闭合检查；加入 IMF 后再检查
完整 Badd。缓存不绕开无散检验。
日志输出 hit/miss 和各 rank 最大初始化耗时（含签名、读取/计算、同步及缓存写入，
不含后续 IMF 加入和完整背景检查）。

现在 Badd 的完整逻辑是：

1. fresh 令 Bind=0；restart 读取并保留 Bind；旧总场检查点拒绝读取。
2. 剩磁关闭则偶极部分为零；开启则按缓存开关读取或由当前模型重算。
3. 同步并验证偶极面通量；仅将偶极部分写入缓存。
4. 将当前 IMF·S 加入 Badd 一次，同步并验证完整背景无散。
5. CT 推进 Bind，物理项使用 Bind+Badd；检查点保存离子守恒量、Bind 与完整 Badd，后者仅供后处理。

新增 `MarsDipoleCacheTest` 验证缓存往返、错误签名/长度、截断、损坏、尾部数据和
非有限值拒绝；双 MPI 旋转接口 smoke 测试覆盖缓存命中、更改模型自动失效、
损坏或缺少单个 rank 缓存时统一重算，以及更换/关闭/重新开启背景、忽略保存的 Badd 和拒绝旧格式。

此次 IMF/感应场分裂回归：Debug 构建、4/4 CTest、双 MPI rank 的 90° 旋转
接口 smoke 测试通过。覆盖同配置续算不重复叠加、剩磁关闭后只保留 IMF、
只改 IMF 仍命中偶极缓存，以及一时间步生成感应场后更换
背景仍保留该感应场。静止均匀等离子体、非均匀偶极背景、Hall 开启时，
单步感应面通量保持在 `1e-14` 绝对容差内；这不是长期稳定性验证。

检查点回归还验证：bin 包含全部 16 个字段；将保存的 Badd 改为 `1e90` 后续算，
输出仍与当前 setup 重建结果逐字节一致；旧 version 1 和缺少 Badd 分量的文件
被拒绝。修改 IMF 后，输出的 Badd 通量差值等于 IMF 投影差值，感应场保持不变。
