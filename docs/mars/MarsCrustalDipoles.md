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

**当前 `9999setup/Dipoles` 已替换为用户提供的 Langlais (2004) M23/−20/14
实测反演模型的格式转换结果，包含全部 4840 个源，不再是合成示例。** 原始解和
用户放在其前面的 README 保存在 `9999setup/2003JE002048_solutionM23_20_14.txt`，
配套作者代码保存在 `9999setup/forward.c`，均未修改。

原始文件（包括前置 README）SHA256：
`6af5244ba33f43c07069cc1ab16139c5ebe1ed223b7c498087c9487a74d75443`。
配套 `forward.c` SHA256：
`79f5831d84946e00bad18b59502d2baec436bc159af604ccfac2039c528d1278`。
这些校验和记录本地用户提供的文件，不代表与出版社原始下载逐字节认证。

转换严格按 `forward.c` 的 `DIP_FORM=1` 路径：

- 数值头 `2.9189 4840 20.000`：平均间距（度）、源数、源深度（km）。
- 六列：东经（度）、纬度（度）、磁化层厚度（km）、Mr、Mtheta、Mphi（A/m）。
- `R_M=3393.5 km`，源半径 `R_M-depth=3373.5 km`，不使用求解器的参考长度替代。
- 作者代码的体积为 `V_km3=(R_M*pi/180*spacing_deg)^2*thickness_km`。
  当前全部厚度为 40 km，故每个源 `V=1195498.035071171 km³`（约 `1.195498e15 m³`）。
  换算采用正演代码这一确切公式，不采用论文摘要的近似体积描述，也不额外乘 π 或 3。
- 将球坐标磁化向量转为笛卡尔，再乘 `V_km3*1e9` 得到 A·m²；core=0。
  输出磁矩未吸收 μ0/(4π)。作者的 `scaling=100` 对应 km、A/m 输入转 nT 的系数，
  不能再乘到 SI 磁矩中。

复现转换（从 MPCNS_Mercury 目录）：

```sh
python3 Z5_Mars/tools/convert_langlais2004.py \
  Z5_Mars/9999setup/2003JE002048_solutionM23_20_14.txt \
  Z5_Mars/9999setup/Dipoles
```

转换器会跳过前置 README，并严格检查指定模型的头、4840 条记录、有限值与坐标范围。
它只转换格式，没有删除 1500 km 以外的源；cutoff 是求值时相对每个观测点的选择。
位置与磁矩保留火星体固坐标，未施加体固到太阳定向坐标的旋转；当前使用相当于将
体固轴与网格轴对齐。默认剩磁开关仍为 `is_Crustal_Bfield=0`，启用时改为 1。
原合成生成器 `make_demo_dipoles.py` 仅保留用于另行生成合成数据，不用于本模型转换。

本次验证直接编译用户提供的 `forward.c`，给其输入去掉 README 的数值副本，使用
`DIP_FORM=1`、作者原有的 1500 km 截断。分别在 200、400、1000 km 高度比较全球
10° 网格（每层 648 点，共 1944 点）的 Br/Btheta/Bphi，与转换后 SI 磁矩的标准
笛卡尔偶极正演比较，最大绝对差小于 `5.0e-6 nT`，即作者输出五位小数的舍入误差。
生产 `DipoleModel` 成功读取全部源，其 A 的数值 curl 在另选 9 点也符合该误差范围。
这验证了转换的单位、体积及方向；不表示跨 cutoff 球面的 CT 面通量与作者逐点截断
磁场完全相同，区别见下文。验证程序只在临时目录运行，未加入求解器构建。

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
当前转换文件将上述体固轴直接对齐模拟轴，属于人为指定的姿态。

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
源求和使用补偿累加；未配置 cutoff 时完整叠加全部源，当前示范 setup 使用 1500 km 截断。
初始化复杂度为 O(边数×源数)，含 ghost；无需在每个时间步重复计算。

可在 `3_MHD_Control.txt` 配置：

```text
double  Crustal_Dipoles_Cutoff_km  1500.0
```

`0` 表示不截断，参数缺省时取 0；当前示范 setup 按用户设置为 1500。正值为源位置到 Edge Gauss 积分点的三维直线距离上限，
单位 km；负值和非有限值报错。距离判定不包含数值软化 core，不是离地高度或球面弧长。
Langlais (2004) 第 17 段讨论的是 1500 km 半径的球形邻域；提供的 `forward.c`
明确使用 `thresh=1500`，按 `E2 <= thresh*thresh` 选择源，与本参数的距离定义相同。

此开关截断的是每个源的**磁矢势贡献**，之后统一同步 Edge 并由 CT 生成 Face 通量，
所以仍保持离散无散。它并不逐点复现作者直接对磁场贡献使用 cutoff 的结果：
硬截断 A 在截止球面不连续，会在相交网格面上引入额外磁场，且边积分可能对
网格分辨率敏感。若要求完整偶极外场，使用 `0`；不能将这个开关解释为不改变
物理场的纯加速选项。不要在 CT 后将超范围的 Face 磁场清零。
更改 cutoff 会使偶极缓存失效并重算；IMF 不受该开关影响。

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

## 编译与已有验证

当前仓库默认只构建 Mars 求解器，已移除独立测试程序与 CTest 注册：

```sh
cmake -S . -B build-mars -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_BUILD_TYPE=Release -DMPCNS_BUILD_MERCURY=OFF
cmake --build build-mars -j4
```

生成 `build-mars/Z5_Mars/MarsZ5`。显式关闭 Mercury 也适用于曾启用 Mercury 的旧
构建目录；新构建目录默认已关闭。`Z5_Mars/tools/make_demo_dipoles.py` 保留为
合成 Dipoles 输入的生成工具，不参与编译；当前实测模型转换使用
`Z5_Mars/tools/convert_langlais2004.py`，同样不参与编译。

移除测试程序前，已完成单偶极解析磁矢势、量纲归一化、源叠加与生产 CT 闭合检查；
当时的合成 4840 源倾斜网格的相对闭合残差为 `2.03e-17`。双 MPI 旋转接口检查覆盖背景切换、
缓存命中与失效、感应场续算保留、忽略保存的 Badd 及拒绝旧格式。
这些历史检查不覆盖球面极点或长期推进稳定性；本次真实模型转换与作者正演的比较见“文献与数据来源”。

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
cutoff 距离以及算法版本。缓存文件头包含标识、配置指纹、元素数量和数据校验和；读取时检查
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
