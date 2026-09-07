# Mars 续算边界与 Backup/Archive 诊断

## 技术结论

零初始化数据中的来流 IMF 与配置一致，未发现 28 核文件缺失、步数错位或初始 `By` 投影错误。续算后偏离的最强解释是旧 Mars 将 Farfield 流体、磁场和电场全部当作零梯度复制边界：内部扰动会进入来流虚网格，重启又会完整恢复这些历史虚网格值。

修复采用 Mercury 的边界结构：固定太阳风流体虚网格、固定 `B_cell` 来流、令边界边电场等于 `-u_inf × B_inf` 的线积分。为兼容旧检查点，Farfield 面磁通只在读取检查点后的第一次 `Bface` 同步中恢复为配置 IMF，之后仍由 CT 演化，不在每一步强制覆盖。

Archive 参数此前已经被读取，但 Mars 时间推进器从未调用归档判定和写出函数，所以功能实际上没有运行。现已补齐与 Mercury 同构的 Archive 与 `DATA_backup` 输出路径。

## 零初始化的 IMF 参数与网格投影正确

`build/Z5_Mars/Para/DATA` 中全部 28 个 rank 的重启文件头一致：magic 为 `MPCNSRST`、版本 1、step 500、物理时间 `0.0027474036964523774`、重启场数 13。rank 0 和 14 各有两个 block，其余 rank 各有一个，没有发现缺 rank 或检查点代次混杂。

INITIAL 给出的 IMF 除以参考磁场后为 `(-0.56, 0.83, 0)`。按 CASE 边界映射筛选 Farfield 外边界，共得到 972 个 Tecplot 节点：

| 分量 | 目标值 | 最小值 | 最大值 | 均值 | 平均绝对误差 | 最大绝对误差 |
|---|---:|---:|---:|---:|---:|---:|
| Bx | -0.56 | -0.5645875931 | -0.5563586354 | -0.5598954321 | 0.0012545422 | 0.0045875931 |
| By | 0.83 | 0.8234012127 | 0.8333641887 | 0.8271451865 | 0.0029488296 | 0.0065987873 |
| Bz | 0 | -0.0066437549 | 0.0058529465 | -0.00001453755 | 0.0005672344 | 0.0066437549 |

`By` 的均值偏差约为 `-0.00285`，最大绝对偏差约为 `0.00660`，小于目标幅值的 0.8%。这支持“从零初始化正确”的测试观察，也排除了 INITIAL 读值和基本曲面投影是主要原因。

## 偏差来源位于运行期 Farfield 策略

原 Mars 为所有物理边界场预注册 `DefaultPhysicalCopy`，却没有像 Mercury 一样覆盖 Farfield 的 `U_H`、`B_cell` 和 `Eedge`。对太阳风入口而言，这相当于把最靠近边界的内部状态复制到虚网格，而不是持续给定来流。

影响链条为：

1. 内部流体和磁场扰动到达外边界。
2. copy 边界把扰动复制进来流虚网格。
3. 边界电场没有保持来流 `-u×B` 的零旋度构造，CT 可以逐步改变边界磁通。
4. 重启文件保存场的完整分配区，续算会恢复已有的虚网格和面磁通偏差。
5. 续算初始化立即同步 `Bface` 并重构 `B_cell`，因此旧边界历史会继续参与计算。

新实现采用以下处理：

- H⁺ Farfield 使用固定太阳风密度、速度和压力；其余离子保持 Mars 零初始化时的低密度静止状态。
- `B_cell` Farfield 虚网格固定为配置 IMF。
- `E_xi/E_eta/E_zeta` Farfield 使用 Mercury 的边线积分，使边界电场对应 `-u_inf × B_inf`。
- 旧检查点可能已经含有漂移的法向面磁通，因此续算后的第一次 `Bface` 同步会按 `B_imf·S` 校正 Farfield 面磁通；完成后关闭该修复，后续保持 Mercury 的 Bface no-op/CT 演化策略。

## Archive 已配置但此前未被执行

当前 CASE 设置为：

- `Archive_Output_Time = 0.5`：每 0.5 小时墙钟时间归档。
- `Archive_Output_PhyTime = 2.0`：物理时间间隔参数。
- `Archive_Output_By_PhyTime = 0`：当前不按物理时间触发，而按墙钟触发。

公共 I/O 模块原本能够读取并打印这些设置，但旧 Mars 的 `Advance()` 只在普通输出点写 `DATA`，从未调用 `IsArchiveDue()`、`ArchiveDataDirectory()` 或 `BackupDataDirectory()`。因此当前 Para 目录只有 `DATA`、没有 `DATA_archive` 和 `DATA_backup`，原因是 Mars 驱动代码路径缺失，不是设置文件解析失败。

修复后：

- 每次普通检查点写完后，将完整检查点镜像到 `DATA_backup`。
- Archive 到期时直接写独立历史目录，不依赖 `DATA` 或 `DATA_backup`。
- Archive 与普通输出在同一步到期时，复用已序列化的 Archive 检查点，避免重复写二进制场文件。
- 续算读入的物理时间被设为下一轮物理时间归档的起点，与 Mercury 一致。

## 数据范围与方法

分析范围是 `build/Z5_Mars/Para/DATA` 中 28 组 rank 文件以及对应 CASE。Farfield 统计只包含边界映射标记为 `Farfield` 的外边界节点；Outflow 节点不计入该表。

重启文件按公共 I/O 的二进制写入格式核对头字段。Tecplot 文件按仓库 TDV112 writer 的变量和 zone 顺序解析，再结合 28 核物理边界映射筛选节点。代码侧逐项比较了 Mars 与 Mercury 的驱动输出通道、初始化/续算顺序和 Farfield handler 注册。

## 局限与验证建议

现有 DATA 是零初始化后的 step 500 快照，没有包含发生 `By` 明显偏离时的续算 Tecplot 文件。因此当前证据能证明零初始化检查点一致，并能证明旧边界策略允许污染，但不能直接量化偏差第一次出现在哪一步或哪个边界单元。单一时间截面也不适合绘制趋势图。

更新后的集群验证应检查：

1. 用旧检查点续算，在首次派生量更新后检查 Farfield `By`，应恢复到约 `0.83`。
2. 前进一步后再次检查 Farfield `By` 和 `|divB|max`，确认一次性面磁通修复没有引入不可接受的散度。
3. 普通输出后确认生成 `DATA_backup`。
4. 累计 0.5 小时墙钟后确认生成新的 `DATA_archive` 子目录。
5. 若 `By` 仍快速偏离，保存偏离前后相邻两个输出时刻，进一步区分 CT 边界演化与块接口 oriented face-owner 变换。

## 证据来源

- `build/Z5_Mars/Para/DATA`
- `build/Z5_Mars/Para/CASE/setup/0_General_Setting.txt`
- `build/Z5_Mars/Para/CASE/setup/1_Initial_condition.txt`
- `build/Z5_Mars/Para/CASE/geometry/boundary_condition`
- `Z5_Mars/src/MarsSolver.cpp`
- `Z5_Mars/src/MarsSolverAdvance.cpp`
- `Z4_Mercury/src/MercurySolver_Driver.cpp`
- `Z4_Mercury/src/boundary/Boundary_Condition.cpp`
- `Z4_Mercury/src/boundary/Boundary_Condition_Fluid.cpp`
