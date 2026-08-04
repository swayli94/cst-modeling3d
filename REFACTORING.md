# 重构记录（refactor 分支）

本文档记录本轮对 cst-modeling3d 仓库的整体重构：每一项改动、原因、对兼容性和数值结果的影响，以及验证方法。供逐条检查。

---

## 1. 目标与验证方法

目标：解决早期编程习惯留下的结构问题，删除旧的 `surface.py`（`surface2.py` 转正），修复已确认的 bug，保证 `example/` 中示例可复现，同步更新文档。

验证方法（每个阶段执行一次）：

1. 重构前先在原始代码上运行全部 32 个示例脚本，保存全部输出文件作为**基线**；
2. 每个阶段结束后重跑全部示例，用数值容差（abs 1e-9 / rel 1e-7）逐文件对比基线；
3. 阶段 1–3 的要求是**零数值差异**，阶段 4 只新增之前无法运行的示例的输出。

验证脚本已收入仓库，可复现整个对比过程：

```bash
./scripts/run_all_examples.sh <示例树目录>          # 运行全部示例，写 _results.txt
python scripts/compare_outputs.py <基线目录> <新目录> # 数值对比输出
```

最终结果：**基线中的 111 个输出文件全部零差异**；通过示例从 24 个增加到 **28 个**；剩余 4 个失败示例均因仓库之外的外部依赖（见第 6 节）。

注意：示例依赖仓库内的 `cst_modeling`。系统里 pip 安装的是 0.2.4（仓库代码已升至 1.0.0），运行脚本已用 `PYTHONPATH` 指向仓库；手动跑示例时请 `pip install -e .` 或设置 `PYTHONPATH` 为仓库根目录。

---

## 2. 模块结构变化

| 变化 | 说明 |
|---|---|
| `surface2.py` → `surface.py` | 新的放样（Lofting）架构成为正式的 `cst_modeling.surface` 模块 |
| 旧 `surface.py` 删除 | 其中的经典 `Surface`/`OpenSurface` 类（线性截面插值架构）**原样移入 `basic.py`**；`surf_axisymmetric` 移入 `math.py` |
| `BasicSection`、`interpolate_basic_sections`：`basic.py` → `section.py` | 所有 section 类集中在一个模块；`basic.py` 保留再导出，`from cst_modeling.basic import BasicSection` 仍然可用 |
| `__init__.py` 重写 | 删除 `__name__ = "cst_modeling"` 这种无意义赋值；坐标系约定文本移入模块 docstring；新增带 `__all__` 的公共 API 导出 |

顶层公共 API（`from cst_modeling import ...`）导出的 `Surface`/`OpenSurface`/`BasicSurface` 是**新架构**的类；经典架构的同名类通过 `cst_modeling.basic` 导入。

重构后的分层（无循环依赖）：

```text
math ← section ← basic（经典曲面）;  math ← operation;  io 独立
surface（新架构）← section, operation, io;  foil ← section, math;  tools/* 在最上层
```

此前 `section.py` 依赖 `basic.py`（为了 `BasicSection`），而 `basic.py` 又包含曲面类，层次是倒置的；`BasicSection` 移入 `section.py` 后层次理顺。

### 导入速查（旧写法 → 新写法）

| 旧 | 新 |
|---|---|
| `from cst_modeling.surface2 import Surface` | `from cst_modeling.surface import Surface` |
| `from cst_modeling.surface import Surface`（经典架构） | `from cst_modeling.basic import Surface` |
| `from cst_modeling.surface import surf_axisymmetric` | `from cst_modeling.math import surf_axisymmetric` |
| `from cst_modeling.foil import unify_foil` | `from cst_modeling.section import normalize_foil`（该函数早已改名，本次只是把引用改对） |
| `from cst_modeling.basic import BasicSection` | 仍可用（兼容再导出），新代码建议 `from cst_modeling.section import BasicSection` |

### ⚠️ 需要明确知晓的不兼容点

**`from cst_modeling.surface import Surface` 的语义变了**：以前得到旧架构的类，现在得到新架构的类。仓库内所有引用（示例、docs、tools/auxiliary.py）已全部迁移；如果仓库之外的代码用了旧写法，需要按上表改导入。这正是"删除旧 surface.py"的必然结果，请确认可以接受。

其余兼容处理：`basic.py` 继续再导出 `rotate`/`transform_foil`/`stretch_fixed_point`/`to_cylinder`/`BasicSection`/`interpolate_basic_sections`（blwf 示例和 docs 脚本在用），`foil.py` 继续再导出 `cst_foil`/`cst_foil_fit`。

---

## 3. Bug 修复（逐条）

全部修复后 111 个基线输出文件**零差异**——包括理论上可能改变行为的两处（第 1、12 条），现有示例均未触发被修复的错误分支。

1. **`math.transform_curve`：`rot_x`/`rot_y` 分支误用 `rot_z` 角度旋转转轴**（复制粘贴错误）。修复后仅当 `rot_axis` 与 `rot_x`/`rot_y` 同时非零时行为才变化——现有示例没有这种用法（已逐一核对），但这是**真实的行为修正**，如果你在外部代码中同时用过这几个参数，结果会变（变为几何上正确的结果）。
2. **`section.foil_increment`**：下表面分支的守卫误写为 `cst_u is not None`，导致只传 `cst_l` 时被静默忽略、只传 `cst_u` 时崩溃。已修正为 `cst_l is not None`。
3. **`operation.GuideCurve.update_section_with_value`**：`sections=None` 默认路径必然 TypeError（元组乘整数用作索引）；修复全范围路径，并处理了全范围时 `smooth_omega_shape_function` 产生 0/0 NaN 的潜在问题。
4. **`basic.BasicSurface.bend`**：`point_ = point.append(chord)`（`append` 返回 `None`），该分支必然崩溃。已修复。
5. **`basic.BasicSurface.linear_interpolate_z`**：`z` 恰好等于最后一个控制截面的 `zLE` 时必然抛 IndexError。已修复，内部情形行为不变。
6. **`io.read_tecplot`**：用 `line[0] in 'Variables='` 这类子串包含判断识别段落头（单字符即可误判）。改为 `startswith` 关键字比较；对本包写出的文件解析结果不变。
7. **分类字符串的 `in` 反模式**全部改为精确比较：`math.rotate` 的 `axis in 'X'`、`basic.add_section` 的 `axis in 'Y'`、两处 `plot` 的 `type in 'wireframe'`；顺带发现 **`tools/xfoil.py` 的平台判断 `platform.system() in 'Windows'`** ——空字符串也会判真，是真实的误分类，已改为 `==`。
8. **`basic.BasicSurface.smooth`**：`_yz` 靠上一次循环迭代的变量泄漏才能在最后一个点工作。改为显式计算，已验证输出逐位一致。
9. **`foil.check_validity`** 会**原地修改调用者的 `x` 数组**（把首末点写成 0/1）。改为不修改；两个调用类（`FoilGeoFeatures`/`FoilModification`）本来依赖该副作用，已在调用处显式归一化。⚠️ 如果外部代码依赖这个副作用，需要自己归一化。
10. **新 `surface.py`**：过时桩方法 `add_section` 的"print 再 raise"改为直接 `raise NotImplementedError`；`split_surface` 的报错信息把弦向误写为展向，已改正文案。
11. **可变默认参数**（`origin=[0,0,0]`、`section_s_loc=[0.0,1.0]`、`zone_id=[]` 等，遍布 math/operation/surface/tools）：全部改为 `None` 哨兵。当前无人变异它们，行为不变。
12. **`section.RoundTipSection.__init__`**：(a) 厚度重缩放分支把实际单位的 `dy_` 混进单位弦长相对坐标，已统一为单位弦长约定；(b) 链式比较 `np.max(np.abs(cst_l))>1E-6 > 0` 化简（原本恰好等价，纯清理）。⚠️ (a) 是真实行为修正，但该分支在仓库示例中从不触发（示例只用 `RoundTipSection` 的静态方法，不实例化）；如果你在外部用过"指定厚度重缩放"构造路径，结果会变（变为单位正确的结果）。

---

## 4. 清理与去重

- **删除死代码**：`math.rotation_3d`（全仓库无调用）、`math.clustcos`（已被矢量化的 `dist_clustcos` 取代，`section.py` 对它的导入一并删除）、`math.py` 末尾的 `__main__` 绘图演示块、多处注释掉的调试代码块（`foil.py`/`math.py`/`surface.py`）。
- **`scaled_sigmoid` 去重**：`basic.py` 中的静态方法副本删除，统一用 `math.scaled_sigmoid`（公式逐字相同）。
- **`GuideCurve.generate_by_spline` / `update_by_spline` 去重**：提取公共的样条构建核心 `_build_spline`，两个公开方法行为不变。
- **6 份 `read_setting` 解析器合一**：`basic.py` 的 `BasicSurface`/经典 `Surface`/`经典 OpenSurface` 与新 `surface.py` 的两个类各自维护一份 ~120 行的手工逐行解析，现在共用 `io.py` 的解析器 `read_surface_setting_blocks` + 布局应用辅助函数 `apply_surface_layout`（后续命名清理中从 `basic.py` 移至 `io.py` 并转为公开函数）；每个类保留自己的语义差异（例如经典 `Surface` 把 `specified_thickness <= 0` 视为 `None`，新架构只在有第 6 列时按 `tt >= 0` 判断——**有意保留了这个差异**，未强行统一）。已用新旧包做解析结果差分测试（24 个用例全部一致）。
- **文档字符串**：补上 `transform_curve` 的 `dz`、`io.output_curve` 的 `fname`、`GuideCurve.__call__` 说明；修正 `output_plot3d_for_parts` 文档中引用不存在函数的问题。

有意未做（避免过度重构）：`flip`/`translate`/`split`/`plot` 在经典（basic.py）与新（surface.py）架构间仍有近似重复，两架构的实现路径不同，合并收益小风险大；`CST_flip:` 设置块的解析结果历来被丢弃（未完成的功能），保留了解析跳过逻辑以免破坏含该块的设置文件。

---

## 5. 示例修复（重构前就跑不通的）

| 示例 | 问题 | 修复 |
|---|---|---|
| `blade/test-blade-simple.py` | 调用 `smooth(isec0=, isec1=)`，v2.0 时参数已改名 `i_sec0/i_sec1`（语义相同），示例没跟上 | 改示例 kwargs |
| `fairing/test-fairing.py` | Windows 路径 `'..\fuselage\tube-curve.dat'` | 改为 `'../fuselage/tube-curve.dat'`（该数据文件本就提交在 fuselage/ 下） |
| `fairing/test-fairing-detail.py` | 同上，读当前目录的 `tube-curve.dat` | 同上 |
| `blwf/surface2blwf.py` | 导入早已改名的 `unify_foil`；还有一个没用的 `from matplotlib.pyplot import axis` | 改用 `section.normalize_foil`，删除多余 import |
| `format-transfer/grd2igs.py` | 期待当前目录有个 `wing.xyz`，仓库里没有 | 改为自包含：先用 `surf_axisymmetric` 生成演示 plot3d 文件再转 IGES |

这 4 个示例（blwf 除外）现在全部跑通；它们的输出是"新增"的，没有旧基线可对比——因为它们在重构前就无法运行。

---

## 6. 仍然无法复现的示例（外部依赖，请你判断）

| 示例 | 缺失 |
|---|---|
| `blwf/run-blwf58.py` | `input.txt`（应由外层优化驱动提供）、`blwf58.exe`（Windows 可执行文件，仓库中没有）、私有包 `cfdpost` |
| `blwf/surface2blwf.py` | `surface.dat`（CFD 解文件）、`blwf58.exe` |
| `xfoil/test-xfoil-aoas.py`、`test-xfoil-cls.py` | 需要 xfoil 可执行文件。仓库只带了 Windows 版 `xfoil.exe`；Linux 下 `tools/xfoil.run_xfoil` 会在 PATH 上找 `xfoil.x`/`xfoil`。装上 xfoil 后这两个示例应当能跑（未验证） |

这 4 个脚本本身的代码已修到"只缺外部输入"的状态。

另外两件值得知道的事：

- `wing-basic/simple-wing-old.py`（经典 API）与 `simple-wing.py`（新 API）输出**同名文件**（`Wing-basic.dat` 等），同一目录里先后运行会相互覆盖。这是重构前就有的情况，未改动；如果不想要这个遗留示例，删掉它即可（docs 的 wing 教程页有提及，需同步删一句）。
- 仓库里提交的部分 PNG（如 airfoil-features、basis-functions）与当前环境重新生成的图有字节差异——这是历史提交图与当前 matplotlib 版本的渲染差异，与本次重构无关（重构前后的重新生成结果完全一致）。

---

## 7. 文档更新

- `docs/cst_modeling.rst`：删除 `surface2` 的 automodule；补上此前漏掉的 `tools.xfoil`、`tools.blwf` automodule。
- `docs/source/index.rst`：项目结构树与模块简介同步更新。
- `docs/source/tutorial/surface.rst`：该页讲的本来就是经典架构，现在显式标注类来自 `cst_modeling.basic`，并在页首加了指向新架构 `cst_modeling.surface` 的说明。
- `wing.rst`/`propeller.rst`/`airfoil.rst` 中 `cst_modeling.surface2.*` 交叉引用全部改为 `cst_modeling.surface.*`。
- `docs/tutorial/others.rst`：删掉了介绍不存在函数 `read_curves` 的小节（`io` 中没有等价物）。
- `docs/tutorial/*.py` 与 `docs/example/*.py` 脚本：修复全部失效导入——既有本次迁移造成的（surface→basic/math），也有**重构前就坏的**（`curve_curvature`、`intersect_*`、`from_cylinder`、`output_curve` 等从没被 `basic.py` 再导出，全部改从 `math`/`io` 导入）。
- `README.md`：两个代码片段的导入更新（新架构 `cst_modeling.surface`、经典架构 `cst_modeling.basic`），顺带修了叶片片段里同样过时的 `isec0/isec1` kwargs。
- 文档构建验证：全新 Sphinx 构建通过；仅剩 10 条 `n_point` 交叉引用歧义警告，**重构前就存在**（纯净代码同样 10 条），与本次改动无关。起因是 `Lofting2Profile` 和 `LoftingRevolution` 都文档化了同名属性。RTD 配置了 `fail_on_warning: true`，如果线上构建因此失败，这是个预存在问题，可以后续单独处理。

---

## 8. 改动规模

```text
cst_modeling/  13 个文件，+2061/-2419 行（净减约 360 行，另含 surface2→surface 的整体迁移）
example/       17 个示例脚本（多数仅改导入，其中 5 个另有实质修复，见第 5 节）
docs/ + README 17 个文件
scripts/       新增示例回归验证脚本（run_all_examples.sh, compare_outputs.py）
```

---

## 9. 命名规范统一（第二轮，commit 816f56e 之后）

对全部公开 API 做了一轮命名审查与改名，57 个文件、约 ±1700 行，**纯改名无逻辑变化**，回归验证 111 个输出文件零差异，Sphinx 构建无新增警告。改名原则：统一两套架构的词汇、消除 camelCase、消除同名异义；保留成熟的领域约定。

### 9.1 两套架构词汇统一（经典 basic.py → 与新架构一致）

| 旧名 | 新名 |
|---|---|
| `n_sec`（构造参数/属性） | `n_section` |
| `.secs` | `.sections` |
| `.surfs` | `.surfaces`（含 operation.py 里新架构自己用旧名的 `Lofting.surfs` 等） |
| `l2d` | `is_2d` |
| `layout_center()` | `update_layout_center()` |
| 新架构的 `update_section()` | `update_sections()`（与经典架构统一为复数） |
| `zLEs` 属性 | `spanwise_locations` |
| `half_span` | `half_size` |
| `add_sec` | `add_section` |
| `section2surf` / `section_surf_axisymmetric` | `sections_to_surface` / `sections_to_axisymmetric_surface` |
| `interp_basic_sec` | `interpolate_basic_sections` |

### 9.2 camelCase → snake_case

| 旧名 | 新名 |
|---|---|
| `fromCylinder` / `toCylinder` | `from_cylinder` / `to_cylinder` |
| `lTwistAroundLE` | `twist_around_LE`（保留 LE 后缀约定） |
| `TwoD`（output_section 参数） | `is_2d` |
| `RTE` / `RTE_`（smooth_axisymmetric 参数） | `r_TE` / `r_TE_inner` |
| `dX, dY, dZ`（translate 参数） | `dx, dy, dz` |
| `interpolate_IDW` | `interpolate_idw` |
| `iLine0` | `i_line0` |
| `plane_3points(P0,P1,P3,...)` | `(p0,p1,p2,...)`（原来跳过了 P2） |
| `intersect_vec_plane(V0,V1,P0,P1,P3,...)` | `(v0,v1,p0,p1,p2,...)` |
| `aLE` / `aTE`（RoundTipSection 参数） | `angle_LE` / `angle_TE` |
| `rLE_new` | `r_LE_new` |
| `foil_for_XFoil` | `foil_for_xfoil` |
| `run_xfoil(AoAs,Cls,Minf,Re,nNodes,iterVis,MReDependence)` | `(alphas,cls,minf,re,n_nodes,iter_viscous,re_mach_dependence)` |
| `naca.py`: `NACA_series`/`finite_TE`/`linspace`/`interpolate`/参数`np`/`queryPoints` | `naca_series`/`finite_te`/`_linspace`/`_interpolate`/`n_point`/`query_points`（消除了对 numpy/scipy 的名字遮蔽） |
| `extrapolate2sym` | `extrapolate_to_sym` |
| `DeflectSurf` 及属性 `LE0/LE1/TE0/TE1/AX0/AX1` | `DeflectSurface` 及 `le_0/le_1/te_0/te_1/ax_0/ax_1` |
| `Lofting_2Profile` / `Lofting_Revolution` | `Lofting2Profile` / `LoftingRevolution` |

### 9.3 同名异义消除

- **三个 section 类的 `section()` 方法 → `geo()`**：`BasicSection.section()`（2D→3D 变换）与 `Section.section(cst_u,cst_l,...)`（CST 生成+变换）签名不同却同名，且方法名与类名相同。改名后与曲面的 `Surface.geo()` 形成一致的"生成几何"词汇。
- **`math.transform` → `transform_foil`**：与 `transform_curve` 对应，不再占用过于通用的名字。
- **删除 `BasicSection` 的别名属性 `scale`（=`chord`）和 `rot_z`（=`twist`）**：一个概念只留一个名字，调用处改用 `.chord`/`.twist`。注意 GuideCurve 字典里的 `'rot_z'` 键是另一回事，未动。
- **删除 `GuideCurve.__call__`**：用 `guide('x')` 做字典访问太隐晦，内部改为显式 `guide.global_guide_curve['x']`。
- **`Lofting2Profile` 的字典属性 `guide_curve` → `guide_curve_data`**：`Lofting.guide_curve` 是 GuideCurve 对象，两者同名不同型。
- **io 写出函数的 `ID` 参数 → `zone_id`**：它身兼"追加标志"与"zone 名"两职，文档已写明 `zone_id=0` 表示新建文件。
- **xfoil 结果字典的键 `'AoAs'/'Minfs'/'CLs'` 未改**：它们是结果数据的字段名，不是 `run_xfoil` 的参数。

### 9.4 顺手修正（非改名）

- `_read_surface_setting_blocks`/`_apply_surface_layout` 从 basic.py 移到 io.py 并改为公开函数 `read_surface_setting_blocks`/`apply_surface_layout`——消除跨模块引用私有成员的问题（io.py 无包内依赖，无循环）。
- `surface.BasicSurface.__init__` 补上 `self.section_s_loc = None`（原本在方法里才创建的属性）。
- `Section.tail` 文档修正：实际存的是**相对弦长**值（read_setting 会除以弦长），原文档误写为绝对值 (m)；`RoundTipSection.tail` 才是绝对值。
- `LoftingRevolution.n_section = n_profile + 1`（旋转面首尾相接）加了注释说明，与 `GuideCurve.n_section` 的差一语义不再是个暗坑。
- `io.output_curve` 文档字符串误写为"输出翼型数据"，实为通用曲线输出，已改正。

### 9.5 有意保留的名字（判断依据）

- `xLE/yLE/zLE` 及 LE/TE 后缀：全库一致的几何约定，文档有说明。
- `cst_foil/cst_curve/cst_u/cst_l` 等 `cst_` 前缀：CST 方法的标准术语。
- `nn`/`ns`：点数约定的简写，200+ 处使用，文档有说明；`n_spanwise` 保留在 GuideCurve/Lofting 家族内部。
- `geo()`：所有示例和 README 的核心方法，改名收益小而破坏性大。
- `dist_clustcos`：Kulfan CST 文献中的标准函数名。
- `tools/blwf.py` 的大写属性（`FNS/ZLE/CHORD/ITH/...`）：与 BLWF 求解器输入文件的令牌一一对应，改名反而会丢失对照关系。
- `FoilGeoFeatures`/`FoilModification`：`Geo` 是无歧义的常见缩写，且在教程中大量使用。
- `example/blade/test-fromCylinder.py` 文件名未改（只改标识符，不改文件名）。

### 9.6 第二轮验证

- 全部 32 个示例重跑：28 通过（与改名前完全相同），4 个外部依赖失败不变；111 个基线输出文件**零差异**。
- Sphinx 全新构建通过，警告数与改名前一致（10 条预存在的 `n_point` 歧义警告）。
- 改名规模：57 个文件，+1672/-1700 行。
