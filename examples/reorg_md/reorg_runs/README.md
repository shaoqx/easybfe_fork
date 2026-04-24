# `reorg_runs`：推荐使用的“easybfe 建目录 + tleap 修正拓扑 + 回灌”流程

这份说明描述的不是“原始 easybfe setup 流程本身”，而是当前这个体系**推荐实际采用的完整工作流**：

1. 先让 `easybfe` 生成完整的 ABFE/reorg 工作目录
2. 再用 `tleap` 基于用户提供的 Amber 参数重新生成**正确的系统拓扑文件**
3. 将新的 `prmtop/inpcrd` 放回 `easybfe` 生成的工作目录中
4. 从头运行 `lambda*` 目录中的 MD

这样做的原因是：

> 对当前这个含有 custom amino acids / 金属中心 bonded model 的体系，直接由目前 easybfe 流程生成的 `system.prmtop` 会缺失某些非常规氨基酸与主链相连时的 **peptide bond 相关参数**；因此更稳妥的做法是保留 easybfe 生成的工作目录结构和 alchemical 输入，但把系统拓扑替换成 `tleap` 重建并经 HMR 处理后的版本。

---

## 1. 用户原始输入

用户需要提供两类输入：

### 1.1 复合物 PDB

- `examples/reorg_md/pdbs/3pzw_sub.pdb`

这个文件中包含：

- 蛋白
- 配体 `LIG`
- 溶剂和离子
- 金属中心 `FE1`
- 与金属中心/蛋白主链相连的非常规氨基酸残基

### 1.2 非常规组分的 Amber 参数

位于：

- `examples/reorg_md/pdbs/ncaa_lib_14/*.mol2`
- `examples/reorg_md/pdbs/ncaa_lib_14/*.frcmod`

当前实际涉及的非常规组分包括：

- `LIG`
- `AN1`
- `FE1`
- `HD1`
- `HD2`
- `HD3`
- `IE1`
- `O11`

这些文件是后续 `tleap` 修正系统拓扑的核心依据。

---

## 2. 第一步：先让 easybfe 生成完整工作目录

首先运行：

```bash
cd /home/shaoq1/bin/easybfe/examples/reorg_md
./run_setup_reorg.sh
```

这一步的作用不是“得到最终可信的拓扑”，而是**先得到 easybfe 擅长生成的整套工作目录骨架**。

### 2.1 这一步具体产出什么

`run_setup_reorg.sh` 会做两件事：

1. 运行 `prepare_inputs.py`
2. 运行 `easybfe cli abfe setup ./config_reorg_5ns.yaml`

最终生成：

- `examples/reorg_md/reorg_runs/LIG/complex/system.prmtop`
- `examples/reorg_md/reorg_runs/LIG/complex/system.inpcrd`
- `examples/reorg_md/reorg_runs/LIG/complex/system.pdb`
- `examples/reorg_md/reorg_runs/LIG/complex/lambda0 ... lambda15`
- 每个 `lambdaN` 下的
  - `01.em`
  - `02.heat`
  - `03.pres`
  - `04.pre_prod`
  - `05.prod`
  - 以及对应 `run.sh` / `run.submit`

### 2.2 为什么仍然要先跑这一步

虽然当前 easybfe 直接生成的 `system.prmtop` 不足以作为最终拓扑，但这一步仍然非常重要，因为它会正确生成：

- 整个 `lambda0`–`lambda15` 目录树
- 每个阶段的 `*.in` 文件
- MBAR / lambda schedule
- `timask1/scmask1` 等 alchemical mask
- ligand positional restraint 等运行设置

也就是说：

> 我们需要 easybfe 生成**工作目录和 alchemical 输入脚本**，但不完全信任它生成的最终 `system.prmtop`。

---

## 3. 第二步：为什么要用 tleap 重建拓扑

当前体系的难点在于：

- 蛋白中包含 custom amino acids
- 这些残基既参与蛋白主链，也参与金属中心 bonded model
- 仅依赖当前 easybfe 的 `extra_ff` / residue template 路径，最终导出的 `system.prmtop` 中可能会丢失某些 **custom residue 与相邻主链之间的 peptide bond 相关参数**

因此需要换一种策略：

### 核心思想

- 让 easybfe 负责：
  - 目录组织
  - lambda 输入
  - alchemical 设置
  - MD workflow 脚本
- 让 `tleap` 负责：
  - 最终 `system.prmtop`
  - 最终 `system.inpcrd`

这样可以同时保留：

- easybfe 在 workflow 自动生成上的便利
- Amber/tleap 在复杂 bonded 拓扑落地上的可靠性

---

## 4. 第三步：基于 easybfe 生成的 `system.pdb` 用 tleap 重建系统拓扑

完成第 1 步后，不直接使用 `reorg_runs/LIG/complex/system.prmtop`，而是使用下面这套替代流程：

```bash
cd /home/shaoq1/bin/easybfe/examples/reorg_md_tleap
SKIP_EASYBFE_SETUP=1 ./run_setup_reorg_tleap.sh
```

这一步会把 `reorg_runs/LIG/complex/system.pdb` 当作输入，用 `tleap` 重新建系统。

### 4.1 为什么以 `system.pdb` 为起点

因为 `system.pdb` 已经是 easybfe 建好的：

- ligand 在最前面
- protein / ncaa / metal center / solvent / ions 顺序固定
- 盒子、水、离子都已经放好

因此只要基于这份 `system.pdb` 重新用 `tleap` 建 `prmtop/inpcrd`，就能最大程度保留与 easybfe 工作目录的一致性。

### 4.2 `reorg_md_tleap` 做了什么

`examples/reorg_md_tleap/build_tleap_system.py` 会：

1. 读取 `reorg_runs/LIG/complex/system.pdb`
2. 对 PDB 做必要规范化
   - 规范 `NA` / `CL` 离子原子名
   - 将 `HIS` 识别并映射为 `HID/HIE/HIP`
   - 规范某些 N 端 `H/H1` 命名差异
3. 载入 `ncaa_lib_14` 中所有需要的 `mol2 + frcmod`
4. 用 `tleap` 生成：
   - `outputs/system_tleap.prmtop`
   - `outputs/system_tleap.inpcrd`

这里的关键是：

> 这一步的系统拓扑是由 `tleap` 直接根据 Amber 参数库构建出来的，因此能够更好地保留 custom amino acids 的 peptide-bond 相关项。

---

## 5. 第四步：对 tleap 拓扑做与 easybfe 一样的 HMR

因为 `config_reorg_5ns.yaml` 中设置了：

- `do_hmr: true`
- `hydrogen_mass: 3.024`

而现有 `lambda*` 输入文件中的后期阶段使用：

- `dt = 0.004`

所以如果只用裸 `tleap` 拓扑直接替换，时间步长假设会不成立。

为了解决这个问题，`reorg_md_tleap/build_tleap_system.py` 会在 `tleap` 完成后，自动调用与 easybfe 相同逻辑的 HMR，对 `tleap` 产物再处理一次，输出：

- `outputs/system_tleap_hmr.prmtop`
- `outputs/system_tleap_hmr.inpcrd`

因此最终真正建议回灌到 easybfe 工作目录中的，不是裸 `system_tleap.*`，而是：

- `system_tleap_hmr.prmtop`
- `system_tleap_hmr.inpcrd`

---

## 6. 第五步：核验 tleap/HMR 后的拓扑能否安全接到 easybfe 工作目录

在替换回去之前，需要检查两件事。

### 6.1 原子数是否一致

脚本会检查：

- 原 easybfe `system.prmtop` vs `system.inpcrd`
- 新 `system_tleap.prmtop` vs `system_tleap.inpcrd`
- 新 `system_tleap_hmr.prmtop` vs `system_tleap_hmr.inpcrd`

确保 `NATOM` 全部一致。

### 6.2 非溶剂原子顺序是否与 easybfe 工作目录兼容

因为 easybfe 自动生成的 alchemical mask 是按**原子索引**写入 `lambda*/0*.in` 的，例如：

- `timask1 = '@1-51'`
- `scmask1 = '@1-51'`

所以必须确认重新建出的拓扑中：

- 配体仍然位于前 51 个原子
- 非溶剂部分的原子顺序基本一致

当前已经完成的核验结果是：

- 在将 `HIS/HID/HIE/HIP` 视为同义后
- 非溶剂原子顺序比较只剩 1 个差异：
  - `MET2:H` vs `MET2:H1`

这意味着：

> 对现有 easybfe 工作目录而言，`system_tleap_hmr.prmtop` 已经足够接近，可作为回灌候选拓扑。

核验结果 CSV 位于：

- `examples/reorg_md_tleap/outputs/prmtop_atom_order_mismatch_nonsolvent.csv`

---

## 7. 第六步：将 tleap/HMR 拓扑放回 easybfe 工作目录

最终要放回 `easybfe` 工作目录的是：

- `examples/reorg_md_tleap/outputs/system_tleap_hmr.prmtop`
- `examples/reorg_md_tleap/outputs/system_tleap_hmr.inpcrd`

目标位置是：

- `examples/reorg_md/reorg_runs/LIG/complex/system.prmtop`
- `examples/reorg_md/reorg_runs/LIG/complex/system.inpcrd`

推荐操作顺序：

1. 先备份原始 easybfe 拓扑
   - `system.prmtop.bak_easybfe`
   - `system.inpcrd.bak_easybfe`
2. 再把 `system_tleap_hmr.*` 复制覆盖到 `complex/system.*`
3. 保留现有 `lambda*` 目录和 `*.in/*.sh` 不变

### 为什么只替换 `system.prmtop/inpcrd`

因为 easybfe 的工作目录在运行期真正依赖的是：

- `system.prmtop`
- `system.inpcrd`
- 各 `lambda*/0*.in`
- 各阶段之间的 `rst7`

而 `ligands/*`、`ffxml/*`、`prepare_inputs.py` 等，主要只在 **setup 阶段** 有用。

所以回灌时不需要重写整套 `lambda*` 输入文件，只需要让它们引用新的系统拓扑即可。

---

## 8. 第七步：替换后如何运行

替换完成后，必须注意：

### 8.1 必须从头开始跑

不能复用旧拓扑对应的：

- `01.em.rst7`
- `02.heat.rst7`
- `03.pres.rst7`
- `04.pre_prod.rst7`
- `05.prod.rst7`

也不能从旧的某个中间阶段继续接着跑。

原因是：

> 一旦 `system.prmtop` 变了，体系哈密顿量就变了；旧 `rst7` 属于旧拓扑，不应再与新拓扑混用。

因此，替换 `system.prmtop/inpcrd` 后，应当：

- 清理各 `lambda*` 中旧的运行产物
- 从 `01.em` 开始完整重跑整套流程

### 8.2 为什么可以保留原来的 `lambda*.in`

因为现在验证表明：

- ligand 仍然在最前面的 51 个原子
- 非溶剂原子顺序与原 easybfe 拓扑高度一致

因此 easybfe 写好的：

- `timask1 = '@1-51'`
- `scmask1 = '@1-51'`
- lambda schedule
- MBAR 设置
- 各阶段控制参数

仍然可以继续使用。

---

## 9. 这套推荐流程的本质

可以把当前推荐流程概括成一句话：

> **easybfe 负责生成工作目录和 alchemical 输入；tleap 负责生成最终可信的系统拓扑；然后把 tleap/HMR 得到的 `system.prmtop/inpcrd` 回灌到 easybfe 工作目录中再正式运行。**

这样做的核心收益是：

- 保留 easybfe 自动生成 workflow 的便利
- 避免直接使用当前 easybfe 生成的、有 custom amino acid peptide-bond 参数缺失风险的 `system.prmtop`
- 同时保持与现有 `lambda*` 输入的索引兼容性
- 继续满足 HMR + `dt=0.004` 的运行前提

---

## 10. 推荐的实际执行顺序

按顺序执行时，推荐这样操作：

### Step A：生成 easybfe 工作目录

```bash
cd /home/shaoq1/bin/easybfe/examples/reorg_md
./run_setup_reorg.sh
```

### Step B：生成 tleap + HMR 拓扑并核验

```bash
cd /home/shaoq1/bin/easybfe/examples/reorg_md_tleap
SKIP_EASYBFE_SETUP=1 ./run_setup_reorg_tleap.sh
```

### Step C：备份并替换 `reorg_runs/LIG/complex/system.*`

推荐替换：

- `system_tleap_hmr.prmtop` → `reorg_runs/LIG/complex/system.prmtop`
- `system_tleap_hmr.inpcrd` → `reorg_runs/LIG/complex/system.inpcrd`

### Step D：清理旧运行产物，从头启动 MD

- 删除旧 `lambda*/**/*.rst7`、`*.out`、`*.log`、`done.tag` 等
- 从 `lambda*/run.sh` 重新开始完整运行

---

## 11. 当前推荐替换对象

对于当前这个项目，**推荐放回 easybfe 工作目录的文件是**：

- `examples/reorg_md_tleap/outputs/system_tleap_hmr.prmtop`
- `examples/reorg_md_tleap/outputs/system_tleap_hmr.inpcrd`

而不是：

- `examples/reorg_md/reorg_runs/LIG/complex/system.prmtop`（原 easybfe 版）
- `examples/reorg_md_tleap/outputs/system_tleap.prmtop`（未做 HMR 的裸 tleap 版）

因为：

- 原 easybfe 版存在 custom amino acid peptide-bond 参数缺失风险
- 裸 tleap 版不满足当前 `dt=0.004` 的 HMR 假设
- `tleap + HMR` 版同时解决了这两个问题

