# Reorg MD (TLEaP system build) example

这个示例是 `examples/reorg_md` 的**独立替代流程**：

- 仍然使用 `reorg_md` 的目标体系（`3pzw_sub.pdb`，含 `LIG` + FE1 metal center）。
- 先用 `reorg_md/run_setup_reorg.sh` 生成 easybfe 标准工作目录。
- 再用 `tleap` + `pdbs/ncaa_lib_14/*.mol2/*.frcmod` 从 `system.pdb` 重新生成
  `system_tleap.prmtop/system_tleap.inpcrd`。
- 然后对 `tleap` 产物施加与 easybfe 相同逻辑的 HMR，得到
  `system_tleap_hmr.prmtop/system_tleap_hmr.inpcrd`。
- 最后自动核验和 easybfe 原始 `system.prmtop` 的一致性（默认比较非溶剂原子顺序）。

## 目录结构

- `build_tleap_system.py`：用 tleap 生成 `system_tleap.prmtop/inpcrd`，并自动生成 HMR 版
- `check_prmtop_atom_order.py`：比较两个 prmtop 原子顺序并导出 CSV
- `run_setup_reorg_tleap.sh`：一键执行（可选重跑 easybfe setup）

## 一键运行

```bash
cd /home/shaoq1/bin/easybfe/examples/reorg_md_tleap
./run_setup_reorg_tleap.sh
```

运行完成后关键输出：

- `outputs/system_tleap.prmtop`
- `outputs/system_tleap.inpcrd`
- `outputs/system_tleap_hmr.prmtop`
- `outputs/system_tleap_hmr.inpcrd`
- `outputs/prmtop_atom_order_mismatch_nonsolvent.csv`

## 仅做比较（不重建 easybfe 工作目录）

```bash
SKIP_EASYBFE_SETUP=1 ./run_setup_reorg_tleap.sh
```


说明：比较脚本默认把 `HIS/HID/HIE/HIP` 视为同义残基。
如果结果仍显示不一致，通常是因为 `system.pdb` 到 tleap 重建过程中
其它命名/原子类型（例如 `H`→`H1`、GAFF 类型大小写）被规范化。
请以 `outputs/prmtop_atom_order_mismatch_nonsolvent.csv` 为准定位差异。
