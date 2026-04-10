# Reorg MD Example (Single-Leg Decoupling)

This example prepares a single-leg complex-decoupling workflow for protein
conformational reorganization estimation via reweighting/MBAR
(analysis stage is not included here).

## What this setup does

- Uses only the `complex` leg (no `solvent` or `restraint` legs).
- Disables Boresch restraints with `boresch: null`.
- Applies positional restraints to ligand heavy atoms (`ntr`) so ligand remains
  localized relative to the simulation box while interactions are decoupled.
- Keeps decoupling masks and lambda schedule from the ABFE engine.

## Files

- `config_reorg_5ns.yaml`: single-leg reorg config.
- `prepare_inputs.py`:
  - splits `pdbs/3pzw_sub.pdb` into protein/ligand helper PDBs and generates a fixed protein (`3pzw_sub_std_no_lig_fixed.pdb`) via `ProteinFixer`;
  - builds `ligands/LIG` directly from `pdbs/ncaa_lib_14/LIG_AM1BCC-GAFF.{mol2,frcmod}`
    via `tleap`, then writes `LIG.sdf`, `LIG.pdb`, `LIG.xml`, `LIG.prmtop`, `LIG.inpcrd`.
- `run_setup_reorg.sh`: one-command setup entry.

## How to run

```bash
cd /home/shaoq1/bin/easybfe/examples/reorg_md
./run_setup_reorg.sh
```

## Note on other nonstandard residues

This example uses `3pzw_sub_std_no_lig_fixed.pdb` (standard protein residues only, termini fixed), so
additional NCAA force fields are not required for the protein in this workflow.
If you later switch to the full protein PDB that includes nonstandard residues,
you will need to convert their Amber params to OpenMM XML and add them as
`extra_ff` entries in the ABFE config.
