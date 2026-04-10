#!/usr/bin/env bash
set -euo pipefail

source ~/bin/miniconda3/bin/activate easybfe

cd "$(dirname "$0")"
repo_root="$(cd ../.. && pwd)"
export PYTHONPATH="${repo_root}:${PYTHONPATH:-}"

cli() {
  python -c 'from easybfe.cli import main; main()' "$@"
}

# 1) Split protein/ligand PDB and build ligands/LIG from Amber mol2+frcmod.
python prepare_inputs.py

# 2) Create single-leg reorg workflow under ./reorg_runs/LIG/complex.
cli abfe setup ./config_reorg_5ns.yaml

echo "Setup finished. Submit job(s), e.g.:"
echo "  cd reorg_runs/LIG/complex && sbatch run.sh -A ... -p ... --gres=gpu:A100:4 ..."
