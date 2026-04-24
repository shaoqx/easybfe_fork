#!/usr/bin/env bash
set -euo pipefail

set +u
source ~/bin/miniconda3/bin/activate easybfe
set -u

THIS_DIR="$(cd "$(dirname "$0")" && pwd)"
REORG_MD_DIR="/home/shaoq1/bin/easybfe/examples/reorg_md"
OUT_DIR="$THIS_DIR/outputs"
mkdir -p "$OUT_DIR"

if [[ "${SKIP_EASYBFE_SETUP:-0}" != "1" ]]; then
  echo "[1/4] Running easybfe baseline setup in $REORG_MD_DIR"
  (cd "$REORG_MD_DIR" && ./run_setup_reorg.sh)
else
  echo "[1/4] SKIP_EASYBFE_SETUP=1, skip baseline setup"
fi

echo "[2/4] Build tleap system prmtop/inpcrd"
python "$THIS_DIR/build_tleap_system.py" \
  --reorg-md-dir "$REORG_MD_DIR" \
  --out-dir "$OUT_DIR"

REF_PRMTOP="$REORG_MD_DIR/reorg_runs/LIG/complex/system.prmtop"
REF_INPCRD="$REORG_MD_DIR/reorg_runs/LIG/complex/system.inpcrd"
RAW_PRMTOP="$OUT_DIR/system_tleap.prmtop"
RAW_INPCRD="$OUT_DIR/system_tleap.inpcrd"
NEW_PRMTOP="$OUT_DIR/system_tleap_hmr.prmtop"
NEW_INPCRD="$OUT_DIR/system_tleap_hmr.inpcrd"

echo "[3/4] Quick NATOM cross-check with inpcrd"
python - <<PY
from pathlib import Path

def natom_prmtop(path: Path) -> int:
    lines = path.read_text().splitlines()
    cur = None
    vals = []
    for ln in lines:
        if ln.startswith('%FLAG '):
            cur = ln[6:].strip()
            continue
        if ln.startswith('%FORMAT'):
            continue
        if cur == 'POINTERS':
            for i in range(0, len(ln), 8):
                s = ln[i:i+8].strip()
                if s:
                    vals.append(int(s))
            if vals:
                return vals[0]
    raise RuntimeError('POINTERS not found')

def natom_inpcrd(path: Path) -> int:
    with path.open() as f:
        f.readline()
        return int(f.readline().split()[0])

pairs = [
    ('ref', Path('$REF_PRMTOP'), Path('$REF_INPCRD')),
    ('tleap_raw', Path('$RAW_PRMTOP'), Path('$RAW_INPCRD')),
    ('tleap_hmr', Path('$NEW_PRMTOP'), Path('$NEW_INPCRD')),
]
for tag, ptop, crd in pairs:
    np = natom_prmtop(ptop)
    nc = natom_inpcrd(crd)
    print(f"{tag}: NATOM(prmtop)={np}, NATOM(inpcrd)={nc}")
    if np != nc:
        raise SystemExit(f"NATOM mismatch for {tag}")
PY

echo "[3b/4] Hydrogen-mass sanity check"
python - <<PY
from pathlib import Path

def read_flags(path):
    flags = {}
    cur = None
    for raw in Path(path).read_text().splitlines():
        if raw.startswith('%FLAG '):
            cur = raw[6:].strip()
            flags[cur] = []
            continue
        if raw.startswith('%FORMAT'):
            continue
        if cur is not None:
            flags[cur].append(raw)
    return flags

def chunks(line, width):
    return [line[i:i+width] for i in range(0, len(line), width)]

def a4(lines):
    out = []
    for ln in lines:
        for c in chunks(ln, 4):
            s = c.strip()
            if s:
                out.append(s)
    return out

def e16(lines):
    out = []
    for ln in lines:
        for c in chunks(ln, 16):
            s = c.strip()
            if s:
                out.append(float(s.replace('D', 'E')))
    return out

for tag, path in [('ref', Path('$REF_PRMTOP')), ('tleap_hmr', Path('$NEW_PRMTOP'))]:
    flags = read_flags(path)
    names = a4(flags['ATOM_NAME'])
    masses = e16(flags['MASS'])
    hmasses = sorted({round(m, 6) for n, m in zip(names, masses) if n.startswith('H')})
    print(f"{tag}: unique H masses = {hmasses}")
PY

echo "[4/4] Compare prmtop atom order (non-solvent)"
set +e
python "$THIS_DIR/check_prmtop_atom_order.py" \
  "$REF_PRMTOP" "$NEW_PRMTOP" \
  --exclude-solvent \
  --output-csv "$OUT_DIR/prmtop_atom_order_mismatch_nonsolvent.csv" \
  --max-print 20
cmp_status=$?
set -e

if [[ $cmp_status -eq 0 ]]; then
  echo "[OK] Non-solvent atom order is identical."
else
  echo "[WARN] Non-solvent atom order is NOT identical; see CSV for details."
fi

echo
if [[ -s "$OUT_DIR/prmtop_atom_order_mismatch_nonsolvent.csv" ]]; then
  echo "Done. Check mismatch CSV: $OUT_DIR/prmtop_atom_order_mismatch_nonsolvent.csv"
fi
