#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

DEFAULT_SOLVENT_RESNAMES = {
    "WAT", "HOH", "SOL", "TIP3", "TIP3P", "OPC", "OPC3", "SPC", "SPCE",
    "NA", "Na+", "K", "K+", "CL", "Cl-", "MG", "MG2", "MG2+", "CA", "CA2", "CA2+", "ZN", "ZN2", "ZN2+",
}

HIS_EQUIV = {"HIS", "HID", "HIE", "HIP"}


@dataclass
class TopAtom:
    index: int
    atom_name: str
    atom_type: str
    charge: float
    mass: float
    residue_name: str
    residue_index: int


def _chunks(line: str, width: int) -> List[str]:
    text = line.rstrip("\n")
    return [text[i : i + width] for i in range(0, len(text), width)]


def _read_flags(prmtop_path: Path) -> Dict[str, List[str]]:
    flags: Dict[str, List[str]] = {}
    cur: Optional[str] = None
    with prmtop_path.open("r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\n")
            if line.startswith("%FLAG "):
                cur = line[6:].strip()
                flags[cur] = []
                continue
            if line.startswith("%FORMAT"):
                continue
            if cur is not None:
                flags[cur].append(line)
    return flags


def _a4(lines: List[str]) -> List[str]:
    out: List[str] = []
    for ln in lines:
        for c in _chunks(ln, 4):
            s = c.strip()
            if s:
                out.append(s)
    return out


def _i8(lines: List[str]) -> List[int]:
    out: List[int] = []
    for ln in lines:
        for c in _chunks(ln, 8):
            s = c.strip()
            if s:
                out.append(int(s))
    return out


def _e16(lines: List[str]) -> List[float]:
    out: List[float] = []
    for ln in lines:
        for c in _chunks(ln, 16):
            s = c.strip()
            if s:
                out.append(float(s.replace("D", "E")))
    return out


def _atoms(prmtop_path: Path) -> List[TopAtom]:
    flags = _read_flags(prmtop_path)
    natom = _i8(flags.get("POINTERS", []))[0]
    atom_names = _a4(flags.get("ATOM_NAME", []))
    atom_types = _a4(flags.get("AMBER_ATOM_TYPE", []))
    charges = _e16(flags.get("CHARGE", []))
    masses = _e16(flags.get("MASS", []))
    resnames = _a4(flags.get("RESIDUE_LABEL", []))
    resptr = _i8(flags.get("RESIDUE_POINTER", []))

    if not (len(atom_names) == len(atom_types) == len(charges) == len(masses) == natom):
        raise ValueError("Invalid prmtop atom arrays")

    out: List[TopAtom] = []
    ridx = 0
    for i in range(1, natom + 1):
        while ridx + 1 < len(resptr) and i >= resptr[ridx + 1]:
            ridx += 1
        out.append(
            TopAtom(
                index=i,
                atom_name=atom_names[i - 1],
                atom_type=atom_types[i - 1],
                charge=charges[i - 1] / 18.2223,
                mass=masses[i - 1],
                residue_name=resnames[ridx],
                residue_index=ridx + 1,
            )
        )
    return out


def _sig(a: TopAtom, strict_type: bool) -> Tuple[str, ...]:
    residue_name = "HIS" if a.residue_name in HIS_EQUIV else a.residue_name
    if strict_type:
        return (residue_name, a.atom_name, a.atom_type, str(a.residue_index))
    return (residue_name, a.atom_name, str(a.residue_index))


def compare(
    prmtop_a: Path,
    prmtop_b: Path,
    output_csv: Optional[Path],
    max_print: int,
    exclude_solvent: bool,
    solvent_resnames: Sequence[str],
    strict_type: bool,
) -> int:
    all_a = _atoms(prmtop_a)
    all_b = _atoms(prmtop_b)
    solvent = set(solvent_resnames)

    atoms_a = [a for a in all_a if a.residue_name not in solvent] if exclude_solvent else all_a
    atoms_b = [b for b in all_b if b.residue_name not in solvent] if exclude_solvent else all_b

    print(f"A: {prmtop_a}")
    print(f"B: {prmtop_b}")
    print(f"NATOM(A)={len(atoms_a)}, NATOM(B)={len(atoms_b)}")
    if exclude_solvent:
        print(f"Filtered non-solvent atoms: A={len(atoms_a)} B={len(atoms_b)}")

    min_n = min(len(atoms_a), len(atoms_b))
    rows: List[Dict[str, object]] = []
    for i in range(min_n):
        a = atoms_a[i]
        b = atoms_b[i]
        if _sig(a, strict_type=strict_type) != _sig(b, strict_type=strict_type):
            rows.append(
                {
                    "atom_index": i + 1,
                    "a_atom_index": a.index,
                    "b_atom_index": b.index,
                    "a_resname": a.residue_name,
                    "a_resid": a.residue_index,
                    "a_atom": a.atom_name,
                    "a_type": a.atom_type,
                    "a_charge": a.charge,
                    "a_mass": a.mass,
                    "b_resname": b.residue_name,
                    "b_resid": b.residue_index,
                    "b_atom": b.atom_name,
                    "b_type": b.atom_type,
                    "b_charge": b.charge,
                    "b_mass": b.mass,
                }
            )

    print(f"Total mismatched atom positions: {len(rows)}")
    if rows:
        print("First mismatches:")
        for row in rows[:max_print]:
            print(
                f"#{row['atom_index']} (A_idx={row['a_atom_index']}, B_idx={row['b_atom_index']}): "
                f"A=({row['a_resname']}{row['a_resid']}:{row['a_atom']}/{row['a_type']}) "
                f"B=({row['b_resname']}{row['b_resid']}:{row['b_atom']}/{row['b_type']})"
            )
    else:
        print("Atom order is fully identical.")

    if output_csv is not None:
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        with output_csv.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=[
                    "atom_index", "a_atom_index", "b_atom_index",
                    "a_resname", "a_resid", "a_atom", "a_type", "a_charge", "a_mass",
                    "b_resname", "b_resid", "b_atom", "b_type", "b_charge", "b_mass",
                ],
            )
            writer.writeheader()
            writer.writerows(rows)
        print(f"Mismatch CSV written: {output_csv}")

    return 0 if not rows else 1


def main() -> int:
    ap = argparse.ArgumentParser(description="Compare atom order between two prmtop files")
    ap.add_argument("prmtop_a", type=Path)
    ap.add_argument("prmtop_b", type=Path)
    ap.add_argument("--output-csv", type=Path, default=Path("prmtop_atom_order_mismatch.csv"))
    ap.add_argument("--max-print", type=int, default=20)
    ap.add_argument("--exclude-solvent", action="store_true")
    ap.add_argument("--solvent-resname", action="append", default=[])
    ap.add_argument("--strict-type", action="store_true", help="Include atom type in order signature")
    args = ap.parse_args()

    solvent = set(DEFAULT_SOLVENT_RESNAMES)
    solvent.update(args.solvent_resname)

    return compare(
        prmtop_a=args.prmtop_a.resolve(),
        prmtop_b=args.prmtop_b.resolve(),
        output_csv=args.output_csv.resolve() if args.output_csv else None,
        max_print=max(1, args.max_print),
        exclude_solvent=args.exclude_solvent,
        solvent_resnames=sorted(solvent),
        strict_type=args.strict_type,
    )


if __name__ == "__main__":
    raise SystemExit(main())
