#!/usr/bin/env python3
"""Compare atom ordering between two Amber prmtop files.

Reports per-atom mismatches and can export to CSV.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


DEFAULT_SOLVENT_RESNAMES = {
    "WAT",
    "HOH",
    "SOL",
    "TIP3",
    "TIP3P",
    "OPC",
    "OPC3",
    "SPC",
    "SPCE",
    "NA",
    "Na+",
    "K",
    "K+",
    "CL",
    "Cl-",
    "MG",
    "MG2",
    "MG2+",
    "CA",
    "CA2",
    "CA2+",
    "ZN",
    "ZN2",
    "ZN2+",
}


@dataclass
class TopAtom:
    index: int  # 1-based
    atom_name: str
    atom_type: str
    charge: float
    mass: float
    residue_name: str
    residue_index: int  # 1-based


def _parse_fixed_width_chunks(line: str, width: int) -> List[str]:
    out: List[str] = []
    text = line.rstrip("\n")
    for i in range(0, len(text), width):
        out.append(text[i : i + width])
    return out


def _read_prmtop_flags(prmtop_path: Path) -> Dict[str, List[str]]:
    flags: Dict[str, List[str]] = {}
    current_flag: Optional[str] = None
    with prmtop_path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("%FLAG "):
                current_flag = line[6:].strip()
                flags[current_flag] = []
                continue
            if line.startswith("%FORMAT"):
                continue
            if current_flag is not None:
                flags[current_flag].append(line)
    return flags


def _parse_a4(lines: List[str]) -> List[str]:
    vals: List[str] = []
    for ln in lines:
        for chunk in _parse_fixed_width_chunks(ln, 4):
            s = chunk.strip()
            if s:
                vals.append(s)
    return vals


def _parse_i8(lines: List[str]) -> List[int]:
    vals: List[int] = []
    for ln in lines:
        for chunk in _parse_fixed_width_chunks(ln, 8):
            s = chunk.strip()
            if s:
                vals.append(int(s))
    return vals


def _parse_e16(lines: List[str]) -> List[float]:
    vals: List[float] = []
    for ln in lines:
        for chunk in _parse_fixed_width_chunks(ln, 16):
            s = chunk.strip()
            if s:
                vals.append(float(s.replace("D", "E")))
    return vals


def _build_atom_table(prmtop_path: Path) -> List[TopAtom]:
    flags = _read_prmtop_flags(prmtop_path)

    pointers = _parse_i8(flags.get("POINTERS", []))
    if not pointers:
        raise ValueError(f"Failed to parse POINTERS from {prmtop_path}")
    natom = pointers[0]

    atom_names = _parse_a4(flags.get("ATOM_NAME", []))
    atom_types = _parse_a4(flags.get("AMBER_ATOM_TYPE", []))
    charges_raw = _parse_e16(flags.get("CHARGE", []))
    masses = _parse_e16(flags.get("MASS", []))
    residue_labels = _parse_a4(flags.get("RESIDUE_LABEL", []))
    residue_pointer = _parse_i8(flags.get("RESIDUE_POINTER", []))

    if len(atom_names) != natom:
        raise ValueError(f"ATOM_NAME length {len(atom_names)} != NATOM {natom}")
    if len(atom_types) != natom:
        raise ValueError(f"AMBER_ATOM_TYPE length {len(atom_types)} != NATOM {natom}")
    if len(charges_raw) != natom:
        raise ValueError(f"CHARGE length {len(charges_raw)} != NATOM {natom}")
    if len(masses) != natom:
        raise ValueError(f"MASS length {len(masses)} != NATOM {natom}")
    if not residue_labels or not residue_pointer:
        raise ValueError("Failed to parse RESIDUE_LABEL / RESIDUE_POINTER")

    atoms: List[TopAtom] = []
    res_idx = 0
    for atom_idx in range(1, natom + 1):
        while res_idx + 1 < len(residue_pointer) and atom_idx >= residue_pointer[res_idx + 1]:
            res_idx += 1
        charge_e = charges_raw[atom_idx - 1] / 18.2223
        atoms.append(
            TopAtom(
                index=atom_idx,
                atom_name=atom_names[atom_idx - 1],
                atom_type=atom_types[atom_idx - 1],
                charge=charge_e,
                mass=masses[atom_idx - 1],
                residue_name=residue_labels[res_idx],
                residue_index=res_idx + 1,
            )
        )
    return atoms


def _atom_signature(atom: TopAtom) -> Tuple[str, str, str, int]:
    return (atom.residue_name, atom.atom_name, atom.atom_type, atom.residue_index)


def compare_atom_order(
    prmtop_a: Path,
    prmtop_b: Path,
    exclude_solvent: bool = False,
    solvent_resnames: Optional[Sequence[str]] = None,
    output_csv: Optional[Path] = None,
    max_print: int = 20,
) -> int:
    all_atoms_a = _build_atom_table(prmtop_a)
    all_atoms_b = _build_atom_table(prmtop_b)

    solvent_set = set(solvent_resnames or [])
    if exclude_solvent:
        atoms_a = [a for a in all_atoms_a if a.residue_name not in solvent_set]
        atoms_b = [b for b in all_atoms_b if b.residue_name not in solvent_set]
    else:
        atoms_a = all_atoms_a
        atoms_b = all_atoms_b

    natom_a = len(atoms_a)
    natom_b = len(atoms_b)
    print(f"A: {prmtop_a}")
    print(f"B: {prmtop_b}")
    print(f"NATOM(A)={natom_a}, NATOM(B)={natom_b}")
    if exclude_solvent:
        print(
            "Filtered non-solvent atoms: "
            f"A={len(atoms_a)} (excluded {len(all_atoms_a) - len(atoms_a)}), "
            f"B={len(atoms_b)} (excluded {len(all_atoms_b) - len(atoms_b)})"
        )
        print(f"Excluded residue names: {','.join(sorted(solvent_set))}")

    min_n = min(natom_a, natom_b)
    mismatches: List[Dict[str, object]] = []
    for i in range(min_n):
        a = atoms_a[i]
        b = atoms_b[i]
        if _atom_signature(a) != _atom_signature(b):
            mismatches.append(
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

    if natom_a != natom_b:
        longer, key = (atoms_a, "a") if natom_a > natom_b else (atoms_b, "b")
        for i in range(min_n, len(longer)):
            atom = longer[i]
            row = {
                "atom_index": i + 1,
                "a_atom_index": "",
                "b_atom_index": "",
                "a_resname": "",
                "a_resid": "",
                "a_atom": "",
                "a_type": "",
                "a_charge": "",
                "a_mass": "",
                "b_resname": "",
                "b_resid": "",
                "b_atom": "",
                "b_type": "",
                "b_charge": "",
                "b_mass": "",
            }
            row[f"{key}_resname"] = atom.residue_name
            row[f"{key}_resid"] = atom.residue_index
            row[f"{key}_atom"] = atom.atom_name
            row[f"{key}_type"] = atom.atom_type
            row[f"{key}_charge"] = atom.charge
            row[f"{key}_mass"] = atom.mass
            row[f"{key}_atom_index"] = atom.index
            mismatches.append(row)

    print(f"Total mismatched atom positions: {len(mismatches)}")
    if not mismatches:
        print("Atom order is fully identical.")
        return 0

    print("\nFirst mismatches:")
    for row in mismatches[:max_print]:
        idx = row["atom_index"]
        print(
            f"#{idx} (A_idx={row['a_atom_index']}, B_idx={row['b_atom_index']}): "
            f"A=({row['a_resname']}{row['a_resid']}:{row['a_atom']}/{row['a_type']}) "
            f"B=({row['b_resname']}{row['b_resid']}:{row['b_atom']}/{row['b_type']})"
        )

    if output_csv is not None:
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        with output_csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "atom_index",
                    "a_atom_index",
                    "b_atom_index",
                    "a_resname",
                    "a_resid",
                    "a_atom",
                    "a_type",
                    "a_charge",
                    "a_mass",
                    "b_resname",
                    "b_resid",
                    "b_atom",
                    "b_type",
                    "b_charge",
                    "b_mass",
                ],
            )
            writer.writeheader()
            writer.writerows(mismatches)
        print(f"\nMismatch CSV written: {output_csv}")

    return 1


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Check whether two prmtop files have identical atom ordering.",
    )
    parser.add_argument("prmtop_a", type=Path, help="Reference prmtop file")
    parser.add_argument("prmtop_b", type=Path, help="Target prmtop file")
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=Path("analysis/prmtop_atom_order_mismatch.csv"),
        help="CSV path for mismatched atom entries",
    )
    parser.add_argument(
        "--exclude-solvent",
        action="store_true",
        help="Compare only non-solvent atoms (exclude common water/ion residue names).",
    )
    parser.add_argument(
        "--solvent-resname",
        action="append",
        default=[],
        help="Additional residue name to exclude when --exclude-solvent is set (repeatable).",
    )
    parser.add_argument(
        "--max-print",
        type=int,
        default=20,
        help="Maximum mismatch lines printed to stdout",
    )
    args = parser.parse_args()

    solvent_resnames = set(DEFAULT_SOLVENT_RESNAMES)
    solvent_resnames.update(args.solvent_resname)

    return compare_atom_order(
        prmtop_a=args.prmtop_a.resolve(),
        prmtop_b=args.prmtop_b.resolve(),
        exclude_solvent=args.exclude_solvent,
        solvent_resnames=sorted(solvent_resnames),
        output_csv=args.output_csv.resolve() if args.output_csv else None,
        max_print=max(1, args.max_print),
    )


if __name__ == "__main__":
    raise SystemExit(main())
