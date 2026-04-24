#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import parmed


NONSTD = ["AN1", "FE1", "HD1", "HD2", "HD3", "IE1", "O11", "LIG"]


def _read_cryst1_box(pdb_path: Path) -> list[float]:
    with pdb_path.open() as handle:
        for line in handle:
            if line.startswith("CRYST1"):
                a = float(line[6:15])
                b = float(line[15:24])
                c = float(line[24:33])
                alpha = float(line[33:40])
                beta = float(line[40:47])
                gamma = float(line[47:54])
                return [a, b, c, alpha, beta, gamma]
    raise ValueError(f"No CRYST1 record found in {pdb_path}")


def _read_inpcrd_box(inpcrd_path: Path) -> list[float]:
    lines = inpcrd_path.read_text().splitlines()
    if len(lines) < 3:
        raise ValueError(f"{inpcrd_path} is too short to contain box dimensions")
    parts = lines[-1].split()
    if len(parts) != 6:
        raise ValueError(f"Last line of {inpcrd_path} does not contain 6 box values")
    return [float(x) for x in parts]


def _resolve_box(src_pdb: Path, inpcrd: Path) -> list[float]:
    try:
        return _read_cryst1_box(src_pdb)
    except ValueError:
        return _read_inpcrd_box(inpcrd)


def _assert_periodic_prmtop(prmtop_path: Path) -> None:
    lines = prmtop_path.read_text().splitlines()
    in_pointers = False
    pointer_values: list[int] = []
    saw_box_dimensions = any(line.startswith("%FLAG BOX_DIMENSIONS") for line in lines)

    for line in lines:
        if line.startswith("%FLAG "):
            flag = line[6:].strip()
            in_pointers = flag == "POINTERS"
            continue
        if line.startswith("%FORMAT"):
            continue
        if in_pointers:
            for i in range(0, len(line), 8):
                chunk = line[i:i + 8].strip()
                if chunk:
                    pointer_values.append(int(chunk))
            if len(pointer_values) >= 28:
                break

    if len(pointer_values) < 28:
        raise RuntimeError(f"Failed to parse POINTERS from {prmtop_path}")
    if pointer_values[27] == 0:
        raise RuntimeError(f"{prmtop_path} still has IFBOX=0 after box rewrite")
    if not saw_box_dimensions:
        raise RuntimeError(f"{prmtop_path} is missing BOX_DIMENSIONS after box rewrite")


def _rewrite_amber_files_with_box(prmtop: Path, inpcrd: Path, box: list[float]) -> None:
    struct = parmed.load_file(str(prmtop), xyz=str(inpcrd))
    struct.box = box
    struct.save(str(prmtop), overwrite=True)
    struct.save(str(inpcrd), overwrite=True)
    _assert_periodic_prmtop(prmtop)


def _normalize_pdb_for_tleap(src_pdb: Path, dst_pdb: Path) -> None:
    """Normalize names to improve tleap template matching."""
    lines = src_pdb.read_text().splitlines()

    atom_records: list[str] = []
    for ln in lines:
        if ln[:6].strip() in {"ATOM", "HETATM"}:
            atom_records.append(ln)

    residue_atoms: dict[tuple[str, str, str], set[str]] = {}
    for ln in atom_records:
        key = (ln[21:22], ln[22:26], ln[17:20])
        if key not in residue_atoms:
            residue_atoms[key] = set()
        residue_atoms[key].add(ln[12:16].strip())

    out: list[str] = []
    for ln in lines:
        rec = ln[:6].strip()
        if rec in {"ATOM", "HETATM"}:
            chain = ln[21:22]
            resid = ln[22:26]
            resn = ln[17:20].strip()
            atom_name = ln[12:16].strip()
            # Ions
            if resn == "NA" and atom_name != "NA":
                ln = f"{ln[:12]}{'NA':>4}{ln[16:]}"
                atom_name = "NA"
            elif resn == "CL" and atom_name != "CL":
                ln = f"{ln[:12]}{'CL':>4}{ln[16:]}"
                atom_name = "CL"

            # Histidine protonation: HIS -> HID/HIE/HIP by side-chain hydrogens.
            if resn == "HIS":
                atoms_here = residue_atoms.get((chain, resid, "HIS"), set())
                has_hd1 = "HD1" in atoms_here
                has_he2 = "HE2" in atoms_here
                new_resn = "HIS"
                if has_hd1 and has_he2:
                    new_resn = "HIP"
                elif has_hd1:
                    new_resn = "HID"
                elif has_he2:
                    new_resn = "HIE"
                if new_resn != "HIS":
                    ln = f"{ln[:17]}{new_resn:>3}{ln[20:]}"

            # Common N-terminus naming mismatch in exported PDBs: H/H2/H3 -> H1/H2/H3
            atoms_here = residue_atoms.get((chain, resid, ln[17:20]), set())
            if atom_name == "H" and "H2" in atoms_here and "H3" in atoms_here:
                ln = f"{ln[:12]}{'H1':>4}{ln[16:]}"

        out.append(ln)

    dst_pdb.write_text("\n".join(out) + "\n")


def run_tleap(reorg_md_dir: Path, out_dir: Path) -> tuple[Path, Path, Path]:
    ncaa_dir = reorg_md_dir / "pdbs" / "ncaa_lib_14"
    src_pdb = reorg_md_dir / "reorg_runs" / "LIG" / "complex" / "system.pdb"

    if not src_pdb.exists():
        raise FileNotFoundError(f"Missing source PDB: {src_pdb}")

    out_dir.mkdir(parents=True, exist_ok=True)
    tleap_in = out_dir / "tleap_build_system.in"
    tleap_log = out_dir / "tleap_build_system.log"
    tleap_input_pdb = out_dir / "system_tleap_input.pdb"
    out_prmtop = out_dir / "system_tleap.prmtop"
    out_inpcrd = out_dir / "system_tleap.inpcrd"

    _normalize_pdb_for_tleap(src_pdb, tleap_input_pdb)

    lines: list[str] = [
        "source leaprc.protein.ff14SB",
        "source leaprc.gaff2",
        "source leaprc.water.tip3p",
        "set default nocenter on",
    ]

    for resn in NONSTD:
        frcmod = ncaa_dir / f"{resn}_AM1BCC-GAFF.frcmod"
        mol2 = ncaa_dir / f"{resn}_AM1BCC-GAFF.mol2"
        if not frcmod.exists() or not mol2.exists():
            raise FileNotFoundError(f"Missing ncaa_lib_14 params for {resn}: {frcmod} / {mol2}")
        lines.append(f"loadamberparams {frcmod.resolve()}")
        lines.append(f"{resn} = loadmol2 {mol2.resolve()}")

    lines.extend(
        [
            f"sys = loadpdb {tleap_input_pdb.resolve()}",
            f"saveamberparm sys {out_prmtop.resolve()} {out_inpcrd.resolve()}",
            "quit",
        ]
    )
    tleap_in.write_text("\n".join(lines) + "\n")

    proc = subprocess.run(
        ["tleap", "-f", str(tleap_in)],
        cwd=out_dir,
        text=True,
        capture_output=True,
    )
    tleap_log.write_text((proc.stdout or "") + "\n" + (proc.stderr or ""))

    if proc.returncode != 0:
        raise RuntimeError(f"tleap failed, see: {tleap_log}")
    if not out_prmtop.exists() or not out_inpcrd.exists():
        raise RuntimeError(f"tleap finished but outputs missing, see: {tleap_log}")

    box = _resolve_box(src_pdb, out_inpcrd)
    _rewrite_amber_files_with_box(out_prmtop, out_inpcrd, box)

    return out_prmtop, out_inpcrd, tleap_log


def apply_hmr(
    reorg_md_dir: Path,
    in_prmtop: Path,
    in_inpcrd: Path,
    out_dir: Path,
    hydrogen_mass: float = 3.024,
    do_hmr_water: bool = False,
) -> tuple[Path, Path]:
    repo_root = reorg_md_dir.parents[1]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

    from easybfe.amber.prep_utils import hydrogen_mass_repartition

    struct = parmed.load_file(str(in_prmtop), xyz=str(in_inpcrd))
    hydrogen_mass_repartition(
        struct,
        hydrogen_mass=hydrogen_mass,
        dowater=do_hmr_water,
    )

    out_prmtop = out_dir / "system_tleap_hmr.prmtop"
    out_inpcrd = out_dir / "system_tleap_hmr.inpcrd"
    struct.save(str(out_prmtop), overwrite=True)
    struct.save(str(out_inpcrd), overwrite=True)
    return out_prmtop, out_inpcrd


def main() -> int:
    parser = argparse.ArgumentParser(description="Build tleap system.prmtop/inpcrd for reorg_md target system")
    parser.add_argument(
        "--reorg-md-dir",
        type=Path,
        default=Path("/home/shaoq1/bin/easybfe/examples/reorg_md"),
        help="Path to the original reorg_md example directory",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "outputs",
        help="Output directory for tleap files",
    )
    parser.add_argument(
        "--apply-hmr",
        action="store_true",
        default=True,
        help="Apply easybfe-style HMR to tleap outputs and write system_tleap_hmr.*",
    )
    parser.add_argument(
        "--no-apply-hmr",
        dest="apply_hmr",
        action="store_false",
        help="Skip post-tleap HMR",
    )
    parser.add_argument(
        "--hydrogen-mass",
        type=float,
        default=3.024,
        help="Target hydrogen mass for HMR",
    )
    parser.add_argument(
        "--do-hmr-water",
        action="store_true",
        help="Also apply HMR to water hydrogens",
    )
    args = parser.parse_args()

    prmtop, inpcrd, log = run_tleap(args.reorg_md_dir.resolve(), args.out_dir.resolve())
    print(f"Wrote: {prmtop}")
    print(f"Wrote: {inpcrd}")
    print(f"Log:   {log}")

    if args.apply_hmr:
        hmr_prmtop, hmr_inpcrd = apply_hmr(
            reorg_md_dir=args.reorg_md_dir.resolve(),
            in_prmtop=prmtop,
            in_inpcrd=inpcrd,
            out_dir=args.out_dir.resolve(),
            hydrogen_mass=args.hydrogen_mass,
            do_hmr_water=args.do_hmr_water,
        )
        print(f"Wrote: {hmr_prmtop}")
        print(f"Wrote: {hmr_inpcrd}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
