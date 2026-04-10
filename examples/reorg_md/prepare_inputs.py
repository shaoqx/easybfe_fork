#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import openmm.app as app
import parmed

STD_PROTEIN = {
    'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','HID','HIE','HIP',
    'ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL',
    'ASH','GLH','LYN','CYM','CYX','ACE','NME'
}


def residue_name(line: str) -> str:
    return line[17:20].strip()


def split_pdb(src: Path, out_full: Path, out_std: Path, out_lig: Path, lig_name: str = 'LIG') -> None:
    lines = src.read_text().splitlines()
    full_lines: list[str] = []
    std_lines: list[str] = []
    lig_lines: list[str] = []

    prev_full = False
    prev_std = False
    prev_lig = False

    for ln in lines:
        rec = ln[:6].strip()

        if rec in {'ATOM', 'HETATM'}:
            resn = residue_name(ln)
            is_lig = (resn == lig_name)
            is_std = (resn in STD_PROTEIN) and (not is_lig)
            is_full = not is_lig

            # If a kept std-protein segment is interrupted by a dropped residue,
            # close the segment explicitly to avoid artificial peptide bonds.
            if (not is_std) and prev_std and (len(std_lines) == 0 or std_lines[-1] != 'TER'):
                std_lines.append('TER')

            if is_lig:
                lig_lines.append(ln)
            else:
                full_lines.append(ln)
                if is_std:
                    std_lines.append(ln)

            prev_full = is_full
            prev_std = is_std
            prev_lig = is_lig
            continue

        if rec == 'TER':
            if prev_full and (len(full_lines) == 0 or full_lines[-1] != 'TER'):
                full_lines.append('TER')
            if prev_std and (len(std_lines) == 0 or std_lines[-1] != 'TER'):
                std_lines.append('TER')
            if prev_lig and (len(lig_lines) == 0 or lig_lines[-1] != 'TER'):
                lig_lines.append('TER')

            prev_full = False
            prev_std = False
            prev_lig = False
            continue

        if rec == 'END':
            continue

    def finalize(content: list[str]) -> str:
        out = content[:]
        if len(out) > 0 and out[-1] != 'TER':
            out.append('TER')
        out.append('END')
        return '\n'.join(out) + '\n'

    out_full.write_text(finalize(full_lines))
    out_std.write_text(finalize(std_lines))
    out_lig.write_text(finalize(lig_lines))


def convert_lig_mol2_to_sdf(mol2_path: Path, sdf_path: Path) -> None:
    obabel = shutil.which('obabel')
    if obabel is None:
        raise RuntimeError('obabel not found in PATH')
    proc = subprocess.run([obabel, str(mol2_path), '-O', str(sdf_path)], capture_output=True, text=True)
    if proc.returncode != 0 or (not sdf_path.exists()) or sdf_path.stat().st_size == 0:
        msg = (proc.stderr or proc.stdout or '').strip()
        raise RuntimeError(f'Failed to convert mol2 to sdf: {mol2_path}\nobabel error: {msg}')


def prepare_ligand_directory(mol2_path: Path, frcmod_path: Path, out_dir: Path, lig_name: str = 'LIG') -> None:
    repo_root = Path(__file__).resolve().parents[2]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))
    from easybfe.smff.utils import convert_to_xml

    out_dir.mkdir(parents=True, exist_ok=True)

    # 1) SDF for Ligand.from_directory()
    sdf_path = out_dir / f'{lig_name}.sdf'
    convert_lig_mol2_to_sdf(mol2_path, sdf_path)

    # Keep original Amber parameter sources together with generated artifacts.
    shutil.copyfile(mol2_path, out_dir / f'{lig_name}.mol2')
    shutil.copyfile(frcmod_path, out_dir / f'{lig_name}.frcmod')

    # 2) Use tleap to generate AMBER topology from provided mol2+frcmod.
    with tempfile.TemporaryDirectory() as tmp:
        tmpd = Path(tmp)
        leap_in = tmpd / 'leap.in'
        leap_in.write_text(
            '\n'.join([
                'source leaprc.gaff2',
                f'loadamberparams {frcmod_path.resolve()}',
                f'lig = loadmol2 {mol2_path.resolve()}',
                f'set lig name {lig_name}',
                f'set lig.1 name {lig_name}',
                f'saveamberparm lig {lig_name}.prmtop {lig_name}.inpcrd',
                'quit',
            ]) + '\n'
        )

        proc = subprocess.run(['tleap', '-f', 'leap.in'], cwd=tmpd, capture_output=True, text=True)
        if proc.returncode != 0:
            msg = (proc.stderr or proc.stdout or '').strip()
            raise RuntimeError(f'tleap failed when building ligand parameters:\n{msg}')

        prmtop_tmp = tmpd / f'{lig_name}.prmtop'
        inpcrd_tmp = tmpd / f'{lig_name}.inpcrd'

        struct = parmed.load_file(str(prmtop_tmp), xyz=str(inpcrd_tmp))
        for residue in struct.residues:
            residue.name = lig_name

        prmtop_out = out_dir / f'{lig_name}.prmtop'
        inpcrd_out = out_dir / f'{lig_name}.inpcrd'
        pdb_out = out_dir / f'{lig_name}.pdb'
        xml_out = out_dir / f'{lig_name}.xml'

        struct.save(str(prmtop_out), overwrite=True)
        struct.save(str(inpcrd_out), overwrite=True)
        app.PDBFile.writeFile(struct.topology, struct.positions, str(pdb_out), keepIds=True)
        convert_to_xml(struct, xml_out)


def prepare_fixed_protein(std_pdb: Path, fixed_pdb: Path) -> None:
    repo_root = Path(__file__).resolve().parents[2]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))
    from easybfe.protein_prep import ProteinFixer

    fixer = ProteinFixer(str(std_pdb), wizard=False)
    fixer.run(out=str(fixed_pdb))


def main() -> None:
    here = Path(__file__).resolve().parent
    pdb_dir = here / 'pdbs'
    src = pdb_dir / '3pzw_sub.pdb'

    out_full = pdb_dir / '3pzw_sub_full_no_lig.pdb'
    out_std = pdb_dir / '3pzw_sub_std_no_lig.pdb'
    out_lig = pdb_dir / '3pzw_sub_lig_only.pdb'
    out_std_fixed = pdb_dir / '3pzw_sub_std_no_lig_fixed.pdb'

    split_pdb(src, out_full, out_std, out_lig, lig_name='LIG')
    prepare_fixed_protein(out_std, out_std_fixed)

    lig_mol2 = pdb_dir / 'ncaa_lib_14' / 'LIG_AM1BCC-GAFF.mol2'
    lig_frcmod = pdb_dir / 'ncaa_lib_14' / 'LIG_AM1BCC-GAFF.frcmod'

    lig_dir = here / 'ligands' / 'LIG'
    prepare_ligand_directory(lig_mol2, lig_frcmod, lig_dir, lig_name='LIG')

    print(f'Wrote: {out_full}')
    print(f'Wrote: {out_std}')
    print(f'Wrote: {out_lig}')
    print(f'Wrote: {out_std_fixed}')
    print(f'Wrote ligand dir: {lig_dir}')


if __name__ == '__main__':
    main()
