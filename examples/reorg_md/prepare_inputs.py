#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

import openmm.app as app
import parmed
from parmed.amber import AmberParameterSet
from parmed.openmm.parameters import OpenMMParameterSet

STD_PROTEIN = {
    'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','HID','HIE','HIP',
    'ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL',
    'ASH','GLH','LYN','CYM','CYX','ACE','NME'
}

METAL_RESIDUES = ['AN1', 'FE1', 'HD1', 'HD2', 'HD3', 'IE1', 'O11']
METAL_COORDINATION_ATOMS = {
    'HD1': 'NE2',  # Y1
    'HD2': 'NE2',  # Y2
    'HD3': 'NE2',  # Y3
    'AN1': 'OD1',  # Y4
    'IE1': 'OXT',  # Y5
    'O11': 'O',    # Y6
}


def residue_name(line: str) -> str:
    return line[17:20].strip()


def enforce_fe_element_column(pdb_path: Path) -> None:
    """Ensure FE1 residue atoms carry element FE in PDB columns 77-78."""
    lines = pdb_path.read_text().splitlines()
    out: list[str] = []
    for ln in lines:
        rec = ln[:6].strip()
        if rec in {'ATOM', 'HETATM'} and residue_name(ln) == 'FE1':
            padded = ln.ljust(78)
            ln = f"{padded[:76]}FE{padded[78:]}"
        out.append(ln)
    pdb_path.write_text('\n'.join(out) + '\n')


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


def _read_mol2_bonds_by_name(mol2_path: Path) -> list[tuple[str, str]]:
    lines = mol2_path.read_text().splitlines()
    atom_map: dict[int, str] = {}
    bonds: list[tuple[str, str]] = []

    in_atom = False
    in_bond = False
    for ln in lines:
        if ln.startswith('@<TRIPOS>ATOM'):
            in_atom = True
            in_bond = False
            continue
        if ln.startswith('@<TRIPOS>BOND'):
            in_atom = False
            in_bond = True
            continue
        if ln.startswith('@<TRIPOS>'):
            in_atom = False
            in_bond = False
            continue

        if in_atom and ln.strip():
            parts = ln.split()
            if len(parts) >= 2:
                atom_map[int(parts[0])] = parts[1]
        elif in_bond and ln.strip():
            parts = ln.split()
            if len(parts) >= 4:
                a1 = atom_map.get(int(parts[1]))
                a2 = atom_map.get(int(parts[2]))
                if a1 and a2:
                    if a1 <= a2:
                        bonds.append((a1, a2))
                    else:
                        bonds.append((a2, a1))

    return sorted(set(bonds))


def add_nonstd_conect_from_templates(pdb_path: Path, lib_dir: Path, nonstd_resnames: list[str]) -> None:
    lines = pdb_path.read_text().splitlines()

    template_bonds: dict[str, list[tuple[str, str]]] = {}
    for resn in nonstd_resnames:
        mol2 = lib_dir / f'{resn}_AM1BCC-GAFF.mol2'
        template_bonds[resn] = _read_mol2_bonds_by_name(mol2)

    residues = []
    curr = None
    residue_seq = 0

    def flush_curr():
        nonlocal curr
        if curr is not None:
            residues.append(curr)
            curr = None

    for ln in lines:
        rec = ln[:6].strip()
        if rec in {'ATOM', 'HETATM'}:
            serial = int(ln[6:11])
            aname = ln[12:16].strip()
            resn = ln[17:20].strip()
            chain = ln[21:22]
            resid = ln[22:26].strip()
            key = (chain, resid, resn)
            if curr is None or curr['key'] != key:
                flush_curr()
                residue_seq += 1
                curr = {
                    'key': key,
                    'resn': resn,
                    'seq': residue_seq,
                    'atoms': {},
                }
            curr['atoms'][aname] = serial
        elif rec == 'TER':
            flush_curr()
            residue_seq += 1
        elif rec == 'END':
            flush_curr()

    bond_pairs: set[tuple[int, int]] = set()

    # 1) Internal nonstandard-residue bonds from mol2 templates.
    for r in residues:
        resn = r['resn']
        if resn not in template_bonds:
            continue
        amap = r['atoms']
        for a1_name, a2_name in template_bonds[resn]:
            a1 = amap.get(a1_name)
            a2 = amap.get(a2_name)
            if a1 is None or a2 is None:
                continue
            pair = (a1, a2) if a1 < a2 else (a2, a1)
            bond_pairs.add(pair)

    # 2) Backbone peptide links across standard/nonstandard boundaries.
    nonstd_set = set(nonstd_resnames)
    for i in range(len(residues) - 1):
        r1 = residues[i]
        r2 = residues[i + 1]
        if r2['seq'] != r1['seq'] + 1:
            continue

        if (r1['resn'] in nonstd_set) or (r2['resn'] in nonstd_set):
            c_serial = r1['atoms'].get('C')
            n_serial = r2['atoms'].get('N')
            if c_serial is not None and n_serial is not None:
                pair = (c_serial, n_serial) if c_serial < n_serial else (n_serial, c_serial)
                bond_pairs.add(pair)

    # 3) MCPB bonded-model metal coordination links (FE1 to donor atoms).
    by_name: dict[str, list[dict]] = {}
    for r in residues:
        by_name.setdefault(r['resn'], []).append(r)

    fe_residues = by_name.get('FE1', [])
    if len(fe_residues) == 1:
        fe_serial = fe_residues[0]['atoms'].get('FE1')
        if fe_serial is not None:
            for resn, donor_atom in METAL_COORDINATION_ATOMS.items():
                donor_res = by_name.get(resn, [])
                if len(donor_res) != 1:
                    continue
                donor_serial = donor_res[0]['atoms'].get(donor_atom)
                if donor_serial is None:
                    continue
                pair = (fe_serial, donor_serial) if fe_serial < donor_serial else (donor_serial, fe_serial)
                bond_pairs.add(pair)

    out = []
    for ln in lines:
        rec = ln[:6].strip()
        if rec in {'CONECT', 'END'}:
            continue
        out.append(ln)

    for a, b in sorted(bond_pairs):
        out.append(f'CONECT{a:5d}{b:5d}')
        out.append(f'CONECT{b:5d}{a:5d}')

    out.append('END')
    pdb_path.write_text('\n'.join(out) + '\n')


def _run_tleap(input_text: str, workdir: Path, required_outputs: list[str] | None = None) -> None:
    leap_in = workdir / 'leap.in'
    leap_in.write_text(input_text)
    proc = subprocess.run(['tleap', '-f', 'leap.in'], cwd=workdir, capture_output=True, text=True)
    if proc.returncode != 0:
        msg = (proc.stderr or proc.stdout or '').strip()
        raise RuntimeError(f'tleap failed:\n{msg}')
    if required_outputs:
        missing = [x for x in required_outputs if not (workdir / x).is_file()]
        if missing:
            msg = (proc.stderr or proc.stdout or '').strip()
            raise RuntimeError(f'tleap did not produce expected output files {missing}.\n{msg}')


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

    # 2) Build prmtop/inpcrd/xml from provided mol2+frcmod.
    with tempfile.TemporaryDirectory() as tmp:
        tmpd = Path(tmp)
        leap_txt = '\n'.join([
            'source leaprc.gaff',
            f'loadamberparams {frcmod_path.resolve()}',
            f'lig = loadmol2 {mol2_path.resolve()}',
            f'set lig name {lig_name}',
            f'set lig.1 name {lig_name}',
            f'saveamberparm lig {lig_name}.prmtop {lig_name}.inpcrd',
            'quit',
        ]) + '\n'
        _run_tleap(leap_txt, tmpd, required_outputs=[f'{lig_name}.prmtop', f'{lig_name}.inpcrd'])

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


def _guess_element_from_atom_name(atom_name: str, fallback: str) -> str:
    name = atom_name.strip().upper()
    if name.startswith('FE'):
        return 'Fe'
    if name.startswith('CL'):
        return 'Cl'
    if name.startswith('BR'):
        return 'Br'
    if name.startswith('ZN'):
        return 'Zn'
    lead = name[:1]
    mapping = {
        'C': 'C',
        'N': 'N',
        'O': 'O',
        'H': 'H',
        'S': 'S',
        'P': 'P',
        'F': 'F',
        'I': 'I',
    }
    return mapping.get(lead, fallback)


def normalize_type_elements(xml_path: Path) -> None:
    tree = ET.parse(str(xml_path))
    root = tree.getroot()

    type_to_atom: dict[str, str] = {}
    residues = root.find('Residues')
    if residues is not None:
        for residue in residues.findall('Residue'):
            for atom in residue.findall('Atom'):
                aname = atom.get('name')
                tname = atom.get('type')
                if aname and tname:
                    type_to_atom[tname] = aname

    atom_types = root.find('AtomTypes')
    if atom_types is not None:
        for typ in atom_types.findall('Type'):
            tname = typ.get('name')
            elem = typ.get('element', 'C')
            if tname in type_to_atom:
                typ.set('element', _guess_element_from_atom_name(type_to_atom[tname], elem))

    tree.write(str(xml_path))


def add_peptide_external_bonds(xml_path: Path) -> None:
    tree = ET.parse(str(xml_path))
    root = tree.getroot()
    residues = root.find('Residues')
    if residues is None:
        tree.write(str(xml_path))
        return

    for residue in residues.findall('Residue'):
        atom_names = {atom.get('name') for atom in residue.findall('Atom')}
        existing = {eb.get('atomName') for eb in residue.findall('ExternalBond')}
        if ('N' in atom_names) and ('C' in atom_names):
            if 'N' not in existing:
                ET.SubElement(residue, 'ExternalBond', {'atomName': 'N'})
            # If residue has OXT, treat as C-terminus and avoid external C bond.
            if ('OXT' not in atom_names) and ('C' not in existing):
                ET.SubElement(residue, 'ExternalBond', {'atomName': 'C'})

    tree.write(str(xml_path))


def add_metal_external_bonds(xml_path: Path) -> None:
    """Add MCPB coordination ExternalBond tags for metal center residues."""
    coord_external = {
        'HD1': ['NE2'],
        'HD2': ['NE2'],
        'HD3': ['NE2'],
        'AN1': ['OD1'],
        'IE1': ['OXT'],
        'O11': ['O'],
        # FE1 coordinates six donor atoms in this site.
        'FE1': ['FE1', 'FE1', 'FE1', 'FE1', 'FE1', 'FE1'],
    }

    tree = ET.parse(str(xml_path))
    root = tree.getroot()
    residues = root.find('Residues')
    if residues is None:
        tree.write(str(xml_path))
        return

    for residue in residues.findall('Residue'):
        resn = residue.get('name')
        if resn not in coord_external:
            continue
        existing = [eb.get('atomName') for eb in residue.findall('ExternalBond')]
        for atom_name in coord_external[resn]:
            if existing.count(atom_name) < coord_external[resn].count(atom_name):
                ET.SubElement(residue, 'ExternalBond', {'atomName': atom_name})
                existing.append(atom_name)

    tree.write(str(xml_path))


def write_mcpb_parameter_ffxml(frcmod_path: Path, out_xml: Path) -> Path:
    """Convert MCPB frcmod (bonded terms) into standalone OpenMM ffxml."""
    params = AmberParameterSet(str(frcmod_path))
    omm_params = OpenMMParameterSet.from_parameterset(params)
    omm_params.write(str(out_xml))
    return out_xml


def prepare_metal_center_ffxml(lib_dir: Path, out_dir: Path, residue_names: list[str]) -> list[Path]:
    repo_root = Path(__file__).resolve().parents[2]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))
    from easybfe.smff.utils import convert_to_xml

    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    for resn in residue_names:
        mol2 = lib_dir / f'{resn}_AM1BCC-GAFF.mol2'
        frcmod = lib_dir / f'{resn}_AM1BCC-GAFF.frcmod'
        if not mol2.is_file() or not frcmod.is_file():
            raise FileNotFoundError(f'Missing files for residue {resn}: {mol2} / {frcmod}')

        with tempfile.TemporaryDirectory() as tmp:
            tmpd = Path(tmp)
            leap_lines = [
                'source leaprc.protein.ff14SB',
                'source leaprc.gaff',
            ]
            # MCPB-style residues can cross-reference custom atom types across frcmods.
            for load_resn in residue_names:
                leap_lines.append(f'loadamberparams {(lib_dir / f"{load_resn}_AM1BCC-GAFF.frcmod").resolve()}')
            leap_lines.extend([
                f'res = loadmol2 {mol2.resolve()}',
                f'set res name {resn}',
                f'set res.1 name {resn}',
                f'saveamberparm res {resn}.prmtop {resn}.inpcrd',
                'quit',
            ])
            _run_tleap(
                '\n'.join(leap_lines) + '\n',
                tmpd,
                required_outputs=[f'{resn}.prmtop', f'{resn}.inpcrd'],
            )

            struct = parmed.load_file(str(tmpd / f'{resn}.prmtop'), xyz=str(tmpd / f'{resn}.inpcrd'))
            for residue in struct.residues:
                residue.name = resn

            xml_path = out_dir / f'{resn}.xml'
            convert_to_xml(struct, xml_path)
            normalize_type_elements(xml_path)
            add_peptide_external_bonds(xml_path)
            add_metal_external_bonds(xml_path)
            written.append(xml_path)

    return written


def main() -> None:
    here = Path(__file__).resolve().parent
    pdb_dir = here / 'pdbs'
    src = pdb_dir / '3pzw_sub.pdb'

    out_full = pdb_dir / '3pzw_sub_full_no_lig.pdb'
    out_std = pdb_dir / '3pzw_sub_std_no_lig.pdb'
    out_lig = pdb_dir / '3pzw_sub_lig_only.pdb'

    split_pdb(src, out_full, out_std, out_lig, lig_name='LIG')
    enforce_fe_element_column(out_full)
    add_nonstd_conect_from_templates(out_full, pdb_dir / "ncaa_lib_14", METAL_RESIDUES)

    lig_mol2 = pdb_dir / 'ncaa_lib_14' / 'LIG_AM1BCC-GAFF.mol2'
    lig_frcmod = pdb_dir / 'ncaa_lib_14' / 'LIG_AM1BCC-GAFF.frcmod'

    lig_dir = here / 'ligands' / 'LIG'
    prepare_ligand_directory(lig_mol2, lig_frcmod, lig_dir, lig_name='LIG')

    ffxml_dir = here / 'ffxml'
    metal_xmls = prepare_metal_center_ffxml(pdb_dir / 'ncaa_lib_14', ffxml_dir, METAL_RESIDUES)
    mcpb_param_xml = write_mcpb_parameter_ffxml(
        pdb_dir / 'ncaa_lib_14' / 'O11_AM1BCC-GAFF.frcmod',
        ffxml_dir / 'mcpb_metal_params.xml',
    )

    print(f'Wrote: {out_full}')
    print(f'Wrote: {out_std}')
    print(f'Wrote: {out_lig}')
    print(f'Wrote ligand dir: {lig_dir}')
    print('Wrote metal ffxml files:')
    for x in metal_xmls:
        print(f'  - {x}')
    print(f'Wrote metal parameter ffxml: {mcpb_param_xml}')


if __name__ == '__main__':
    main()
