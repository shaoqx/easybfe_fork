#!/usr/bin/env python3
from __future__ import annotations

import copy
import json
import os
import tempfile
from collections import Counter
from pathlib import Path

import parmed
import openmm as mm
import openmm.app as app
import openmm.unit as unit

import sys
REPO = Path('/home/shaoq1/bin/easybfe').resolve()
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from easybfe.config import read_file
from easybfe.config.amber.abfe import AmberAbfeConfig
from easybfe.core.ligand import Ligand
from easybfe.core.protein import Protein
from easybfe.amber.prep_ligand_abfe import (
    _build_residue_templates,
    _collect_residue_names_from_ffxml,
    _patch_missing_bond_types,
)
from easybfe.amber.prep_utils import computeBoxVectorsWithPadding, shiftToBoxCenter


def atom_label(atom) -> str:
    return f"{atom.residue.name}{atom.residue.idx+1}:{atom.name}({getattr(atom,'type',None)})"


def omm_top_atom_label(atom) -> str:
    return f"{atom.residue.name}{int(atom.residue.id)}:{atom.name}"


def is_custom_residue_name(name: str, custom_residue_names: set[str]) -> bool:
    if name in custom_residue_names:
        return True
    if name in {"HID", "HIE", "HIP"}:
        return True
    return False


def interaction_kind_for_bond(a1, a2) -> str:
    names = {a1.name, a2.name}
    resnames = {a1.residue.name, a2.residue.name}
    if names == {"C", "N"}:
        return "peptide_boundary"
    if "FE1" in resnames:
        return "metal_coordination"
    return "other_inter_residue"


def summarize_missing(parmed_struct, custom_residue_names: set[str]):
    summary = {}

    missing_bonds = [b for b in parmed_struct.bonds if b.type is None]
    missing_angles = [a for a in parmed_struct.angles if a.type is None]
    missing_dihedrals = [d for d in parmed_struct.dihedrals if d.type is None]

    def bond_is_boundary(b):
        return (b.atom1.residue is not b.atom2.residue) and (
            is_custom_residue_name(b.atom1.residue.name, custom_residue_names)
            or is_custom_residue_name(b.atom2.residue.name, custom_residue_names)
        )

    def angle_is_boundary(a):
        residues = {a.atom1.residue, a.atom2.residue, a.atom3.residue}
        return len(residues) > 1 and any(
            is_custom_residue_name(x.name, custom_residue_names) for x in residues
        )

    def dihedral_is_boundary(d):
        residues = {d.atom1.residue, d.atom2.residue, d.atom3.residue, d.atom4.residue}
        return len(residues) > 1 and any(
            is_custom_residue_name(x.name, custom_residue_names) for x in residues
        )

    boundary_bonds = [b for b in missing_bonds if bond_is_boundary(b)]
    boundary_angles = [a for a in missing_angles if angle_is_boundary(a)]
    boundary_dihedrals = [d for d in missing_dihedrals if dihedral_is_boundary(d)]

    summary["missing_bonds_total"] = len(missing_bonds)
    summary["missing_angles_total"] = len(missing_angles)
    summary["missing_dihedrals_total"] = len(missing_dihedrals)
    summary["missing_boundary_bonds"] = len(boundary_bonds)
    summary["missing_boundary_angles"] = len(boundary_angles)
    summary["missing_boundary_dihedrals"] = len(boundary_dihedrals)

    detail = {
        "boundary_bonds": [],
        "boundary_angles": [],
        "boundary_dihedrals": [],
    }

    for b in boundary_bonds:
        detail["boundary_bonds"].append({
            "kind": interaction_kind_for_bond(b.atom1, b.atom2),
            "a1": atom_label(b.atom1),
            "a2": atom_label(b.atom2),
        })

    for a in boundary_angles:
        detail["boundary_angles"].append({
            "atoms": [atom_label(a.atom1), atom_label(a.atom2), atom_label(a.atom3)],
            "residues": [a.atom1.residue.name, a.atom2.residue.name, a.atom3.residue.name],
        })

    for d in boundary_dihedrals:
        detail["boundary_dihedrals"].append({
            "atoms": [atom_label(d.atom1), atom_label(d.atom2), atom_label(d.atom3), atom_label(d.atom4)],
            "residues": [d.atom1.residue.name, d.atom2.residue.name, d.atom3.residue.name, d.atom4.residue.name],
            "improper": bool(d.improper),
        })

    return summary, detail


def inspect_openmm_forces(modeller: app.Modeller, system: mm.System, custom_residue_names: set[str]):
    top_atoms = list(modeller.topology.atoms())
    top_bonds = list(modeller.topology.bonds())

    bond_force = None
    angle_force = None
    torsion_force = None
    for force in system.getForces():
        if isinstance(force, mm.HarmonicBondForce):
            bond_force = force
        elif isinstance(force, mm.HarmonicAngleForce):
            angle_force = force
        elif isinstance(force, mm.PeriodicTorsionForce):
            torsion_force = force

    bond_pairs = Counter()
    if bond_force is not None:
        for idx in range(bond_force.getNumBonds()):
            a, b, _, _ = bond_force.getBondParameters(idx)
            bond_pairs[tuple(sorted((int(a), int(b))))] += 1

    angle_triples = Counter()
    if angle_force is not None:
        for idx in range(angle_force.getNumAngles()):
            a, b, c, _, _ = angle_force.getAngleParameters(idx)
            angle_triples[(int(a), int(b), int(c))] += 1

    torsion_quads = Counter()
    if torsion_force is not None:
        for idx in range(torsion_force.getNumTorsions()):
            a, b, c, d, _, _, _ = torsion_force.getTorsionParameters(idx)
            torsion_quads[(int(a), int(b), int(c), int(d))] += 1

    boundary_report = []
    for bond in top_bonds:
        a1, a2 = bond.atom1, bond.atom2
        if a1.residue == a2.residue:
            continue
        if not (
            is_custom_residue_name(a1.residue.name, custom_residue_names)
            or is_custom_residue_name(a2.residue.name, custom_residue_names)
        ):
            continue
        i = a1.index
        j = a2.index
        entry = {
            "kind": interaction_kind_for_bond(a1, a2),
            "bond": [omm_top_atom_label(a1), omm_top_atom_label(a2)],
            "bond_term_in_system": bond_pairs.get(tuple(sorted((i, j))), 0),
            "angle_terms_involving_boundary": 0,
            "torsion_terms_involving_boundary": 0,
        }
        for (x, y, z), count in angle_triples.items():
            if tuple((x, y)) == (i, j) or tuple((y, z)) == (i, j) or tuple((x, y)) == (j, i) or tuple((y, z)) == (j, i):
                entry["angle_terms_involving_boundary"] += count
        for (x, y, z, w), count in torsion_quads.items():
            pairs = [(x, y), (y, z), (z, w)]
            if (i, j) in pairs or (j, i) in pairs:
                entry["torsion_terms_involving_boundary"] += count
        boundary_report.append(entry)
    return boundary_report


def main():
    reorg_dir = Path('/home/shaoq1/bin/easybfe/examples/reorg_md').resolve()
    out_dir = Path('/tmp/easybfe_openmm_debug')
    out_dir.mkdir(parents=True, exist_ok=True)

    os.chdir(reorg_dir)
    cfg_dict = read_file(reorg_dir / 'config_reorg_5ns.yaml')
    config = AmberAbfeConfig.model_validate(cfg_dict)
    leg_cfg = config.complex

    ligand = Ligand.from_directory(reorg_dir / 'ligands' / 'LIG')
    protein = Protein.from_pdb(reorg_dir / 'pdbs' / '3pzw_sub_full_no_lig.pdb', name='3pzw_sub_full_no_lig')

    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        ligand.dump(tmpdir)
        ligand_pdb = app.PDBFile(str(tmpdir / 'LIG.pdb'))
        ff = app.ForceField(*leg_cfg.forcefields, str(tmpdir / 'LIG.xml'))
        ligand_resname = list(ligand_pdb.topology.residues())[0].name
        custom_residue_names = _collect_residue_names_from_ffxml(leg_cfg.extra_ff)

        modeller = app.Modeller(app.Topology(), [])
        modeller.add(ligand_pdb.topology, ligand_pdb.positions)
        protein_openmm = protein.to_openmm()
        modeller.add(protein_openmm.topology, protein_openmm.positions)

        buffer = leg_cfg.buffer / 10 * unit.nanometers
        box_vectors = computeBoxVectorsWithPadding(modeller.positions, buffer, leg_cfg.box_shape)
        modeller.positions = shiftToBoxCenter(modeller.positions, box_vectors)
        modeller.topology.setPeriodicBoxVectors(box_vectors)
        residue_templates = _build_residue_templates(modeller.topology, custom_residue_names, ligand_resname)
        modeller.addSolvent(
            forcefield=ff,
            model=leg_cfg.water_model,
            neutralize=True,
            ionicStrength=leg_cfg.ionic_strength * unit.molar,
            residueTemplates=residue_templates,
        )

        residue_templates = _build_residue_templates(modeller.topology, custom_residue_names, ligand_resname)
        system = ff.createSystem(
            modeller.topology,
            nonbondedMethod=app.PME,
            constraints=None,
            rigidWater=False,
            residueTemplates=residue_templates,
        )

        parmed_struct_pre = parmed.openmm.load_topology(modeller.topology, system, xyz=modeller.positions)
        summary_pre, detail_pre = summarize_missing(parmed_struct_pre, custom_residue_names)

        parmed_struct_post = copy.deepcopy(parmed_struct_pre)
        patch_error = None
        try:
            _patch_missing_bond_types(parmed_struct_post)
        except Exception as exc:
            patch_error = repr(exc)
        summary_post, detail_post = summarize_missing(parmed_struct_post, custom_residue_names)

        omm_boundary = inspect_openmm_forces(modeller, system, custom_residue_names)

        report = {
            "custom_residue_names": sorted(custom_residue_names),
            "forcefield_files": leg_cfg.forcefields + [str(tmpdir / 'LIG.xml')],
            "summary_pre_patch": summary_pre,
            "summary_post_patch": summary_post,
            "patch_error": patch_error,
            "openmm_boundary_terms": omm_boundary,
            "missing_detail_pre_patch": detail_pre,
            "missing_detail_post_patch": detail_post,
        }

        (out_dir / 'report.json').write_text(json.dumps(report, indent=2))

        lines = []
        lines.append('=== easybfe OpenMM boundary debug ===')
        lines.append(f"custom residues: {sorted(custom_residue_names)}")
        lines.append('')
        lines.append('--- Missing terms before _patch_missing_bond_types ---')
        for k, v in summary_pre.items():
            lines.append(f"{k}: {v}")
        lines.append('')
        lines.append('--- Missing terms after _patch_missing_bond_types ---')
        for k, v in summary_post.items():
            lines.append(f"{k}: {v}")
        if patch_error:
            lines.append(f"patch_error: {patch_error}")
        lines.append('')
        lines.append('--- Boundary bond terms directly present in OpenMM system ---')
        for item in omm_boundary[:40]:
            lines.append(json.dumps(item, ensure_ascii=False))
        lines.append('')
        lines.append('--- Example missing boundary angles (pre-patch) ---')
        for item in detail_pre['boundary_angles'][:40]:
            lines.append(json.dumps(item, ensure_ascii=False))
        lines.append('')
        lines.append('--- Example missing boundary dihedrals (pre-patch) ---')
        for item in detail_pre['boundary_dihedrals'][:60]:
            lines.append(json.dumps(item, ensure_ascii=False))
        (out_dir / 'report.txt').write_text('\n'.join(lines) + '\n')

        try:
            xml_text = mm.XmlSerializer.serialize(system)
            (out_dir / 'system_openmm.xml').write_text(xml_text)
        except Exception as exc:
            (out_dir / 'system_openmm_xml_error.txt').write_text(repr(exc) + '\n')

    print(out_dir / 'report.txt')
    print(out_dir / 'report.json')
    if (out_dir / 'system_openmm.xml').exists():
        print(out_dir / 'system_openmm.xml')


if __name__ == '__main__':
    main()
