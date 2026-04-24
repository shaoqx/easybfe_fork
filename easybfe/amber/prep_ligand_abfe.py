import os
from pathlib import Path
import logging
from typing import Optional
import xml.etree.ElementTree as ET

import numpy as np
import openmm.app as app
import openmm.unit as unit
import parmed

from .prep_utils import *
from ..config import AmberFepSimulationConfig, AmberWtSettings
from ..config.amber.abfe import AmberAbfeConfig, BoreschRestraintGeneratorConfig
from .workflow import Step, Workflow, create_script_for_workflows
from ..core import Ligand, Protein
from .boresch import BORESCH_FINDER_REGISTRY, BoreschRestraint, compute_boresch_energy
from ..parallel import run_func_parallel


logger = logging.getLogger(__name__)


def _collect_residue_names_from_ffxml(ff_paths: list[str]) -> set[str]:
    """Collect residue template names from local ffxml files."""
    residue_names: set[str] = set()
    for ff_path in ff_paths:
        xml_path = Path(ff_path).expanduser()
        if not xml_path.is_file():
            continue
        try:
            root = ET.parse(str(xml_path)).getroot()
        except Exception:
            continue
        residues = root.find('Residues')
        if residues is None:
            continue
        for residue in residues.findall('Residue'):
            name = residue.get('name')
            if name:
                residue_names.add(name)
    return residue_names


def _build_residue_templates(
    topology: app.Topology,
    custom_residue_names: set[str],
    ligand_resname: str | None = None,
) -> dict[app.topology.Residue, str]:
    """Build explicit OpenMM residueTemplates mapping."""
    residue_templates: dict[app.topology.Residue, str] = {}
    has_hd_histidine_templates = any(name in custom_residue_names for name in ('HD1', 'HD2', 'HD3'))
    has_asn_like_templates = any(name.startswith('AN') for name in custom_residue_names)
    has_ile_like_templates = any(name.startswith('IE') for name in custom_residue_names)

    def _terminal_variant(base: str, residue: app.topology.Residue) -> str:
        atom_names = {atom.name for atom in residue.atoms()}
        if 'OXT' in atom_names:
            return f'C{base}'
        if ('H1' in atom_names) or ('H2' in atom_names) or ('H3' in atom_names):
            return f'N{base}'
        return base

    def _his_base(residue: app.topology.Residue) -> str:
        atom_names = {atom.name for atom in residue.atoms()}
        has_hd1 = 'HD1' in atom_names
        has_he2 = 'HE2' in atom_names
        if has_hd1 and has_he2:
            return 'HIP'
        if has_he2:
            return 'HIE'
        return 'HID'

    for residue in topology.residues():
        if (ligand_resname is not None) and (residue.name == ligand_resname):
            residue_templates[residue] = residue.name
        elif residue.name in custom_residue_names:
            residue_templates[residue] = residue.name
        elif has_hd_histidine_templates and residue.name == 'HIS':
            residue_templates[residue] = _terminal_variant(_his_base(residue), residue)
        elif has_asn_like_templates and residue.name == 'ASN':
            residue_templates[residue] = _terminal_variant('ASN', residue)
        elif has_ile_like_templates and residue.name == 'ILE':
            residue_templates[residue] = _terminal_variant('ILE', residue)
    return residue_templates


def _patch_missing_bond_types(parmed_struct: parmed.Structure) -> None:
    """Fill missing bond types caused by custom residue boundary bonds."""
    missing_bonds = [bond for bond in parmed_struct.bonds if bond.type is None]
    if not missing_bonds:
        return

    type_lookup: dict[tuple[str, str], parmed.BondType] = {}
    for bond in parmed_struct.bonds:
        if bond.type is None:
            continue
        key = tuple(sorted((bond.atom1.type, bond.atom2.type)))
        type_lookup.setdefault(key, bond.type)

    # Peptide-like C-N reference present in this system.
    peptide_ref = type_lookup.get(tuple(sorted(('C4', 'N1'))))
    fallback_pairs = {tuple(sorted(('C4', 'N2'))), tuple(sorted(('C6', 'N1')))}

    for bond in missing_bonds:
        key = tuple(sorted((bond.atom1.type, bond.atom2.type)))
        bond_type = type_lookup.get(key)
        if (bond_type is None) and (key in fallback_pairs):
            bond_type = peptide_ref
        if bond_type is None:
            a1, a2 = bond.atom1, bond.atom2
            raise RuntimeError(
                f"Missing bonded parameter for {a1.residue.name}:{a1.name}({a1.type}) - "
                f"{a2.residue.name}:{a2.name}({a2.type})"
            )
        bond.type = bond_type

    logger.warning(
        "Patched %d missing bond type(s) at custom residue boundaries.",
        len(missing_bonds),
    )


def setup_ligand_abfe_leg(
    ligand: Ligand, 
    protein: Protein | None,
    config: AmberFepSimulationConfig,
    wdir: os.PathLike,
    duplicate_ligand: bool = False,
    restraints: BoreschRestraint | None = None,
    basename: str | None = None,
):  
    # setup workding dir
    wdir = Path(wdir).expanduser().resolve()
    wdir.mkdir(exist_ok=True)
    basename = wdir.stem if not basename else basename

    ligand.dump(wdir)
    ligand_pdb = app.PDBFile(str(wdir / f'{ligand.name}.pdb'))

    # charges
    ligand_charge = compute_net_charge_from_openmm_system(
        app.ForceField(str(wdir / f'{ligand.name}.xml')).createSystem(ligand_pdb.topology)
    )

    # force field initialization
    ff = app.ForceField(*config.forcefields, str(wdir / f'{ligand.name}.xml'))
    ligand_resname = list(ligand_pdb.topology.residues())[0].name
    custom_residue_names = _collect_residue_names_from_ffxml(config.extra_ff)

    # setup systems
    modeller = app.Modeller(app.Topology(), [])
    modeller.add(ligand_pdb.topology, ligand_pdb.positions)
    
    # use for resolve restraints in protein-ligand complex
    if duplicate_ligand:
        modeller.add(ligand_pdb.topology, ligand_pdb.positions)

    if protein:
        protein_openmm = protein.to_openmm()
        modeller.add(protein_openmm.topology, protein_openmm.positions)
    
    buffer = config.buffer / 10 * unit.nanometers
    box_vectors = computeBoxVectorsWithPadding(modeller.positions, buffer, config.box_shape)
    modeller.positions = shiftToBoxCenter(modeller.positions, box_vectors)
    modeller.topology.setPeriodicBoxVectors(box_vectors)
    assert not config.gas_phase, 'Gas-phase ABFE are ill-defined!'
    residue_templates = _build_residue_templates(modeller.topology, custom_residue_names, ligand_resname)
    modeller.addSolvent(
        forcefield=ff,
        model=config.water_model,
        neutralize=True,
        ionicStrength=config.ionic_strength * unit.molar,
        residueTemplates=residue_templates
    )

    # generate masks
    mode = 'abfe_restr' if duplicate_ligand else 'abfe'
    num_ligand_atoms = len(list(ligand_pdb.topology.atoms()))
    mask = generate_amber_mask(num_ligand_atoms, -1, {}, mode=mode)

    # alchemical water
    alchem_waters, rst_settings = [], []
    if ligand_charge != 0:
        logger.info(f"ABFE with for ligand {ligand.name} with net charge {int(ligand_charge)}")
        if duplicate_ligand:
            fix_excess_charge(modeller, ligand_charge)
        elif config.use_charge_change and (not config.gas_phase):
            scIndices = list[int](range(ligand_pdb.topology.getNumAtoms()))
            coion_info = create_alchemical_ions(modeller, ligand_charge, 0, scIndices, method=config.charge_change_method)
            if config.add_restraint_for_alchem_water:
                rst_settings = set_alchemical_water_restraints(modeller, scIndices, coion_info)
            alchem_waters = [] if config.use_settle_for_alchemical_water else coion_info['alchemical_water_residues']
            mask = generate_amber_mask(num_ligand_atoms, -1, {}, coion_info, mode=mode)
        else:
            logger.warning("Charge change not enabled. Results are not trustworthy.")
    
    # apply boresch restraints
    if restraints:
        assert protein is not None, 'Boresch restraints must be used with protein'
        ligand_pos_angstrom = np.array([[p.x, p.y, p.z] for p in ligand_pdb.positions]) * 10
        protein_pos_angstrom = np.array([[p.x, p.y, p.z] for p in protein_openmm.positions]) * 10
        restraints.compute_rst_vals(protein_pos_angstrom, ligand_pos_angstrom)
        rst_settings += restraints.make_rst(
            offset=num_ligand_atoms if not duplicate_ligand else num_ligand_atoms*2
        )

    residue_templates = _build_residue_templates(modeller.topology, custom_residue_names, ligand_resname)
    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        constraints=None,
        rigidWater=False,
        residueTemplates=residue_templates,
    )
    parmed_struct = parmed.openmm.load_topology(modeller.topology, system, xyz=modeller.positions)
    _patch_missing_bond_types(parmed_struct)
    
    # Handle Amber special SETTLE water 
    sanitize_water(parmed_struct, 'ALW', alchem_waters)

    # HMR
    if config.do_hmr:
        hydrogen_mass_repartition(parmed_struct, config.hydrogen_mass, config.do_hmr_water)
    
    # output
    parmed_struct.save(str(wdir / f'{basename}.inpcrd'), overwrite=True)
    parmed_struct.save(str(wdir / f'{basename}.prmtop'), overwrite=True)
    parmed_struct.save(str(wdir / f'{basename}.pdb'), overwrite=True)

    # setup workflow
    workflows = []
    for n, clambda in enumerate(config.lambdas):
        steps = []
        for step_template in config.workflow:
            # Build a per-lambda AmberStepConfig with updated cntrl / wt / rst
            update = mask.copy()
            update.update({
                "clambda": clambda, 
                'ntwx': 10*step_template.cntrl.ntwx if (config.reduce_storage and n > 0 and n < (len(config.lambdas)-1)) else step_template.cntrl.ntwx
            })
            base = step_template.cntrl.model_dump()
            base.update(update)
            step_lambda_cntrl = step_template.cntrl.__class__.model_validate(base)
            rst = step_template.rst + rst_settings
            wt = step_template.wt + [AmberWtSettings(type="DUMPFREQ", istep1=step_lambda_cntrl.ofreq)]

            step_config = step_template.model_copy()
            step_config.cntrl = step_lambda_cntrl
            step_config.rst = rst
            step_config.wt = wt

            step = Step(config=step_config)
            steps.append(step)

        lambda_dir = wdir / f'lambda{n}'
        prmtop = wdir / f'{basename}.prmtop'
        inpcrd = wdir / f'{basename}.inpcrd'
        wf = Workflow(wdir=lambda_dir, prmtop=prmtop, inpcrd=inpcrd, steps=steps)
        wf.create()
        workflows.append(wf)
    
    # Use groupfile and MPI to run steps except energy minimization
    create_script_for_workflows(workflows, wdir, config.num_procs)

    return True


def setup_ligand_abfe(
    ligand: Ligand,
    protein: Protein,
    leg_configs: dict[str, AmberFepSimulationConfig],
    output_dir: os.PathLike,
    restraints: BoreschRestraint | BoreschRestraintGeneratorConfig | None = None,
    auto_find_boresch: bool = True,
):
    """Set up ABFE/reorg legs for a single ligand.

    Parameters
    ----------
    leg_configs
        Mapping of leg name -> config. Supported legs are ``solvent``, ``complex``,
        and ``restraint``. Any subset is allowed.
    restraints
        Boresch restraints object/config. If ``None`` and ``auto_find_boresch`` is
        True, a default Boresch restraint will be generated when needed.
    auto_find_boresch
        Backward-compatible switch: keep legacy behavior for direct function calls.
        Config-driven setup passes ``False`` so ``boresch: null`` truly disables it.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    unknown_legs = set(leg_configs) - {'solvent', 'complex', 'restraint'}
    if unknown_legs:
        raise ValueError(f"Unsupported ABFE leg(s): {sorted(unknown_legs)}")
    if len(leg_configs) == 0:
        raise ValueError('No ABFE legs configured')

    needs_boresch = ('complex' in leg_configs) or ('restraint' in leg_configs)
    resolved_restraints: BoreschRestraint | None = None

    if needs_boresch:
        if isinstance(restraints, BoreschRestraint):
            resolved_restraints = restraints
        elif isinstance(restraints, BoreschRestraintGeneratorConfig):
            boresch_config = restraints
            finder_kwargs = {
                'protein': protein,
                'ligand': ligand,
                'wts': tuple(boresch_config.rst_wts),
                **boresch_config.options,
            }
            finder = BORESCH_FINDER_REGISTRY.create(boresch_config.algorithm, **finder_kwargs)
            resolved_restraints = finder.find()
        elif restraints is None and auto_find_boresch:
            boresch_config = BoreschRestraintGeneratorConfig()
            finder_kwargs = {
                'protein': protein,
                'ligand': ligand,
                'wts': tuple(boresch_config.rst_wts),
                **boresch_config.options,
            }
            finder = BORESCH_FINDER_REGISTRY.create(boresch_config.algorithm, **finder_kwargs)
            resolved_restraints = finder.find()

    if 'restraint' in leg_configs and resolved_restraints is None:
        raise ValueError('restraint leg requires Boresch restraints, but none were provided')

    if 'solvent' in leg_configs:
        setup_ligand_abfe_leg(
            ligand,
            None,
            leg_configs['solvent'],
            output_dir / 'solvent',
            basename='system'
        )

    if 'complex' in leg_configs:
        setup_ligand_abfe_leg(
            ligand,
            protein,
            leg_configs['complex'],
            output_dir / 'complex',
            restraints=resolved_restraints,
            basename='system'
        )

    if 'restraint' in leg_configs:
        setup_ligand_abfe_leg(
            ligand,
            protein,
            leg_configs['restraint'],
            output_dir / 'restraint',
            duplicate_ligand=True,
            restraints=resolved_restraints,
            basename='system'
        )

    # Only write standard-state correction when both complex/restraint legs exist.
    if ('complex' in leg_configs) and ('restraint' in leg_configs) and (resolved_restraints is not None):
        boresch_fe = compute_boresch_energy(resolved_restraints.rst_vals, resolved_restraints.rst_wts)
        (output_dir / 'boresch.dat').write_text(str(boresch_fe))


def _resolve_ligand_directory(ligand_base: Path | None, component: os.PathLike) -> Path:
    """Resolve a ligand directory under optional ``ligand_base`` or as an absolute path."""
    comp = Path(component)
    if ligand_base is not None:
        return (Path(ligand_base).expanduser().resolve() / comp).resolve()
    return comp.expanduser().resolve()


def _setup_ligand_abfe_one(
    ligand_path: Path,
    protein: Optional[Protein],
    leg_configs: dict[str, AmberFepSimulationConfig],
    boresch_config: BoreschRestraintGeneratorConfig | None,
    output_dir: Path,
) -> None:
    """Load ligand from path and call setup_ligand_abfe (used for batch runs)."""
    ligand = Ligand.from_directory(ligand_path)
    setup_ligand_abfe(
        ligand=ligand,
        protein=protein,
        leg_configs=leg_configs,
        restraints=boresch_config,
        auto_find_boresch=False,
        output_dir=output_dir,
    )


def setup_ligand_abfe_from_config(
    config: AmberAbfeConfig,
    num_procs: Optional[int] = None,
) -> None:
    """Run setup_ligand_abfe from an :class:`AmberAbfeConfig`.

    **Ligand directories**

    If :attr:`~AmberAbfeConfig.ligand_base` is set, ligand paths are resolved as
    ``ligand_base / relative_path``; otherwise they are treated as full paths.

    **Output**

    Batch mode requires :attr:`~AmberAbfeConfig.output_base`; each run writes under
    ``output_base / ligand_stem``. Single-ligand mode uses
    ``output_base / ligand.name`` when ``output_base`` is set, otherwise
    :attr:`~AmberAbfeConfig.output_dir` (required in that case).
    """
    assert config.protein is not None, "AmberAbfeConfig.protein must be set"

    leg_configs = {leg: getattr(config, leg) for leg in config.active_legs}
    protein = Protein.from_pdb(config.protein, name=config.protein.stem)
    lig_base = Path(config.ligand_base).expanduser().resolve() if config.ligand_base is not None else None

    if config.ligand_batch is not None and config.ligand is not None:
        raise ValueError(
            "AmberAbfeConfig must set either ligand or ligand_batch, not both"
        )

    use_batch = config.ligand_batch is not None

    if use_batch:
        if config.output_base is None:
            raise ValueError(
                "AmberAbfeConfig.output_base is required for batch mode (non-empty ligand_batch)"
            )
        out_base = Path(config.output_base).expanduser().resolve()
        out_base.mkdir(parents=True, exist_ok=True)
        nprocs = num_procs if num_procs is not None else -1
        args_list = [
            (
                _resolve_ligand_directory(lig_base, path),
                protein,
                leg_configs,
                config.boresch,
                out_base / Path(path).stem,
            )
            for path in config.ligand_batch
        ]
        run_func_parallel(
            _setup_ligand_abfe_one,
            args_list,
            nprocs=nprocs,
            unpack_args=True,
            desc="setup_ligand_abfe",
        )
        return

    if config.ligand is None:
        raise ValueError("AmberAbfeConfig must set either ligand or ligand_batch")

    ligand_dir = _resolve_ligand_directory(lig_base, config.ligand)

    if config.output_base is not None:
        run_out = (
            Path(config.output_base).expanduser().resolve()
            / Path(config.ligand).name
        )
    elif config.output_dir is not None:
        run_out = Path(config.output_dir).expanduser().resolve()
    else:
        raise ValueError(
            "Set output_base or output_dir for single-ligand ABFE setup "
            "(output_dir is required when output_base is not set)"
        )

    run_out.mkdir(parents=True, exist_ok=True)
    ligand = Ligand.from_directory(ligand_dir)
    setup_ligand_abfe(
        ligand=ligand,
        protein=protein,
        leg_configs=leg_configs,
        restraints=config.boresch,
        auto_find_boresch=False,
        output_dir=run_out,
    )
