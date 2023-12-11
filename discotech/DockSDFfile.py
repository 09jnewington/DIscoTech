""" Runs docking on a single SDF file using autodock Vina. """
import yaml
from collections import defaultdict
from enum import Enum
from typing import List
from pathlib import Path

from openbabel import openbabel
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem
from vina import Vina

from discotech.utils import setup_logger


class ConfigDefaults(Enum):
    """Default settings for Vina. By default, assumes the binding pocket is centered at the origin."""

    center: List[float] = [0.0, 0.0, 0.0]
    box_size: List[float] = [25, 25, 25]
    exhaustiveness: int = 8
    n_poses: int = 10


def load_vina_config(config_file: str):
    """Loads a Vina configuration file (yaml format)."""
    with open(config_file, "r") as f:
        config = yaml.load(f, Loader=yaml.SafeLoader)

    # add defaults if not present
    for key, value in ConfigDefaults.__members__.items():
        if key not in config:
            config[key] = value.value

    return config


def convert_pdb_to_pdbqt(file: str):
    """Converts a PDB file to a PDBQT file."""
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("pdb", "pdbqt")

    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, file)

    # Add hydrogens
    mol.AddHydrogens()

    # Add Gasteiger charges
    mol.AddPolarHydrogens()
    for atom in openbabel.OBMolAtomIter(mol):
        atom.SetPartialCharge(0.0)
    openbabel.OBChargeModel.FindType("gasteiger").ComputeCharges(mol)

    output_file = str(Path(file).with_suffix(".pdbqt"))
    obConversion.WriteFile(mol, output_file)

    remove_tags = {"ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF"}
    with open(output_file) as file:
        filtered_lines = []
        for line in file.readlines():
            for tag in remove_tags:
                if tag in line:
                    break
            else:
                filtered_lines.append(line)

    with open(output_file, "w") as file:
        file.writelines(filtered_lines)

    return output_file


def run(sdf_file: str, receptor_file: str, config_file: str, output_dir: str):
    """Runs docking on a single SDF file."""
    # Convert the SDF file to PDBQT strings
    logger = setup_logger("DockSDFfile.log")

    cfg = load_vina_config(config_file)
    if cfg["center"] == [0.0, 0.0, 0.0]:
        logger.warning(
            f"Using default location of {cfg['center']} for the binding pocket, "
            f"specify 'center' in the config file to change this."
        )

    if Path(receptor_file).suffix == ".pdb":
        logger.info("Converting receptor PDB file to PDBQT file..")
        receptor_file = convert_pdb_to_pdbqt(receptor_file)
    elif Path(receptor_file).suffix != ".pdbqt":
        raise ValueError(
            f"Receptor file must be a PDB or PDBQT file, got {receptor_file}."
        )

    logger.info("Converting SDF file to PDBQT file..")
    mol_supplier = Chem.SDMolSupplier(sdf_file, removeHs=False)
    mol_preparer = MoleculePreparation()
    mol_pdbqt_strings = defaultdict(list)
    for i, mol in enumerate(mol_supplier):
        # Prepare the PDBQT inputs
        mol_setups = mol_preparer.prepare(mol)
        for setup in mol_setups:
            pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(setup)
            if is_ok:
                mol_pdbqt_strings[i].append(pdbqt_string)
            else:
                logger.error(
                    f"Could not prepare molecule {i+1}, got error: {error_msg}."
                )

    all_pdbqt_strings = [s for i in mol_pdbqt_strings for s in mol_pdbqt_strings[i]]

    # Define the output path for the docked complex
    docked_output_path = Path(f"{output_dir}/{Path(sdf_file).stem}_docked.pdbqt")

    # Run Vina for docking
    logger.info(f"Running Vina on: {receptor_file}...")
    v = Vina()
    v.set_receptor(receptor_file)
    v.set_ligand_from_string(all_pdbqt_strings)

    v.compute_vina_maps(center=cfg["center"], box_size=cfg["box_size"])
    v.dock(exhaustiveness=cfg["exhaustiveness"], n_poses=cfg["n_poses"])
    v.write_poses(str(docked_output_path), n_poses=cfg["n_poses"])


if __name__ == "__main__":
    run("sdf_output.sdf", "5f19.pdbqt", "vina_config.yaml", ".")
