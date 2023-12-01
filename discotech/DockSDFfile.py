import os
import subprocess
from pathlib import Path

from .utils import setup_logger


def run(sdf_file: str, receptor_file: str, config_file: str, output_dir: str):
    """Runs docking on a single SDF file."""
    # Convert the SDF file to individual PDB files
    logger = setup_logger("DockSDFfile.log")

    pdb_path = os.path.join(output_dir, Path(sdf_file).with_suffix(".pdb").name)
    logger.info("Converting SDF file to PDB files...")
    subprocess.run(["obabel", sdf_file, "-O", pdb_path, "-m"])

    pdbqt_path = Path(pdb_path).with_suffix(".pdbqt")

    # Run the obabel command
    logger.info(f"Running open babel on: {pdb_path}")
    subprocess.run(["obabel", pdb_path, "-O", str(pdbqt_path), "--gen3d", "best", "-p"])

    # Define the output path for the docked complex
    docked_output_path = pdbqt_path.with_suffix(".pdbqt")

    # Run Vina for docking
    logger.info(f"Running Vina on: {str(pdbqt_path)}...")
    subprocess.run(
        [
            "vina",
            "--receptor",
            receptor_file,
            "--ligand",
            str(pdbqt_path),
            "--out",
            str(docked_output_path),
            "--config",
            config_file,
        ]
    )
