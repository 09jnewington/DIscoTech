import os
import subprocess
import AutoDockTools
from AutoDockTools.Utilities24 import prepare_ligand4
from AutoDockTools.Utilities24 import prepare_receptor4
from rdkit import Chem
script_path = prepare_ligand4.__file__


def extract_docking_score(pdb):
    for line in pdb:
        if line.startswith("REMARK VINA RESULT"):
            parts = line.split()
            score = parts[3]  # The score is typically the fourth item in this line
            return float(score)
    return None

def add_scores_to_existing_sdf(scores, sdf_file):
    """Assign docking scores from a list to molecules in an SDF file and overwrite the file."""
    # Load the SDF file
    suppl = Chem.SDMolSupplier(sdf_file)
    mols = [mol for mol in suppl if mol is not None]

    # Check if the number of scores matches the number of molecules
    if len(scores) != len(mols):
        raise ValueError("The number of scores does not match the number of molecules in the SDF file.")

    # Assign scores to molecules
    for score, mol in zip(scores, mols):
        mol.SetProp("DockingScore", str(score))

    # Overwrite the SDF file with the modified molecules
    writer = Chem.SDWriter(sdf_file)
    for mol in mols:
        writer.write(mol)
    writer.close()
    print("scores added to SDF file")

def run(): 

    # Define paths and directories
    sdf_file = r"sdf_output.sdf"
    output_dir = r"pdb_files"
    scores = []

    for pdb in os.listdir(output_dir):
        if pdb.endswith(".pdb"):
            pdb_path = os.path.join(output_dir, pdb)
            pdbqt_path = os.path.join(output_dir, pdb.replace('.pdb', '.pdbqt'))

            # Run the obabel command
            subprocess.run(['python',script_path, "-l", pdb_path, "-o", pdbqt_path])

            # Define receptor path and config file path
            receptor_path = r"5LXT_removesolvent_removehet_hadd.pdbqt"
            config_path = r"conf.txt"

            # Define the output path for the docked complex
            docked_output_path = os.path.join(output_dir, pdb.replace('.pdb', '_docked.pdbqt'))

            # Run Vina for docking
            subprocess.run(["vina", "--receptor", receptor_path, "--ligand", pdbqt_path, "--out", docked_output_path, "--config", config_path])
            
            with open(docked_output_path, "r") as docked_output:
                score = extract_docking_score(docked_output)
                scores.append(score)
    add_scores_to_existing_sdf(scores,sdf_file)




if __name__ == '__main__':
    run()




