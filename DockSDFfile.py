import os
import subprocess

def run():

    input2 = input("ready to proceed with docking?")

    # Define paths and directories
    sdf_file = r"sdf_output.sdf"
    output_dir = r"pdb_files"

    # Convert the SDF file to individual PDB files
    subprocess.run(["obabel", sdf_file, "-O", os.path.join(output_dir, "molecule_.pdb"), "-m"])

    for pdb in os.listdir(output_dir):
        if pdb.endswith(".pdb"):
            pdb_path = os.path.join(output_dir, pdb)
            pdbqt_path = os.path.join(output_dir, pdb.replace('.pdb', '.pdbqt'))

            # Run the obabel command
            subprocess.run(["obabel", pdb_path, "-O", pdbqt_path, "--gen3d", "best", "-p"])

            # Define receptor path and config file path
            receptor_path = r"5LXT_removesolvent_removehet_hadd.pdbqt"
            config_path = r"conf.txt"

            # Define the output path for the docked complex
            docked_output_path = os.path.join(output_dir, pdb.replace('.pdb', '_docked.pdbqt'))

            # Run Vina for docking
            subprocess.run(["vina", "--receptor", receptor_path, "--ligand", pdbqt_path, "--out", docked_output_path, "--config", config_path])

if __name__ == '__main__':
    run()
