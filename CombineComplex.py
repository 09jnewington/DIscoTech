import os
from pymol import cmd
from openbabel import openbabel
from plip.structure.preparation import PDBComplex
import os
import pandas as pd
from io import StringIO
from rdkit import Chem
import re

# Define the path to your PDBQT file
#pdbqt_file_pathway = r"C:\Users\joshi\OneDrive\Desktop\DiscoTech\pdb_files\molecule_10_docked.pdbqt"
directory_path = r"pdb_files"

#Interaction types
interaction_patterns = {
    'Salt Bridge': r"saltbridge\((.*?)\)",
    'Hydrogen Bond': r"hbond\((.*?)\)",
    'Hydrophobic Interaction': r"hydroph_interaction\((.*?)\)",
    'pistacking': r"pistacking\((.*?)\)",
    'pication': r"pication\((.*?)\)",
    'halogen': r"halogen\((.*?)\)",
    'metal': r"metal\((.*?)\)"
}


def run(directory_path):
    # Get a list of .pdbqt files in the directory
    pdbqt_files = [f for f in os.listdir(directory_path) if f.endswith("_docked.pdbqt")]

    # Check if there are any .pdbqt files
    if not pdbqt_files:
        print("No .pdbqt files found in the directory.")
        return

    # Path of the first .pdbqt file
    first_file_path = os.path.join(directory_path, pdbqt_files[0])

    # Extract and store the binding site information of the first file
    first_file_interactions = extract_binding_site(first_file_path)
    
    # Loop over all .pdbqt files in the directory
    for filename in pdbqt_files[1:]:
        pdbqt_file_path = os.path.join(directory_path, filename)
        compare_and_process_file(pdbqt_file_path, first_file_interactions)

def parse_and_print_interactions(output):

    interaction_data = []
    for interaction_type, pattern in interaction_patterns.items():
        matches = re.findall(pattern, output)
        if matches:
            #print(f"{interaction_type}s:")
            #interaction_data.append(f"{interaction_type}s:")
            for match in matches:
                # Parse the necessary components from the match string
                bsatom_orig_idx_match = re.search(r'bsatom_orig_idx=(\d+)', match)
                ligatom_orig_idx_match = re.search(r'ligatom_orig_idx=(\d+)', match)
                reschain_match = re.search(r'reschain=\'(\w)\'', match)
                restype_match = re.search(r'restype=\'(\w+)\'', match)
                resnr_match = re.search(r'resnr=(\d+)', match)

                if not (bsatom_orig_idx_match and ligatom_orig_idx_match and reschain_match and restype_match and resnr_match):
                    continue

                # Extract the necessary components from the match string
                bsatom_orig_idx = bsatom_orig_idx_match.group(1)
                ligatom_orig_idx = ligatom_orig_idx_match.group(1)
                reschain = reschain_match.group(1)
                restype = restype_match.group(1)
                resnr = resnr_match.group(1)

                # append the parsed information
                interaction_data.append(f"reschain={reschain}, restype={restype}, resnr={resnr}")
                print("test3")
    return interaction_data

def extract_binding_site(file_path):
    # This should return the binding site information as a set

    with open(file_path, 'r',  encoding='utf-8') as file:
        content = file.readlines()

# Replace 'UNL' with '7AK' in each line
    new_content = []
    for line in content:
        if 'UNL' in line:
            new_line = line.replace('UNL', '7AK')
            new_content.append(new_line)
        else:
            new_content.append(line)

# Write the modified content back to the file
    with open(file_path, 'w',  encoding='utf-8') as file:
        file.writelines(new_content)

    print("File updated successfully.")

    cmd.reinitialize()        
    cmd.load(file_path, "ligand")
    protein_file = r'C:\Users\joshi\OneDrive\Desktop\DiscoTech\5LXT_no_7AK.pdb'
    complex_file =   r'C:\Users\joshi\OneDrive\Desktop\DiscoTech\complex4.pdb'
    cmd.load(protein_file, "prot")
    cmd.create("complex", "ligand, prot")
    cmd.save(complex_file, "complex")

#extracts interactions of ligand and protein
# Create a PDBComplex object
    mol = PDBComplex()
# Load the PDB file
    mol.load_pdb(complex_file)
# Analyze the loaded PDB file
    mol.analyze()

# Get the ligands information
    longnames = [x.longname for x in mol.ligands]
    bsids = [":".join([x.hetid, x.chain, str(x.position)]) for x in mol.ligands]

    # Find indices of ligands with the name '7AK'
    indices = [j for j, x in enumerate(longnames) if x == '7AK']
    # Print interactions for each '7AK' ligand found
    for idx in indices:
        bsid = bsids[idx]
        interactions = mol.interaction_sets[bsid]
        interactions_output = parse_and_print_interactions(str(interactions.all_itypes))
    interaction_output_set = set(interactions_output)
    return(interaction_output_set)

def compare_and_process_file(file_path, reference_interactions):
    # Process each file and compare its binding site
    current_file_interactions = extract_binding_site(file_path)
  
    print("Current file interactions: ")
    print(current_file_interactions)
    print("Reference interactions: ")
    print(reference_interactions)

    # Calculating the number of common interactions
    common_interactions = reference_interactions.intersection(current_file_interactions)
    num_common_interactions = len(common_interactions)

    # Calculating the total number of interactions for Ligand 1
    total_ligand1_interactions = len(reference_interactions)

    # Calculating the conservation score
    conservation_score = (num_common_interactions / total_ligand1_interactions) * 100

    print(conservation_score)


    if conservation_score < 90:
        print("Less than 90 percent conservation, deleting file")
        os.remove(file_path)

    else:
        print("Interactions is a subset, file saved")
    
if __name__ == '__main__':

    run(directory_path)
