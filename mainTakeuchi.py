import sys
import random
from rdkit import Chem
import os
from rdkit.Chem import RDConfig
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer


#DTX dummy atoms are [6, 9, 18, 23, 32]

#List of common MedChem R groups from the Takeuchi Paper
Takeuchi_list = r'smiles.txt.txt'

def get_random_smiles(filename):
    #picks a random sidechain from a file
    with open(filename, 'r') as file:
        lines = file.readlines()
    return random.choice(lines).strip()

random_sidechain_smiles = Chem.MolFromSmiles(get_random_smiles(Takeuchi_list))


def get_highlighted_atoms(mol, important_indices=[]):
    #Take the user inputted important indices and stores a list of highlighted atoms
    connecting_atoms = set()
    for idx in important_indices:
        for jdx in important_indices:
            if idx != jdx:
                path = Chem.rdmolops.GetShortestPath(mol, idx, jdx)
                connecting_atoms.update(path)
    highlighted_atoms = set(important_indices).union(connecting_atoms)
    return list(highlighted_atoms)

def mark_and_replace_atoms(smiles, important_atom_indices, num_variants=5):
    modified_mols = []
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return [], "Invalid SMILES code."

    highlighted_atoms = get_highlighted_atoms(mol, important_atom_indices)
    mol = reduce_structure(mol, important_atom_indices)
    modified_mols.append(mol)

    for i in range(num_variants):
        fragment = Chem.MolFromSmiles(get_random_smiles(Takeuchi_list))  # Getting a new random fragment for each replacement
        available_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() not in ['H', 'F', 'Cl', 'Br']]
        replacement_connection_point = random.choice(available_atoms)
        mod_mol = Chem.ReplaceSubstructs(fragment, Chem.MolFromSmiles('[*]'), mol, replaceAll=True, replacementConnectionPoint=replacement_connection_point)
        modified_mols.append(mod_mol[0])

    return modified_mols, None

def reduce_structure(mol, important_indices=None):
    if important_indices is None:
        important_indices = []

    # Calculate connecting atoms
    connecting_atoms = set()
    for idx in important_indices:
        for jdx in important_indices:
            if idx != jdx:
                path = Chem.rdmolops.GetShortestPath(mol, idx, jdx)
                connecting_atoms.update(path)

    # Combined set of important and connecting atoms
    combined_important_atoms = set(important_indices).union(connecting_atoms)

    # Get ring info and determine which rings to keep
    ring_info = mol.GetRingInfo()
    rings_to_keep = set()
    
    for ring in ring_info.AtomRings():
        if any(atom_idx in ring for atom_idx in combined_important_atoms):
            rings_to_keep.update(ring)


    # Determine atoms to keep
    keep_atoms = combined_important_atoms.union(rings_to_keep)

    # Remove unneeded atoms
    remove_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIdx() not in keep_atoms]
    for idx in sorted(remove_atoms, reverse=True):
        mol = Chem.RWMol(mol)
        mol.RemoveAtom(idx)

    return mol

##Test

def SAScorer(mol_to_draw):
        l1 = len(mol_to_draw)
        s_0 = 0
        molstoremoveidx = []
        for idx, mol in enumerate(mol_to_draw):
            if mol is not None:
                Chem.GetSSSR(mol)
                if idx == 0:
                    s_0 = sascorer.calculateScore(mol)
                else:
                    s = sascorer.calculateScore(mol)
                    if s > s_0:
                        molstoremoveidx.append(idx)
        for idx, mol in enumerate(mol_to_draw):
            if idx in molstoremoveidx:
                mol_to_draw.remove(mol)
        print("number of molecules removed: ", l1 - len(mol_to_draw))
        return(mol_to_draw)

def run():
    # Prompt the user for inputs
    smiles = input("Enter SMILES: ")
    mol_count = Chem.MolFromSmiles(smiles)
    num_atoms = mol_count.GetNumAtoms()
    important_indices_str = input(f"The number of indices in your molecule is {num_atoms}, enter important atom indices (comma separated): ")
    num_variants_str = input("Number of modified variants: ")
    

    num_variants = 5 if not num_variants_str else int(num_variants_str)
    important_indices = list(map(int, important_indices_str.split(',')))

    modified_molecules, error = mark_and_replace_atoms(smiles, important_indices, num_variants)

    if error:
        print(f"Error: {error}")
        sys.exit(1)

    # Generate the image with original and modified molecules
    mols_to_draw = [Chem.MolFromSmiles(smiles)] + modified_molecules

                
    #mol_to_draw = SAScorer(mols_to_draw)

    with Chem.SDWriter('sdf_output.sdf') as writer:
        for mol in mols_to_draw:
            writer.write(mol)

    # Make pdb files
    for idx, mol in enumerate(mols_to_draw):
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        pdb_filename = f"pdb_files/molecule_{idx}.pdb"
        with open(pdb_filename, 'w') as file:
            file.write(Chem.MolToPDBBlock(mol))

    
 
    print("sdf file saved as sdf_output.sdf")

if __name__ == '__main__':
    run()

