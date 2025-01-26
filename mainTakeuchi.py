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
#C[C@H]1[C@@H](OC(=O)[C@@H]([C@H]1O)C)C[C@@H](/C=C\[C@H](C)[C@@H]([C@@H](C)/C=C(/C)\C[C@H](C)[C@H]([C@H](C)[C@H]([C@@H](C)/C=C\C=C)OC(=O)N)O)O)O

#List of common MedChem R groups from the Takeuchi Paper
Takeuchi_list = r'smiles.txt.txt'

def get_random_smiles(filename):
    #picks a random sidechain from a file
    with open(filename, 'r') as file:
        lines = file.readlines()
    return random.choice(lines).strip()


def mark_and_replace_atoms(smiles, important_atom_indices, num_variants=5):
    modified_mols = []
    mol = Chem.MolFromSmiles(smiles)
    modified_mols.append((mol))
    if mol is None:
        return [], "Invalid SMILES code."
    

    mol_noh = reduce_structure(mol, important_atom_indices)
    modified_mols.append((mol_noh))

    mol_noh = Chem.RemoveHs(mol_noh)
    for i in range(num_variants):
        fragment = Chem.MolFromSmiles(get_random_smiles(Takeuchi_list))  # Getting a new random fragment for each replacement
        # Define the SMARTS pattern for the reaction
        reaction_smarts = '[*:1].[#0]-[!#0:3] >> [*:1]-[!#0:3]'
        # Create the reaction from the SMARTS pattern
        reaction = AllChem.ReactionFromSmarts(reaction_smarts)
        ps = reaction.RunReactants((mol_noh, fragment))
        if ps:
            # Flatten the list of tuples
            flat_ps = [product for sublist in ps for product in sublist]
            # Pick a random product
            mol = random.choice(flat_ps)
        Chem.SanitizeMol(mol)
        modified_mols.append((mol))
        print(Chem.MolToSmiles(mol))

    return modified_mols, None

def prepare_mol(mol):
    """Prepares a molecule for docking, adding hydrogens, partial charges, and 3D coordinates."""
    Chem.SanitizeMol(mol)
    mol = AllChem.AddHs(mol)
    AllChem.ComputeGasteigerCharges(mol, throwOnParamFailure=True)

    # Generate 3D conformation
    AllChem.EmbedMolecule(mol, randomSeed=0)
    AllChem.MMFFOptimizeMolecule(mol)

    return mol

def reduce_structure(mol, important_indices):
    """
    Reduces the structure to only the important atoms and their connecting atoms.

    If no important indices are provided, all atoms are considered important.
    """
    if important_indices is not None:
        # Calculate connecting atoms
        connecting_atoms = set()
        for idx in important_indices:
            for jdx in important_indices:
                if idx != jdx:
                    path = Chem.rdmolops.GetShortestPath(mol, idx, jdx)
                    connecting_atoms.update(path)

        # Combined set of important and connecting atoms
        combined_important_atoms = set(important_indices).union(connecting_atoms)
    else:
        combined_important_atoms = set([atom.GetIdx() for atom in mol.GetAtoms()])

    # Get ring info and determine which rings to keep
    ring_info = mol.GetRingInfo()
    rings_to_keep = set()

    for ring in ring_info.AtomRings():
        if any(atom_idx in ring for atom_idx in combined_important_atoms):
            rings_to_keep.update(ring)

    # Determine atoms to keep
    keep_atoms = combined_important_atoms.union(rings_to_keep)

    # Remove unneeded atoms
    remove_atoms = [
        atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIdx() not in keep_atoms
    ]
    mol = Chem.RWMol(mol)
    for idx in sorted(remove_atoms, reverse=True):
        mol.RemoveAtom(idx)
    mol = mol.GetMol()

    redmol = Chem.MolToSmiles(mol)
    redmol = Chem.MolFromSmiles(redmol)


    print(type(redmol))
    return redmol



def SAScorer(mol_to_draw):
    l1 = len(mol_to_draw)
    s_0 = 0
    molstoremoveidx = []

    # Assuming s_0 is a threshold score, you should define it here
    # For example, s_0 = 5 or any other value that makes sense in your context

    for idx, mol in enumerate(mol_to_draw):
        if mol is not None:
            Chem.GetSSSR(mol)
            s = sascorer.calculateScore(mol)
            mol.SetProp('SAScore', str(s))
            print(s)
            if idx == 0:
                s_0 = s
                print("s_0: ", s_0)
            if s > s_0:
                molstoremoveidx.append(idx)
                print(f"s: {s} > s_0: {s_0} - removing molecule {idx}")

    # Create a new list for the molecules to keep
    mol_to_keep = [mol for idx, mol in enumerate(mol_to_draw) if idx not in molstoremoveidx]

    print("number of molecules removed: ", l1 - len(mol_to_keep))
    return mol_to_keep

def run():
    # Prompt the user for inputs
    smiles = ('C[C@H]1[C@@H](OC(=O)[C@@H]([C@H]1O)C)C[C@@H](/C=C\[C@H](C)[C@@H]([C@@H](C)/C=C(/C)\C[C@H](C)[C@H]([C@H](C)[C@H]([C@@H](C)/C=C\C=C)OC(=O)N)O)O)O')  #input("Enter SMILES: ")
    mol_count = Chem.MolFromSmiles(smiles)
    num_atoms = mol_count.GetNumAtoms()
    important_indices_str = '6,9,18,23,32' ##input(f"The number of indices in your molecule is {num_atoms}, enter important atom indices (comma separated): ")
    num_variants_str = '5' #input("Number of modified variants: ")
    

    num_variants = 5 if not num_variants_str else int(num_variants_str)
    important_indices = list(map(int, important_indices_str.split(',')))

    modified_molecules, error = mark_and_replace_atoms(smiles, important_indices, num_variants)

    if error:
        print(f"Error: {error}")
        sys.exit(1)

    # Generate the image with original and modified molecules
    mols_to_draw = [Chem.MolFromSmiles(smiles)] + modified_molecules

                
    mols_to_draw = SAScorer(mols_to_draw)

    mols_to_draw = [prepare_mol(mol) for mol in mols_to_draw]

    with Chem.SDWriter('sdf_output.sdf') as writer:
        for mol in mols_to_draw:
            writer.write(mol)
            

    # Make pdb files
    for idx, mol in enumerate(mols_to_draw):
        try:
            pdb_filename = f"pdb_files/molecule_{idx}.pdb"
            with open(pdb_filename, 'w') as file:
                file.write(Chem.MolToPDBBlock(mol))
                print(f"pdb file saved as {pdb_filename}")
        except:
            print("Error creating pdb file")
            continue
        
 
    print("sdf file saved as sdf_output.sdf")

if __name__ == '__main__':
    run()
