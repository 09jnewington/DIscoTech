import random
from enum import Enum
from typing import Optional, List

from pkg_resources import resource_filename
from rdkit import Chem
from rdkit.Chem import Mol, AllChem, rdmolops
from discotech.utils import setup_logger


# List of common MedChem R groups from the Takeuchi Paper
with open(resource_filename("discotech", "data/takeuchi_R_groups.txt"), "r") as file:
    TAKEUCHI_LIST = [l.strip() for l in file.readlines()]


class Valences(Enum):
    """Atomic valences for common elements."""

    C = 4
    N = 3
    O = 2


def get_highlighted_atoms(mol: Mol, important_indices: Optional[List[int]] = None):
    # Take the user inputted important indices and stores a list of highlighted atoms
    if important_indices is None:
        important_indices = []

    connecting_atoms = set()
    for idx in important_indices:
        for jdx in important_indices:
            if idx != jdx:
                path = Chem.rdmolops.GetShortestPath(mol, idx, jdx)
                connecting_atoms.update(path)
    highlighted_atoms = set(important_indices).union(connecting_atoms)
    return list(highlighted_atoms)


def mark_and_replace_atoms(
    smiles: str,
    important_atom_indices: Optional[List[int]] = None,
    num_variants: int = 5,
) -> List[Mol]:
    """Takes a SMILES string, marks the important atoms, and replaces them with random fragments."""
    modified_mols = []
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}.")

    mol = reduce_structure(mol, important_atom_indices)
    n_atoms = mol.GetNumAtoms()
    modified_mols.append(mol)

    n_molecules_made = 0
    while n_molecules_made < num_variants:
        fragment = Chem.MolFromSmiles(
            random.choice(TAKEUCHI_LIST)
        )  # Getting a new random fragment for each replacement
        available_atoms = [
            atom.GetIdx()
            for atom in mol.GetAtoms()
            if atom.GetSymbol() in Valences.__members__ and not atom.GetIsAromatic()
        ]
        replacement_atom_index = random.choice(available_atoms)

        combined_mol = Chem.CombineMols(mol, fragment)
        mod_mol = Chem.EditableMol(combined_mol)
        replacement_atom = combined_mol.GetAtomWithIdx(replacement_atom_index)

        first_fragment_atom = n_atoms

        # if the replacement atom valence is full, break a random bond
        replacement_valence = replacement_atom.GetExplicitValence()

        mod_mol.AddBond(
            replacement_atom_index,
            first_fragment_atom,
            order=Chem.BondType.SINGLE,
        )

        if (
            replacement_valence
            >= Valences.__members__[replacement_atom.GetSymbol()].value
        ):
            bond_indices = [
                bond.GetIdx()
                for bond in replacement_atom.GetBonds()
                if bond.GetEndAtomIdx() != first_fragment_atom
            ]
            if not bond_indices:
                raise ValueError(
                    f"Atom {replacement_atom_index} ({replacement_atom.GetSymbol()}) has full valence but "
                    f"no bonds could be found."
                )
            bond_index = random.choice(bond_indices)
            bond_to_break = combined_mol.GetBondWithIdx(bond_index)
            bond_type = bond_to_break.GetBondType()

            mod_mol.RemoveBond(
                replacement_atom_index,
                bond_to_break.GetOtherAtomIdx(replacement_atom_index),
            )

            if bond_type == Chem.rdchem.BondType.DOUBLE:
                mod_mol.AddBond(
                    replacement_atom_index,
                    bond_to_break.GetOtherAtomIdx(replacement_atom_index),
                    order=Chem.BondType.SINGLE,
                )
            elif bond_type == Chem.rdchem.BondType.TRIPLE:
                mod_mol.AddBond(
                    replacement_atom_index,
                    bond_to_break.GetOtherAtomIdx(replacement_atom_index),
                    order=Chem.BondType.DOUBLE,
                )

        mod_mol = mod_mol.GetMol()

        # set the attachment atom to a carbon
        attachment_point_atom = mod_mol.GetAtomWithIdx(first_fragment_atom)
        attachment_point_atom.SetAtomicNum(6)

        Chem.SanitizeMol(mod_mol)
        AllChem.AddHs(mod_mol)
        AllChem.ComputeGasteigerCharges(mod_mol, throwOnParamFailure=True)

        # Generate 3D conformation
        AllChem.EmbedMolecule(mod_mol, randomSeed=0)
        AllChem.MMFFOptimizeMolecule(mod_mol)

        for atom in mod_mol.GetAtoms():
            charge = atom.GetProp("_GasteigerCharge").lower()
            if "nan" in charge or "inf" in charge:
                continue

        modified_mols.append(mod_mol)
        n_molecules_made += 1

    return modified_mols


def reduce_structure(mol: Mol, important_indices: Optional[List[int]] = None) -> Mol:
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

    return mol


def run(
    smiles: str,
    preserve_indices: Optional[List[int]] = None,
    n: int = 5,
    outfile: str = "sdf_output.sdf",
):
    logger = setup_logger("mainTakeuchi.log")

    logger.info("Generating new molecules...")
    # Prompt the user for inputs
    modified_molecules = mark_and_replace_atoms(smiles, preserve_indices, n)

    logger.info(f"Writing new molecules to {outfile}...")
    mols_to_write = [Chem.MolFromSmiles(smiles)] + modified_molecules
    with Chem.SDWriter(outfile) as writer:
        for mol in mols_to_write:
            writer.write(mol)

    logger.info("Molecule generation finished successfully.")


if __name__ == "__main__":
    run("CC(=O)OC1=CC=CC=C1C(=O)O")
