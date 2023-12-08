from rdkit import Chem
from rdkit.Chem import RDConfig
from rdkit.Chem import Descriptors
import os
import sys
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
# now you can import sascore
import sascorer
output_file = r"final_sdf_output.sdf"  # Name of the output SDF file


def run():
    supplier = Chem.SDMolSupplier(output_file)


    for mol in supplier:
            if mol is not None:
                s = sascorer.calculateScore(mol)
                print(s)
            else:
                print("Encountered a None molecule, skipping...")


if __name__ == '__main__':
    run()


