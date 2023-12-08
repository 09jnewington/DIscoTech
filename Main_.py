import mainTakeuchi
import DockSDFfile
import CombineComplex
import FinalSDF

directory_path = r"pdb_files"


## call main Takeuchi
mainTakeuchi.run()
## takes the output of the Mains and runs DockSDFFile on it
DockSDFfile.run()
## Run combinecomplex
CombineComplex.run(directory_path)
##Collect mols and add to final SDF file
FinalSDF.run()


