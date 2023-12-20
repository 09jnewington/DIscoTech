import mainTakeuchi
import DockPDBFiles
import CombineComplex
directory_path = r"pdb_files"



## call main Takeuchi
mainTakeuchi.run()
## takes the output of the Mains and runs DockSDFFile on it
DockPDBFiles.run()
## Run combinecomplex
CombineComplex.run(directory_path)
