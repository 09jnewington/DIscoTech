import mainTakeuchi
import DockSDFfile
import CombineComplex

directory_path = r"C:\Users\joshi\OneDrive\Desktop\DiscoTech\pdb_files"


## call main Takeuchi
mainTakeuchi.run()
## takes the output of the Mains and runs DockSDFFile on it
DockSDFfile.run()
## Run combinecomplex
CombineComplex.run(directory_path)

