# DiscoTech
DiscoTech Biosciences drug discovery platform

DiscoTech is a platfom that takes a natural product drug as a user input, along with the atoms in that natural compound that are necessary for binding and produces
a set of drug-like small molecules that contain only the necessry atoms from the natural compound (and the atoms needed to connect them together) augmented  with common R 
groups in medicinal chemistry (Takeuchi K, Kunimoto R, Bajorath J. Systematic mapping of R-group space enables the generation of an R-group replacement system for medicinal chemistry. Eur J Med Chem. 2021;225:113771. doi:10.1016/j.ejmech.2021.113771). This set is then docked and the subset that preserve all the protein interactions as the starting natural product is returned as a final sdf file.

**To run this code you will need to :**
- Install the necessary python modules
- Download AutoDock Vina and have vina.exe and vina_split.exe in the directory
- conf.txt which contains the docking details (binding site and number of modes) provided here for discodermolide in tubulin (5LXT)
- A pdb file of your protein with the ligand removed (provided here for discodermolide in tubulin (5LXT)) This can be done in PyMol.

  Run Main_.py to run start the program.


