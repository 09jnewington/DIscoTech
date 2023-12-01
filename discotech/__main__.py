from . import mainTakeuchi, DockSDFfile, CombineComplex
import argparse

directory_path = r"C:\Users\joshi\OneDrive\Desktop\DiscoTech\pdb_files"


def main():
    """Runs the entire pipeline."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d",
        "--directory_path",
        type=str,
        required=True,
        help="Directory containing PDB files.",
    )
    args = parser.parse_args()

    ## call main Takeuchi
    mainTakeuchi.run()
    ## takes the output of the Mains and runs DockSDFFile on it
    DockSDFfile.run()
    ## Run combinecomplex
    CombineComplex.run(args.directory_path)
