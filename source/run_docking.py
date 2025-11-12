import re
import argparse
import os
import pandas as pd
import sys, platform
from glob import glob
from prody import *
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit

print("rdkit version:", rdkit.__version__)

# Helper functions
def locate_file(from_path : Path = None, query_path = None, query_name = "query file"):

    if not from_path or not query_path:
        raise ValueError("Must specify from_path and query_path")


    possible_path = list(from_path.rglob(query_path))

    if not possible_path:
        raise FileNotFoundError(f"Cannot find {query_name} from {from_path} by {query_path}")

    return_which = (
        f"using {query_name} at:\n"
        f"{possible_path[0]}\n"
    )
    print(return_which)

    return possible_path[0]


# Commandline scripts
if sys.platform == 'linux':
    path_root = Path("/usr/local/bin")
elif sys.platform == 'win32':
    path_root = Path(sys.executable).parent
else:
    raise ValueError('Unsupported OS.')

scrub = locate_file(from_path = Path.cwd().parent, query_path = "scrub.py", query_name = "scrub.py")
mk_prepare_ligand = locate_file(from_path = path_root, query_path = "mk_prepare_ligand.py", query_name = "mk_prepare_ligand.py")
mk_prepare_receptor = locate_file(from_path = path_root, query_path = "mk_prepare_receptor.py", query_name = "mk_prepare_receptor.py")
mk_export = locate_file(from_path = path_root, query_path = "mk_export.py", query_name = "mk_export.py")

def prepare_ligand(ligand_smiles, ph = 6, skip_tautomer=False, skip_acidbase=False, ligand_name='test_ligand'):
    # Adapted from https://colab.research.google.com/drive/1cHSl78lBPUc_J1IZxLgN4GwD_ADmohVU?usp=sharing#scrollTo=qBQN6WzwvkGB
    args = ""
    if skip_tautomer:
        args += "--skip_tautomer "

    if skip_acidbase:
        args += "--skip_acidbase "

    ligand_name = re.sub(r'-', '_', ligand_name)
    ligandSDF = f"{ligand_name}_scrubbed.sdf"
    # Scrub the molecule 
    os.system(f'python {scrub} "{ligand_smiles}" -o {ligandSDF} --ph {ph} {args}')

    # Runs meeko mk_prepare_ligand with the following arguments
    os.system(f'python {mk_prepare_ligand} -i ../data/{ligandSDF} -o ../data/{ligand_name}.pdbqt')

def dock_ligands(args):
    ligands = pd.read_csv(args.ligands_path, sep='\t')
    pdbs = list(Path(args.pdbs_path).rglob('*.pdb'))
    for name, smiles in zip(ligands['name'], ligands['smiles']):
        prepare_ligand(smiles, ligand_name=name)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_path', default='../data/p2rank_output', help='Output directory path')
    parser.add_argument('--pdbs_path', default='../data/pdbs', help='Path to the directory where .pdb files are stored')
    parser.add_argument('--ligands_path', default='../data/ligands.csv', help='Path to ligands')
    parser.add_argument('--ph', default=6, type=int, help='pH for ligand preparation')
    parser.add_argument('--skip_tautomer', action='store_true', help='Skip tautomers in ligand preparation')
    parser.add_argument('--skip_acidbase', action='store_true', help='Skip acidbase in ligand preparation')

    args = parser.parse_args()
    dock_ligands(args)