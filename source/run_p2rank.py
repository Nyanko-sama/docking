import os
import argparse
import glob
from molscrub import Scrub
from rdkit import Chem
import sys
import pandas as pd
import re

def create_ds_file(pdbs_path):
    paths = glob.glob(f'{pdbs_path}/*.pdb')
    with open('../data/pdbs.ds', 'w') as f:
        f.write('\n'.join(paths))

def run_p2rank(p2rank_path, ds_path, outputd_dir):
    if not os.path.exists(outputd_dir):
        os.mkdir(outputd_dir)

    os.system(f'{p2rank_path}\\prank predict {ds_path} -o {outputd_dir}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--p2rank_path', default='..\\p2rank', help='p2rank root directory')
    parser.add_argument('--output_path', default='../data/p2rank_output', help='Output directory path')
    parser.add_argument('--pdbs_path', default='../data/pdbs', help='Path to the directory where .pdb files are stored')
    parser.add_argument('--ligands_path', default='../data/ligands.csv', help='Path to ligands')
    parser.add_argument('--ph', default=6, type=int, help='pH for ligand preparation')
    parser.add_argument('--skip_tautomer', action='store_true', help='Skip tautomers in ligand preparation')
    parser.add_argument('--skip_acidbase', action='store_true', help='Skip acidbase in ligand preparation')

    args = parser.parse_args()
    create_ds_file(args.pdbs_path)
    run_p2rank(args.p2rank_path, args.pdbs_path, args.output_path)
    