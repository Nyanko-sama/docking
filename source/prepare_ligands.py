import argparse
import re
import sys
import os
import subprocess
import pandas as pd
from pathlib import Path
from utils import locate_file, get_path_root

mk_prepare_ligand = locate_file(from_path = get_path_root(), query_path = "mk_prepare_ligand.py", query_name = "mk_prepare_ligand.py")
scrub = locate_file(from_path = Path.cwd().parent, query_path = "scrub.py", query_name = "scrub.py")

def prepare_ligand(ligand_smiles, ph = 7, skip_tautomer=False, skip_acidbase=False, ligand_name='test_ligand') -> Path:
    # Adapted from https://colab.research.google.com/drive/1cHSl78lBPUc_J1IZxLgN4GwD_ADmohVU?usp=sharing#scrollTo=qBQN6WzwvkGB
    args = " "
    if skip_tautomer:
        args += "--skip_tautomers "

    if skip_acidbase:
        args += "--skip_acidbase "

    ligand_name = re.sub(r'-', '_', ligand_name)
    ligandSDF = f"../data/temp/{ligand_name}_scrubbed.sdf"
    output_path = f'../data/prepared_ligands/{ligand_name}.pdbqt'
    # Scrub the molecule
    subprocess.run([
        'python', str(scrub), ligand_smiles, '-o', ligandSDF, '--ph', str(ph), args
    ], check=True)

    # Runs meeko mk_prepare_ligand with the following arguments
    subprocess.run([
        'python', str(mk_prepare_ligand), '-i', f'../data/{ligandSDF}', '-o', output_path
    ], check=True)

    return Path(output_path)


def standardize_name(ligand_name):
    ligand_name = re.sub(r'\'"', '', ligand_name)
    return re.sub(r'[-, ]', '_', ligand_name)

def prepare_ligands(args):
    ligands = pd.read_csv(args.ligands_path, sep='\t')
    ligands = ligands.reindex(['smiles', 'name'], axis=1)
    ligands['name'] = ligands['name'].apply(standardize_name)

    batch = False
    if ligands.shape[0] > 1:
        batch = True

    if not os.path.exists('../data/prepared_ligands'):
        os.makedirs('../data/prepared_ligands')

    if not os.path.exists('../data/temp'):
        os.makedirs('../data/temp')
    
    if batch:
        smi_path = '../data/temp/ligands.smi'
        sdf_path = '../data/temp/ligand_batch.sdf'
        ligands.to_csv(smi_path, sep=' ', header=False, index=False)
        output_path = f'../data/prepared_ligands/'
        scrub_additional_args = ''
        if args.skip_tautomers:
            scrub_additional_args += '--skip_tautomers '
        if args.skip_acidbase:
            scrub_additional_args += '--skip_acidbase '

        # Scrub the molecules
        print('running scrub')
        scrub_cmd = ['python', str(scrub), smi_path, '-o', sdf_path, '--ph', str(args.ph)]
        if len(scrub_additional_args) > 0:
            scrub_cmd.extend(scrub_additional_args.split(' '))
        
        subprocess.run(scrub_cmd, check=True)
        print('running mk')
        subprocess.run([
            'python', str(mk_prepare_ligand), '-i', sdf_path, '--multimol_outdir', output_path
        ], check=True)
        
    else:
        prepare_ligand(ligands['smiles'], args.ph, args.tautomers, args.skip_acidbase, ligands['name'])

    return ligands['name']

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--ligands_path', default='../data/ligands.csv', help='Path to the ligands csv file. Should contain fields "smiles" and "name".')
    parser.add_argument('--ph', default=7, type=int, help='pH for ligand preparation')
    parser.add_argument('--tol', default=20, type=int, help='Maximum pocket center distance from the centroid of residues matched to the MSA indices for membrane pocket selection. Pockets further away are not considered when determining the best pocket.')
    parser.add_argument('--skip_acidbase', action='store_true', help='Skip acid/base conjugates in ligand preparation')
    parser.add_argument('--skip_tautomers', action='store_true', help='Skip tautomers in ligand preparation.')
    args = parser.parse_args()
    prepare_ligands(args)
