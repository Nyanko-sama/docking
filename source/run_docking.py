import re
import argparse
import os
import pandas as pd
import sys
import numpy as np
from Bio.PDB import PDBParser
from glob import glob
from prody import *
from pathlib import Path
from utils import locate_file
import rdkit
import subprocess

print("rdkit version:", rdkit.__version__)


# Commandline scripts

if sys.platform == 'linux':
    path_root = Path("/usr/local/bin")
elif sys.platform == 'win32':
    path_root = Path(sys.executable).parent
else:
    raise ValueError('Unsupported OS.')

scrub = locate_file(from_path = Path.cwd().parent, query_path = "scrub.py", query_name = "scrub.py")
vina = locate_file(from_path=Path.cwd().parent, query_path=f'vina*', query_name='AutoDock Vina')
mk_prepare_ligand = locate_file(from_path = path_root, query_path = "mk_prepare_ligand.py", query_name = "mk_prepare_ligand.py")
mk_prepare_receptor = locate_file(from_path = path_root, query_path = "mk_prepare_receptor.py", query_name = "mk_prepare_receptor.py")
mk_export = locate_file(from_path = path_root, query_path = "mk_export.py", query_name = "mk_export.py")

def prepare_ligand(ligand_smiles, ph = 6, skip_tautomer=False, skip_acidbase=False, ligand_name='test_ligand') -> Path:
    # Adapted from https://colab.research.google.com/drive/1cHSl78lBPUc_J1IZxLgN4GwD_ADmohVU?usp=sharing#scrollTo=qBQN6WzwvkGB
    args = ""
    if skip_tautomer:
        args += "--skip_tautomer "

    if skip_acidbase:
        args += "--skip_acidbase "

    ligand_name = re.sub(r'-', '_', ligand_name)
    ligandSDF = f"../data/temp/{ligand_name}_scrubbed.sdf"
    output_path = f'../data/prepared_ligands/{ligand_name}.pdbqt'
    # Scrub the molecule
    subprocess.run([
        'python', str(scrub), ligand_smiles, '-o', ligandSDF, '--ph', str(ph) + args
    ], check=True)
    #os.system(f'python {scrub} "{ligand_smiles}" -o {ligandSDF} --ph {ph} {args}')

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
        if args.skip_tautomer:
            scrub_additional_args += '--skip_tautomer '
        if args.skip_acidbase:
            scrub_additional_args += '--skip_acidbase '

        # Scrub the molecules
        subprocess.run([
            'python', str(scrub), smi_path, '-o', sdf_path, '--ph', str(args.ph) + scrub_additional_args
        ], check=True, shell=True)

        subprocess.run([
            'python', str(mk_prepare_ligand), '-i', sdf_path, '--multimol_outdir', output_path
        ], check=True, shell=True)
        
    else:
        prepare_ligand(ligands['smiles'], ligand_name=ligands['name'])

    return ligands['name']

def prepare_receptor(pdb_path : Path, center_coords, box_sizes) -> list[Path]:
    # Export receptor atoms
    center_x, center_y, center_z = center_coords
    size_x, size_y, size_z = box_sizes

    command = [
        "python",
        str(mk_prepare_receptor),
        "-i", str(pdb_path),
        "-o", f'../data/docking_files/{pdb_path.stem}_prepared',
        "-p",  # Generate PDBQT file
        "-v",  # Generate Vina config
        "--box_center", str(center_x), str(center_y), str(center_z),
        "--box_size", str(size_x), str(size_y), str(size_z)
    ]

    subprocess.run(command, check=True, shell=True)

    return (Path(f'../data/docking_files/{pdb_path.stem}_prepared.pdbqt'), Path(f'../data/docking_files/{pdb_path.stem}_prepared.box.txt'))

def find_best_pocket(pdb, pockets):
    return pockets.sort_values('score').iloc[0]

def prepare_receptors(args) -> list[tuple[Path, Path]]:
    pdbs = list(Path(args.pdbs_path).rglob('*.pdb'))
    parser = PDBParser()
    out = []
    for pdb in pdbs:
        name = pdb.stem
        try:
            pockets = pd.read_csv(f'{args.pocket_preds_path}/{name}.pdb_predictions.csv', skipinitialspace=True)
        except FileNotFoundError:
            raise FileExistsError(f'Predictions for {name} not found. Use run_p2rank.py to generate the pocket predictions first.')
        
        best = find_best_pocket(pdb, pockets)
        molecule = parser.get_structure(pdb.stem, pdb)
        residues = list(molecule.get_residues())
        center = best['center_x'], best['center_y'], best['center_z']
        center_np = np.asarray(center)
        max_dist = 0

        # Find the furthest atom from the center for correct box size
        for residue in best['residue_ids'].split(' '):
            idx = int(residue[2:])
            max_dist_atom = 0
            for atom in residues[idx].get_atoms():
                c_dist = np.sum((atom.get_coord() - center_np) ** 2)
                if c_dist > max_dist_atom:
                    max_dist_atom = c_dist

            if max_dist_atom > max_dist:
                max_dist = max_dist_atom

        box_size = np.sqrt(max_dist) + 2 # 2A padding
        out.append(prepare_receptor(pdb, center, (box_size, box_size, box_size)))
    
    return out

def dock_ligands(receptor_info : list[tuple[Path, Path]], lig_names : list[Path]):
    if not os.path.exists('../data/docking_files'):
        os.makedirs('../data/docking_files')

    if not os.path.exists('../output'):
        os.mkdir('../output')

    for receptor, config in receptor_info:
        for ligand in lig_names:
            for path in glob(f'../data/prepared_ligands/{ligand}*.pdbqt'):
                lig_stem = Path(path).stem
                out_path = f'../output/{receptor.stem}_{lig_stem}.pdbqt'
                print(f'Docking {lig_stem} to {receptor.stem}...')
                command = ' '.join([
                    str(vina),
                    '--receptor', str(receptor),
                    '--ligand', str(path),
                    '--config', str(config),
                    f'--exhaustiveness={args.exhaustiveness}',
                    '--out', out_path
                ])
                subprocess.run(command, shell=True, check=True)
                print(f'Output saved to {out_path}')

def run_docking(args):
    lig_names = prepare_ligands(args)
    receptor_info = prepare_receptors(args)
    dock_ligands(receptor_info, lig_names)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pocket_preds_path', default='../data/p2rank_output', help='P2rank pocket predictions output path')
    parser.add_argument('--pdbs_path', default='../data/pdbs', help='Path to the directory where .pdb files are stored')
    parser.add_argument('--ligands_path', default='../data/ligands.csv', help='Path to the ligands csv file. Should contain fields "smiles" and "name".')
    parser.add_argument('--ph', default=7, type=int, help='pH for ligand preparation')
    parser.add_argument('--skip_tautomer', action='store_true', help='Skip tautomers in ligand preparation')
    parser.add_argument('--skip_acidbase', action='store_true', help='Skip acid/base conjugates in ligand preparation')
    # parser.add_argument('--box_size', type=int, default=20, help='Box size in Å for docking. Prankweb has a default of 40Å, here 20Å is used, same as default in Vina.')
    parser.add_argument('-e', '--exhaustiveness', type=int, default=32, help='Search exhaustiveness for Vina')

    args = parser.parse_args()
    run_docking(args)