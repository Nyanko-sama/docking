import argparse
import os
import sys
from glob import glob
from prody import *
from pathlib import Path
from utils import locate_file
import rdkit
import subprocess

print("rdkit version:", rdkit.__version__)

vina = locate_file(from_path=Path.cwd().parent, query_path=f'vina*', query_name='AutoDock Vina')


def dock_ligands(receptor_info : list[tuple[Path, Path]], ligands_folder, exhaustiveness=32, num_modes=9):
    if not os.path.exists('../data/docking_files'):
        os.makedirs('../data/docking_files')

    if not os.path.exists('../output'):
        os.mkdir('../output')
    
    if sys.platform == 'win32':
        ligands_path_cmd = ' --batch '.join(list(glob(f'..\data\{ligands_folder}\*.pdbqt')))
        
    else:
        ligands_path_cmd =  Path(f'../data/{ligands_folder}/*.pdbqt')
    
    ligands = [Path(p).stem for p in glob(f'..\data\{ligands_folder}\*.pdbqt')]
    print(f'Found {len(ligands)} ligand(s) to dock: {", ".join(ligands)}')
    print(f'Found {len(receptor_info)} receptor(s) to dock to')
    
    completed_dockings = 0
    skipped_dockings = 0
    
    for receptor, config in receptor_info:
        out_path = f'../output/{receptor.stem}'.removesuffix('_prepared')
        if not os.path.exists(out_path):
            os.makedirs(out_path, exist_ok=True)
        elif all([os.path.exists(f'{out_path}/{lig}_out.pdbqt') for lig in ligands]):
            print(f'Already finished {receptor.stem}, skipping..')
            skipped_dockings += len(ligands)
            continue

        print(f'Docking {len(ligands)} ligand(s) to {receptor.stem}...')
        command = ' '.join([
            str(vina),
            '--receptor', str(receptor),
            '--batch', str(ligands_path_cmd),
            '--config', str(config),
            f'--exhaustiveness={exhaustiveness}',
            f'--num_modes={num_modes}',
            '--cpu', str(os.cpu_count()),
            '--dir', out_path
        ])
        subprocess.run(command, shell=True, check=True)
        print(f'Completed docking {len(ligands)} ligand(s) to {receptor.stem}')
        print(f'Output saved to {out_path}')
        completed_dockings += len(ligands)
    
    total_dockings = completed_dockings + skipped_dockings
    print(f'\n=== Docking Summary ===')
    print(f'Total docking runs completed: {completed_dockings}')
    if skipped_dockings > 0:
        print(f'Total docking runs skipped (already done): {skipped_dockings}')
    print(f'Total docking runs: {total_dockings} ({len(receptor_info)} receptors × {len(ligands)} ligands)')

def run_docking(args):
    if not os.path.exists(args.dock_files_path):
        print(f"ERROR: Docking files directory '{args.dock_files_path}' does not exist.", file=sys.stderr)
        sys.exit(1)
    
    receptors = list(glob(f'{args.dock_files_path}/*.pdbqt'))
    receptors = [Path(r) for r in receptors]
    
    if len(receptors) == 0:
        print(f"WARNING: No receptor .pdbqt files found in '{args.dock_files_path}'", file=sys.stderr)
        print(f"Make sure prepare_receptors.py has been run successfully.", file=sys.stderr)
        return
    
    configs = [f'{args.dock_files_path}/{rec.stem}.box.txt' for rec in receptors]
    receptor_info = list(zip(receptors, configs))
    dock_ligands(receptor_info, args.ligands_path, exhaustiveness=args.exhaustiveness, num_modes=args.num_modes)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--dock_files_path', default='../data/docking_files', help='Path to the folder where the preprared receptor .pdbqt files and vina configs are stored.')
    parser.add_argument('-l', '--ligands_path', default='../data/prepared_ligands', help='Folder containing prepared ligands for docking.')
    parser.add_argument('-e', '--exhaustiveness', type=int, default=32, help='Search exhaustiveness for Vina (lower=faster, less accurate; default=32)')
    parser.add_argument('-n', '--num_modes', type=int, default=9, help='Number of binding modes to generate (lower=faster; default=9)')
    args = parser.parse_args()
    run_docking(args)
