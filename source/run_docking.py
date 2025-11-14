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
import subprocess

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
    ligandSDF = f"{ligand_name}_scrubbed.sdf"
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


def prepare_ligands(args) -> list[Path]:
    ligands = pd.read_csv(args.ligands_path, sep='\t')
    if not os.path.exists('../data/prepared_ligands'):
        os.makedirs('../data/prepared_ligands')
    lig_paths = []
    
    for name, smiles in zip(ligands['name'], ligands['smiles']):
        lig_paths.append(prepare_ligand(smiles, ligand_name=name))

    return lig_paths

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

    subprocess.run(command, check=True)

    return (Path(f'../data/docking_files/{pdb_path.stem}_prepared.pdbqt'), Path(f'../data/docking_files/{pdb_path.stem}_prepared.box.txt'))

def prepare_receptors(args) -> list[tuple[Path, Path]]:
    pdbs = list(Path(args.pdbs_path).rglob('*.pdb'))
    out = []
    for pdb in pdbs:
        name = pdb.stem
        try:
            pockets = pd.read_csv(f'{args.pocket_preds_path}/{name}.pdb_predictions.csv', skipinitialspace=True)
        except FileNotFoundError:
            raise FileExistsError(f'Predictions for {name} not found. Use run_p2rank.py to generate the pocket predictions first.')
        
        best = pockets.sort_values('score').iloc[0]
        out.append(prepare_receptor(pdb, (best['center_x'], best['center_y'], best['center_z']), (args.box_size, args.box_size, args.box_size)))
    
    return out

def dock_ligands(receptor_info : list[tuple[Path, Path]], lig_paths : list[Path]):
    if not os.path.exists('../data/docking_files'):
        os.makedirs('../data/docking_files')

    if not os.path.exists('../output'):
        os.mkdir('../output')

    for receptor, config in receptor_info:
        for ligand in lig_paths:
            subprocess.run([
                str(vina),
                '--receptor', str(receptor),
                '--ligand', str(ligand),
                '--config', str(config),
                f'--exhaustiveness={args.e}',
                '--out', f'../output/{receptor.stem}_{ligand.stem}'
            ])

def run_docking(args):
    lig_paths = prepare_ligands(args)
    receptor_info = prepare_receptors(args)
    dock_ligands(receptor_info, lig_paths)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pocket_preds_path', default='../data/p2rank_output', help='Output directory path')
    parser.add_argument('--pdbs_path', default='../data/pdbs', help='Path to the directory where .pdb files are stored')
    parser.add_argument('--ligands_path', default='../data/ligands.csv', help='Path to ligands')
    parser.add_argument('--ph', default=6, type=int, help='pH for ligand preparation')
    parser.add_argument('--skip_tautomer', action='store_true', help='Skip tautomers in ligand preparation')
    parser.add_argument('--skip_acidbase', action='store_true', help='Skip acidbase in ligand preparation')
    parser.add_argument('--box_size', type=int, default=40, help='Box size in angstroms for docking. Default is the Prankweb value.')
    parser.add_argument('-e', '--exhaustiveness', type=int, default=32, help='Search exhaustiveness for Vina')

    args = parser.parse_args()
    run_docking(args)