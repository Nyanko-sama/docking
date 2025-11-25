import re
import argparse
import os
import pandas as pd
import sys
import numpy as np
from Bio.PDB import PDBParser
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.PDB import Structure
from glob import glob
from prody import *
from pathlib import Path
from utils import locate_file, l2_norm
from pdbfixer import PDBFixer
from openmm.app import PDBFile
import rdkit
import subprocess
#import

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

def prepare_receptor(pdb_path : Path, pocket_id, center_coords, box_sizes) -> list[Path]:
    # Export receptor atoms
    center_x, center_y, center_z = center_coords
    size_x, size_y, size_z = box_sizes
    out_path = f'../data/docking_files/{pdb_path.stem}_{pocket_id}'
    command = [
        "python",
        str(mk_prepare_receptor),
        "-i", str(pdb_path),
        "-o", out_path,
        "-p",  # Generate PDBQT file
        "-v",  # Generate Vina config
        "--box_center", str(center_x), str(center_y), str(center_z),
        "--box_size", str(size_x), str(size_y), str(size_z)
    ]

    subprocess.run(command, check=True, shell=True)

    return (Path(f'{out_path}.pdbqt'), Path(f'{out_path}.box.txt'))

def parse_msa_target_indices(path):
    with open(path, 'r') as f:
        lines = f.readlines()

    res = []
    for line in lines:
        left, right = line.strip().split('-')
        res.extend(list(range(int(left), int(right))))

    return res

def find_best_pockets(struct : Structure, pockets : pd.DataFrame, msa : MultipleSeqAlignment, id_to_msa_index : dict, target_indices, tol=20, mode = 'close_all'):

    if mode == 'best':
        return pockets.sort_values('score').iloc[0]
    
    # List containing a mapping from MSA indices to last preceding sequence index
    ungapped_index = []
    idx = 0
    seq_row = id_to_msa_index[struct.id]
    for i in range(msa.get_alignment_length()):
        ungapped_index.append(idx)
        if msa[seq_row, i] != '-':
            idx += 1

    # Target indices to sequence indices
    residue_indices = []
    for i in target_indices:
        if msa[seq_row, i] != '-':
            residue_indices.append(ungapped_index[i])

    # Calculate a centroid from the target residues
    residues = list(struct.get_residues())
    centroid = np.average([residues[i].center_of_mass() for i in residue_indices], axis=0)

    # Mask pockets according to their distance from the centroid of selected residues 
    mask = []
    for _, pocket in pockets.iterrows():
        min_dist = np.inf
        for residue in pocket['residue_ids'].split(' '):
            idx = int(residue[2:]) # residue ids have a pattern of A_123
            # Iterate through atoms of the residue to find the furthest one
            for atom in residues[idx].get_atoms():
                c_dist = l2_norm(atom.get_coord() - centroid)
                if c_dist < min_dist:
                    min_dist = c_dist
        #center_x, center_y, center_z = pocket['center_x'], pocket['center_y'], pocket['center_z']
        print(min_dist)
        mask.append(min_dist < tol)
    
    # All pockets are far, return the best one
    if np.sum(mask) == 0:
        return pockets.sort_values('score').iloc[0]

    # Return all close pockets
    if mode == 'close_all':
        return pockets[mask]
    
    # Return the pocket from the close pockets that has the highest score
    else:
        return pockets[mask].sort_values('score').iloc[0]
    
def create_msa_index_table(msa : MultipleSeqAlignment):
    res = {}
    for i , prot in enumerate(msa):
        res[prot.id] = i

    return res

def calculate_box_size(residues, center, pocket):
    max_dist = 0
    # Find the furthest atom from the center for correct box size
    for residue in pocket['residue_ids'].split(' '):
        idx = int(residue[2:]) # residue ids have a pattern of A_123
        max_dist_atom = 0
        # Iterate through atoms of the residue to find the furthest one
        for atom in residues[idx].get_atoms():
            c_dist = l2_norm(atom.get_coord() - center)
            if c_dist > max_dist_atom:
                max_dist_atom = c_dist

        if max_dist_atom > max_dist:
            max_dist = max_dist_atom

    return max_dist + 2 # 2A padding

def protonate_pdb(pdb_path : Path, ph=7):
    fixer = PDBFixer(str(pdb_path))
    fixer.addMissingHydrogens(ph)
    out_path = f'../data/temp/{pdb_path.stem}_H.pdb'
    PDBFile.writeFile(fixer.topology, fixer.positions, open(out_path, 'w'))
    return out_path

def prepare_receptors(args) -> list[tuple[Path, Path]]:
    pdbs = list(Path(args.pdbs_path).rglob('*.pdb'))
    parser = PDBParser()
    msa = AlignIO.read(args.msa_path, format='fasta')
    # Mapping of unproit ids to rows in the MSA
    id_to_msa_index = create_msa_index_table(msa)

    # Prepare MSA index residues
    target_indices = parse_msa_target_indices(args.target_idxs_path)

    out = []
    for pdb in pdbs:
        name = pdb.stem
        try:
            pockets = pd.read_csv(f'{args.pocket_preds_path}/{name}.pdb_predictions.csv', skipinitialspace=True)
        except FileNotFoundError:
            raise FileExistsError(f'Predictions for {name} not found. Use run_p2rank.py to generate the pocket predictions first.')

        h_path = protonate_pdb(pdb, args.ph)
        molecule = parser.get_structure(pdb.stem, h_path)

        # Determine the pocket for docking        
        best = find_best_pockets(molecule, pockets, msa, id_to_msa_index, target_indices, tol=args.tol, mode=args.pocket_selection_mode)
        residues = list(molecule.get_residues())
        if best.shape[0] > 0:
            for _, pocket in best.iterrows():
                center = pocket['center_x'], pocket['center_y'], pocket['center_z']
                center_np = np.asarray(center)
                box_size = calculate_box_size(residues, center_np, pocket)
                out.append(prepare_receptor(pdb, f'p{pocket["rank"]}', center, (box_size, box_size, box_size)))
        break
    return out

def dock_ligands(receptor_info : list[tuple[Path, Path]]):
    if not os.path.exists('../data/docking_files'):
        os.makedirs('../data/docking_files')

    if not os.path.exists('../output'):
        os.mkdir('../output')
    
    if sys.platform == 'win32':
        ligands_path_cmd = ' --batch '.join(list(glob('..\data\prepared_ligands\*.pdbqt')))
        
    else:
        ligands_path_cmd =  Path('../data/prepared_ligands/*.pdbqt')

    for receptor, config in receptor_info:
        out_path = f'../output/{receptor.stem}'
        if not os.path.exists(out_path):
            os.makedirs(out_path, exist_ok=True)

        print(f'Docking ligands to {receptor.stem}...')
        command = ' '.join([
            str(vina),
            '--receptor', str(receptor),
            '--batch', ligands_path_cmd,
            '--config', str(config),
            f'--exhaustiveness={args.exhaustiveness}',
            '--cpu', str(os.cpu_count()),
            '--dir', out_path
        ])
        subprocess.run(command, shell=True, check=True)
        print(f'Output saved to {out_path}')

def run_docking(args):
    # prepare_ligands(args)
    receptor_info = prepare_receptors(args)
    dock_ligands(receptor_info)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pocket_preds_path', default='../data/p2rank_output', help='P2rank pocket predictions output path')
    parser.add_argument('--pdbs_path', default='../data/pdbs', help='Path to the directory where .pdb files are stored')
    parser.add_argument('--ligands_path', default='../data/ligands.csv', help='Path to the ligands csv file. Should contain fields "smiles" and "name".')
    parser.add_argument('--ph', default=7, type=int, help='pH for ligand preparation')
    parser.add_argument('--msa_path', default='../data/aligned_sequences.fasta', help='MSA for the receptors, used for selecting pockets outside the membrane.')
    parser.add_argument('--target_idxs_path', default='../data/msa_index_ranges.txt', help='File containing indices to the MSA for membrane pocket selection')
    parser.add_argument('--tol', default=20, type=int, help='Maximum pocket center distance from the centroid of residues matched to the MSA indices for membrane pocket selection. Pockets further away are not considered when determining the best pocket.')
    parser.add_argument('--skip_tautomer', action='store_true', help='Skip tautomers in ligand preparation')
    parser.add_argument('--skip_acidbase', action='store_true', help='Skip acid/base conjugates in ligand preparation')
    # parser.add_argument('--box_size', type=int, default=20, help='Box size in Å for docking.')
    parser.add_argument('-e', '--exhaustiveness', type=int, default=32, help='Search exhaustiveness for Vina')
    parser.add_argument('--pocket_selection_mode', default='close_all', choices=['best', 'close_best', 'close_all'], help='''
                        "best" mode simply takes the highest scoring pocket from the P2rank prediction.
                        "close_best" mode only considers the best pocket from ones that are close to the external part of the protein, determined via the indices from target_idxs_path
                        "close_all" same as above, but docks to all close pockets
                        ''')
    args = parser.parse_args()
    run_docking(args)