import subprocess
import pandas as pd
import numpy as np
import argparse
import os

from utils import locate_file, l2_norm, get_path_root
from pathlib import Path
from Bio.PDB import Structure, PDBParser
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from datetime import datetime

mk_prepare_receptor = locate_file(from_path = get_path_root(), query_path = "mk_prepare_receptor.py", query_name = "mk_prepare_receptor.py")
#reduce = locate_file(from_path = Path(str(get_path_root().parent) + '/lib'), query_path = "reduce.py", query_name = "reduce.py")

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

    subprocess.run(command, check=True)

    return (Path(f'{out_path}.pdbqt'), Path(f'{out_path}.box.txt'))

def parse_msa_target_indices(path):
    with open(path, 'r') as f:
        lines = f.readlines()

    res = []
    for line in lines:
        left, right = line.strip().split('-')
        res.extend(list(range(int(left), int(right))))

    return res

def find_best_pockets(struct : Structure, pockets : pd.DataFrame, msa : MultipleSeqAlignment, id_to_msa_index : dict, target_indices, tol=20, mode = 'close_all',
                      verbose=1, include_best=True):

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

    if len(residue_indices) == 0:
        # No matched residues

        # Log the protein
        if not os.path.exists(f'../temp/no_match_prots.txt'):
            os.makedirs('../temp', exist_ok=True)
        
        with open(f'../temp/no_match_prots.txt', 'r') as f:
            lines = f.readlines()
            lines.append(struct.id)

        with open(f'../temp/no_match_prots.txt', 'w') as f:
            f.writelines(lines)
            
        if mode == 'best':
            return pockets.sort_values('score').iloc[0]

    # Calculate a centroid from the target residues
    residues = list(struct.get_residues())
    centroid = np.average([residues[i].center_of_mass() for i in residue_indices], axis=0)

    # Mask pockets according to their distance from the centroid of selected residues 
    mask = []
    pockets = pockets.sort_values('score')
    if verbose > 0:
        print('Pocket distance from selected residue centroid, order as in the predictions .csv:')
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
        if verbose > 0:
            print(min_dist)
        mask.append(min_dist < tol)
    
    if include_best:
        # Always include the best pocket
        mask[0] = True
    if verbose > 0:
        print(mask)
        
    # All pockets are far, return the best one
    if np.sum(mask) == 0:
        return pockets.iloc[0]

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
    out_path = f'../data/temp/{pdb_path.stem}_H.pdb'

    #subprocess.run([
    #    reduce, '-FLIP', pdb_path, '>', out_path
    #])
    fixer = PDBFixer(str(pdb_path))
    fixer.findNonstandardResidues()
    print(fixer.nonstandardResidues)
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.addMissingHydrogens(ph)

    PDBFile.writeFile(fixer.topology, fixer.positions, open(out_path, 'w'))
    return Path(out_path)

def prepare_receptors(args) -> list[tuple[Path, Path]]:
    pdbs = list(Path(args.pdbs_path).rglob('*.pdb'))
    parser = PDBParser()
    msa = AlignIO.read(args.msa_path, format='fasta')
    # Mapping of uniprot ids to rows in the MSA
    id_to_msa_index = create_msa_index_table(msa)

    # Prepare MSA index residues
    target_indices = parse_msa_target_indices(args.target_idxs_path)

    if not os.path.exists('../data/docking_files/'):
        os.makedirs('../data/docking_files')

    out = []
    for pdb in pdbs:
        name = pdb.stem
        try:
            pockets = pd.read_csv(f'{args.pocket_preds_path}/{name}.pdb_predictions.csv', skipinitialspace=True)
        except FileNotFoundError:
            raise FileExistsError(f'Predictions for {name} not found. Use run_p2rank.py to generate the pocket predictions first.')

        h_path = protonate_pdb(pdb, args.ph)
        molecule = parser.get_structure(pdb.stem, pdb)

        # Determine the pocket for docking        
        best = find_best_pockets(molecule, pockets, msa, id_to_msa_index, target_indices, tol=args.tol, mode=args.pocket_selection_mode, include_best=args.include_best)
        residues = list(molecule.get_residues())
        if best.shape[0] > 0:
            for _, pocket in best.iterrows():
                center = pocket['center_x'], pocket['center_y'], pocket['center_z']
                center_np = np.asarray(center)
                box_size = calculate_box_size(residues, center_np, pocket)
                out.append(prepare_receptor(h_path, f'p{pocket["rank"]}', center, (box_size, box_size, box_size)))

    return out

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pocket_preds_path', default='../data/p2rank_output', help='P2rank pocket predictions output path')
    parser.add_argument('--pdbs_path', default='../data/pdbs', help='Path to the directory where .pdb files are stored')
    parser.add_argument('--ph', default=7, type=int, help='pH for protonation')
    parser.add_argument('--msa_path', default='../data/aligned_sequences.fasta', help='MSA for the receptors, used for selecting pockets outside the membrane.')
    parser.add_argument('--target_idxs_path', default='../data/msa_index_ranges.txt', help='File containing indices to the MSA for membrane pocket selection')
    parser.add_argument('--tol', default=20, type=int, help='Maximum pocket center distance from the centroid of residues matched to the MSA indices for membrane pocket selection. Pockets further away are not considered when determining the best pocket.')
    parser.add_argument('--pocket_selection_mode', default='close_all', choices=['best', 'close_best', 'close_all'], help='''
                        "best" mode simply takes the highest scoring pocket from the P2rank prediction.
                        "close_best" mode only considers the best pocket from ones that are close to the external part of the protein, determined via the indices from target_idxs_path
                        "close_all" same as above, but docks to all close pockets
                        ''')
    parser.add_argument('-v', '--verbose', default=1, type=int, help='Verbosity level')
    parser.add_argument('--include_best', default=True, type=bool, help='Always include the best pocket in the selected pockets. Relevant for the "close" pocket selection mode.')
    args = parser.parse_args()
    prepare_receptors(args)
