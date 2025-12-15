import subprocess
import pandas as pd
import numpy as np
import argparse
import os
import sys

from utils import locate_file, l2_norm, get_path_root
from pathlib import Path
from Bio.PDB import Structure, PDBParser, PDBIO
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
        "-a",  # Allow bad residues (automatically remove problematic ones)
        "--box_center", str(center_x), str(center_y), str(center_z),
        "--box_size", str(size_x), str(size_y), str(size_z)
    ]
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        # Check if there were warnings about problematic residues
        if result.stderr:
            # Filter out expected warnings about template matching failures that are handled
            stderr_lines = result.stderr.split('\n')
            for line in stderr_lines:
                if 'Template matching failed' in line and 'Ignored due to allow_bad_res' in line:
                    print(f"Warning (handled): {line}")
                elif 'RuntimeError' in line or 'Error' in line:
                    # This is a real error, not just a handled warning
                    print(f"Error in mk_prepare_receptor: {line}")
        if result.stdout:
            print(result.stdout)

    except subprocess.CalledProcessError as e:
        print(f"\n{'='*60}")
        print(f"Error preparing receptor {pdb_path.stem}_{pocket_id}:")
        print(f"{'='*60}")
        
        # Check for specific error types and provide helpful messages
        error_output = ""
        if e.stderr:
            error_output = e.stderr
        if e.stdout:
            error_output += "\n" + e.stdout if error_output else e.stdout
        
        if "Expected" in error_output and "paddings" in error_output:
            print("Structural connectivity error detected:")
            print("  The PDB file has excess inter-residue bonds that cannot be resolved.")
            print("  This typically occurs when:")
            print("    - Residues have incorrect connectivity")
            print("    - Missing atoms cause bond mismatches")
            print("    - Chain breaks are not properly handled")
            print("    - The structure is missing the first part of the protein (N-terminal region)")
            print("      causing improper polymer connectivity")
            print("  Suggestion: Try using a different PDB structure or manually fix the connectivity.")
            print("  If the structure is missing N-terminal residues, consider using a complete structure.")
        elif "No template matched" in error_output and "residue_key" in error_output:
            print("Template matching error detected:")
            print("  A residue at the start of a chain could not be matched to a standard template.")
            print("  This often occurs when:")
            print("    - The structure is missing the first part of the protein")
            print("    - The N-terminal residue is missing required atoms (N, CA, C, O)")
            print("    - The first residue is a non-standard or modified residue")
            print("  PDBFixer should add missing atoms, but if the structure starts mid-sequence,")
            print("  template matching may still fail. Consider using a complete structure.")
        elif "Template matching failed" in error_output:
            print("Template matching error (this should be handled by -a flag):")
            print("  Some residues could not be matched to standard templates.")
        
        print(f"\nCommand: {' '.join(command)}")
        print(f"Return code: {e.returncode}")
        if error_output:
            print(f"\nFull error output:\n{error_output}")
        
        # Ensure problem_pdbs.txt directory exists
        os.makedirs('../data', exist_ok=True)
        with open('../data/problem_pdbs.txt', 'a') as f:
            f.write(f"{pdb_path} (pocket {pocket_id})\n")
        print(f"{'='*60}\n")
        
        # Return None to indicate failure, caller should handle this
        return None

    return (Path(f'{out_path}.pdbqt'), Path(f'{out_path}.box.txt'))

def parse_msa_target_indices(path):
    """Parse MSA target indices from a file. Returns empty list if file doesn't exist."""
    if not os.path.exists(path):
        return []
    
    with open(path, 'r') as f:
        lines = f.readlines()

    res = []
    for line in lines:
        line = line.strip()
        if not line:  # Skip empty lines
            continue
        left, right = line.split('-')
        res.extend(list(range(int(left), int(right))))

    return res

def find_best_pockets(struct : Structure, pockets : pd.DataFrame, msa : MultipleSeqAlignment, id_to_msa_index : dict, target_indices, tol=20, mode = 'close_all',
                      verbose=1, include_best=True):
    if verbose > 0:
        print(struct.id)

    pockets = pockets.sort_values('score', ascending=False)
    if mode == 'best':
        return pockets.iloc[0]
    
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
            with open('../temp/no_match_prots.txt', 'w') as f:
                f.write(struct.id)

        else:
            with open(f'../temp/no_match_prots.txt', 'r') as f:
                lines = f.readlines()
                lines.append(struct.id)

            with open(f'../temp/no_match_prots.txt', 'w') as f:
                f.writelines(lines)

        return pockets.iloc[0]

    # Calculate a centroid from the target residues
    residues = list(struct.get_residues())
    centroid = np.average([residues[i].center_of_mass() for i in residue_indices], axis=0)

    # Mask pockets according to their distance from the centroid of selected residues 
    mask = []

    if verbose > 0:
        print('Pocket distance from selected residue centroid, order as in the predictions .csv:')
    for _, pocket in pockets.iterrows():
        min_dist = np.inf
        for residue in pocket['residue_ids'].split(' '):
            idx = int(residue[2:]) # residue ids have a pattern of A_123
            # Iterate through atoms of the residue to find the furthest one
            try:
    	        for atom in residues[idx].get_atoms():
                    c_dist = l2_norm(atom.get_coord() - centroid)
                    if c_dist < min_dist:
                        min_dist = c_dist
            except IndexError as e:
                print(e)
                print(struct.id)
                print(idx)
                print(len(list(struct.get_residues())))
                
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
        return pockets[mask].iloc[0]
    
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

def check_missing_nterminal(pdb_path : Path) -> dict:
    """Check if PDB structure has missing N-terminal regions.
    
    Returns a dict with information about chains that may have missing N-terminal atoms.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('temp', str(pdb_path))
    
    required_backbone_atoms = {'N', 'CA', 'C', 'O'}
    issues = {}
    
    for model in structure:
        for chain in model.get_chains():
            residues = list(chain.get_residues())
            if len(residues) == 0:
                continue
                
            first_res = residues[0]
            resname = first_res.get_resname()
            resnum = first_res.get_id()[1]
            first_res_atoms = {atom.get_id() for atom in first_res.get_atoms()}
            missing_atoms = required_backbone_atoms - first_res_atoms
            
            if missing_atoms:
                chain_key = f"{model.id}:{chain.id}"
                issues[chain_key] = {
                    'residue_name': resname,
                    'residue_number': resnum,
                    'missing_atoms': missing_atoms,
                    'has_all_backbone': False
                }
            else:
                chain_key = f"{model.id}:{chain.id}"
                issues[chain_key] = {
                    'residue_name': resname,
                    'residue_number': resnum,
                    'missing_atoms': set(),
                    'has_all_backbone': True
                }
    
    return issues

def clean_pdb_with_biopython(pdb_path : Path) -> Path:
    """Pre-clean PDB file using BioPython to remove problematic chains/residues.
    
    Also detects and reports missing N-terminal regions that can cause template matching failures.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('temp', str(pdb_path))
    
    # Required backbone atoms for a proper N-terminus
    required_backbone_atoms = {'N', 'CA', 'C', 'O'}
    
    # Remove chains that have very few residues (likely ligands or problematic chains)
    # Also remove chains that start with non-standard residues at position 1
    chains_to_remove = []
    chains_with_missing_nterm = []
    
    for model in structure:
        for chain in model.get_chains():
            residues = list(chain.get_residues())
            # Remove chains with very few residues (likely not protein chains)
            if len(residues) < 3:
                chains_to_remove.append((model.id, chain.id))
                continue
            
            # Check if first residue is problematic (non-standard or missing key atoms)
            if len(residues) > 0:
                first_res = residues[0]
                resname = first_res.get_resname()
                standard_residues = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 
                                    'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 
                                    'THR', 'TRP', 'TYR', 'VAL'}
                
                # Check if first residue has all required backbone atoms
                first_res_atoms = {atom.get_id() for atom in first_res.get_atoms()}
                missing_atoms = required_backbone_atoms - first_res_atoms
                
                if resname in standard_residues:
                    # Standard residue but missing backbone atoms - likely missing N-terminal region
                    if missing_atoms:
                        chains_with_missing_nterm.append((model.id, chain.id, resname, missing_atoms))
                        print(f"Warning: Chain {chain.id} in {pdb_path.stem} appears to have missing N-terminal atoms "
                              f"(residue {resname} at position 1 missing: {missing_atoms}). "
                              f"PDBFixer should add these, but template matching may still fail.")
                elif resname not in {'NME', 'ACE'} and not resname.startswith('H'):
                    # Non-standard residue that's not a cap - likely problematic
                    chains_to_remove.append((model.id, chain.id))
                    continue
    
    # Remove identified problematic chains
    for model_id, chain_id in chains_to_remove:
        model = structure[model_id]
        model.detach_child(chain_id)
        print(f"Removed chain {chain_id} from model {model_id} in {pdb_path.stem} (likely problematic)")
    
    # Write cleaned structure
    os.makedirs('../data/temp', exist_ok=True)
    cleaned_path = f'../data/temp/{pdb_path.stem}_cleaned.pdb'
    io = PDBIO()
    io.set_structure(structure)
    io.save(cleaned_path)
    return Path(cleaned_path)

def protonate_pdb(pdb_path : Path, ph=7):
    out_path = f'../data/temp/{pdb_path.stem}_H.pdb'

    # Check for missing N-terminal regions before cleaning
    nterm_issues = check_missing_nterminal(pdb_path)
    for chain_key, info in nterm_issues.items():
        if not info['has_all_backbone']:
            print(f"Warning: {pdb_path.stem} chain {chain_key} first residue ({info['residue_name']} {info['residue_number']}) "
                  f"is missing backbone atoms: {info['missing_atoms']}. "
                  f"This suggests a missing N-terminal region. PDBFixer will attempt to add missing atoms.")

    # First try to clean with BioPython to remove obviously problematic chains
    try:
        cleaned_path = clean_pdb_with_biopython(pdb_path)
        input_path = cleaned_path
    except Exception as e:
        print(f"Warning: BioPython cleaning failed for {pdb_path}: {e}, using original file")
        input_path = pdb_path

    fixer = PDBFixer(str(input_path))
    fixer.findNonstandardResidues()
    print(fixer.nonstandardResidues)
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    # addMissingAtoms should add missing N-terminal atoms if the structure is missing the first part
    fixer.addMissingAtoms(seed=42)
    # subprocess.run([
    #    reduce, '-FLIP', pdb_path, '>', out_path
    # ])
    fixer.addMissingHydrogens(ph)
    
    # Note: If the structure is missing the first part of the protein, PDBFixer's addMissingAtoms
    # should add the missing N-terminal atoms. However, if the structure starts mid-sequence
    # (e.g., residue 50 instead of residue 1), template matching in meeko may still fail because
    # it expects proper N-terminal capping or connectivity.
    
    # Additional cleaning: remove chains with no residues or problematic chains
    # This helps with template matching failures
    try:
        PDBFile.writeFile(fixer.topology, fixer.positions, open(out_path, 'w'))
    except Exception as e:
        print(f"Warning: Error writing fixed PDB for {pdb_path}: {e}")
        # Try to write without some problematic chains if possible
        raise
    
    return Path(out_path)

def prepare_receptors(args) -> list[tuple[Path, Path]]:
    pdbs = list(Path(args.pdbs_path).rglob('*.pdb'))
    parser = PDBParser()
    msa = AlignIO.read(args.msa_path, format='fasta')
    # Mapping of uniprot ids to rows in the MSA
    id_to_msa_index = create_msa_index_table(msa)

    # Prepare MSA index residues
    target_indices = parse_msa_target_indices(args.target_idxs_path)
    
    # Warn if target_indices are needed but not provided
    if args.pocket_selection_mode in ['close_best', 'close_all'] and len(target_indices) == 0:
        if not os.path.exists(args.target_idxs_path):
            print(f"WARNING: target_idxs_path '{args.target_idxs_path}' not found.", file=sys.stderr)
            print(f"WARNING: MSA index ranges are required for '{args.pocket_selection_mode}' mode.", file=sys.stderr)
            print(f"WARNING: Switching to 'best' mode (using highest scoring pocket only).", file=sys.stderr)
            args.pocket_selection_mode = 'best'
        else:
            print(f"WARNING: target_idxs_path '{args.target_idxs_path}' exists but is empty.", file=sys.stderr)
            print(f"WARNING: Switching to 'best' mode (using highest scoring pocket only).", file=sys.stderr)
            args.pocket_selection_mode = 'best'

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
        residues = {r.get_id()[1] : r for r in list(molecule.get_residues())}
        if len(best.shape) > 1:
            for _, pocket in best.iterrows():
                center = pocket['center_x'], pocket['center_y'], pocket['center_z']
                center_np = np.asarray(center)
                box_size = calculate_box_size(residues, center_np, pocket)
                result = prepare_receptor(h_path, f'p{pocket["rank"]}', center, (box_size, box_size, box_size))
                if result is not None:
                    out.append(result)
                else:
                    print(f"Warning: Failed to prepare receptor {pdb.stem} for pocket {pocket['rank']}, skipping...")
        else:
            center = best['center_x'], best['center_y'], best['center_z']
            center_np = np.asarray(center)
            box_size = calculate_box_size(residues, center_np, best)
            result = prepare_receptor(h_path, f'p{best["rank"]}', center, (box_size, box_size, box_size))
            if result is not None:
                out.append(result)
            else:
                print(f"Warning: Failed to prepare receptor {pdb.stem} for pocket {best['rank']}, skipping...")

    return out

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pocket_preds_path', default='../data/p2rank_output', help='P2rank pocket predictions output path')
    parser.add_argument('--pdbs_path', default='../data/pdbs', help='Path to the directory where .pdb files are stored')
    parser.add_argument('--ph', default=7, type=int, help='pH for protonation')
    parser.add_argument('--msa_path', default='../data/aligned_sequences.fasta', help='MSA for the receptors, used for selecting pockets outside the membrane.')
    parser.add_argument('--target_idxs_path', default='../data/msa_index_ranges.txt', help='File containing indices to the MSA for membrane pocket selection')
    parser.add_argument('--tol', default=20, type=int, help='Maximum pocket center distance from the centroid of residues matched to the MSA indices for membrane pocket selection. Pockets further away are not considered when determining the best pocket.')
    parser.add_argument('--pocket_selection_mode', default='best', choices=['best', 'close_best', 'close_all'], help='''
                        "best" mode simply takes the highest scoring pocket from the P2rank prediction.
                        "close_best" mode only considers the best pocket from ones that are close to the external part of the protein, determined via the indices from target_idxs_path
                        "close_all" same as above, but docks to all close pockets
                        ''')
    parser.add_argument('-v', '--verbose', default=1, type=int, help='Verbosity level')
    parser.add_argument('--include_best', default=True, type=bool, help='Always include the best pocket in the selected pockets. Relevant for the "close" pocket selection mode.')
    args = parser.parse_args()
    prepare_receptors(args)
