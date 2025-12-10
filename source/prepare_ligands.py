import argparse
import re
import sys
import os
import subprocess
import pandas as pd
from pathlib import Path
from utils import locate_file, get_path_root

mk_prepare_ligand = locate_file(from_path=get_path_root(), query_path="mk_prepare_ligand.py", query_name="mk_prepare_ligand.py")
scrub = locate_file(from_path=Path.cwd().parent, query_path="scrub.py", query_name="scrub.py")


def prepare_ligand(ligand_smiles, ph=7, skip_tautomer=False, skip_acidbase=False,
                   ligand_name='test_ligand', timeout=None) -> Path:
    """
    Prepare a single ligand with optional timeout.
    If timeout is exceeded, subprocess.TimeoutExpired will be raised.
    """
    # Adapted from:
    # https://colab.research.google.com/drive/1cHSl78lBPUc_J1IZxLgN4GwD_ADmohVU?usp=sharing#scrollTo=qBQN6WzwvkGB
    args = []
    if skip_tautomer:
        args.append("--skip_tautomers")

    if skip_acidbase:
        args.append("--skip_acidbase")

    ligand_name = re.sub(r'-', '_', ligand_name)
    ligandSDF = f"../data/temp/{ligand_name}_scrubbed.sdf"
    output_path = f'../data/prepared_ligands/{ligand_name}.pdbqt'

    # Scrub the molecule
    scrub_cmd = ['python', str(scrub), ligand_smiles, '-o', ligandSDF, '--ph', str(ph)] + args
    subprocess.run(scrub_cmd, check=True, timeout=timeout)

    # Run meeko mk_prepare_ligand
    mk_cmd = [
        'python', str(mk_prepare_ligand),
        '-i', ligandSDF,
        '-o', output_path
    ]
    subprocess.run(mk_cmd, check=True, timeout=timeout)

    return Path(output_path)


def standardize_name(ligand_name, max_length=200):
    """
    Sanitize ligand name for use as filename.
    Removes/replaces problematic characters and truncates if too long.
    """
    if pd.isna(ligand_name):
        return 'unknown_ligand'
    
    ligand_name = str(ligand_name)
    
    # Remove quotes
    ligand_name = re.sub(r'[\'"]', '', ligand_name)
    
    # Replace problematic characters with underscores
    # This includes: parentheses, brackets, semicolons, slashes, colons, etc.
    ligand_name = re.sub(r'[()\[\];/\\:<>|?*]', '_', ligand_name)
    
    # Replace spaces, dashes, commas, and other separators with underscores
    ligand_name = re.sub(r'[-\s,]+', '_', ligand_name)
    
    # Remove multiple consecutive underscores
    ligand_name = re.sub(r'_+', '_', ligand_name)
    
    # Remove leading/trailing underscores
    ligand_name = ligand_name.strip('_')
    
    # If empty after sanitization, use a default name
    if not ligand_name:
        ligand_name = 'unknown_ligand'
    
    # Truncate if too long, but try to preserve meaningful part
    if len(ligand_name) > max_length:
        # Try to truncate at a word boundary (underscore)
        truncated = ligand_name[:max_length]
        last_underscore = truncated.rfind('_')
        if last_underscore > max_length * 0.7:  # If we can keep at least 70% of the name
            ligand_name = truncated[:last_underscore]
        else:
            ligand_name = truncated
    
    return ligand_name


def prepare_ligands(args):
    ligands = pd.read_csv(args.ligands_path)
    ligands = ligands.reindex(['smiles', 'name'], axis=1)
    ligands['name'] = ligands['name'].apply(standardize_name)

    if not os.path.exists('../data/prepared_ligands'):
        os.makedirs('../data/prepared_ligands')

    if not os.path.exists('../data/temp'):
        os.makedirs('../data/temp')

    skipped_ligands = []

    # Process ligands one by one, with timeout and error handling
    for idx, row in ligands.iterrows():
        smiles = row['smiles']
        name = row['name']

        print(f"\n=== Preparing ligand {idx + 1}/{len(ligands)}: {name} ===", flush=True)

        try:
            prepare_ligand(
                ligand_smiles=smiles,
                ph=args.ph,
                skip_tautomer=args.skip_tautomers,
                skip_acidbase=args.skip_acidbase,
                ligand_name=name,
                timeout=args.timeout
            )
        except subprocess.TimeoutExpired:
            print(f"TIMEOUT: Skipping ligand '{name}' (took longer than {args.timeout} seconds).", file=sys.stderr, flush=True)
            skipped_ligands.append(name)
        except subprocess.CalledProcessError as e:
            print(f"ERROR: Skipping ligand '{name}' due to subprocess error: {e}", file=sys.stderr, flush=True)
            skipped_ligands.append(name)
        except Exception as e:
            print(f"UNEXPECTED ERROR for ligand '{name}': {e}", file=sys.stderr, flush=True)
            skipped_ligands.append(name)

    # Summary of skipped ligands
    if skipped_ligands:
        print("\n======================================")
        print("The following ligands were SKIPPED:")
        for name in skipped_ligands:
            print(f" - {name}")
        print("======================================")
    else:
        print("\nAll ligands were processed successfully.")

    return ligands['name']


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--ligands_path',
        default='../data/ligands.csv',
        help='Path to the ligands csv file. Should contain fields "smiles" and "name".'
    )
    parser.add_argument(
        '--ph', default=7, type=int,
        help='pH for ligand preparation'
    )
    parser.add_argument(
        '--tol', default=20, type=int,
        help=('Maximum pocket center distance from the centroid of residues matched to the MSA indices '
              'for membrane pocket selection. Pockets further away are not considered when determining '
              'the best pocket.')
    )
    parser.add_argument(
        '--skip_acidbase', action='store_true',
        help='Skip acid/base conjugates in ligand preparation'
    )
    parser.add_argument(
        '--skip_tautomers', action='store_true',
        help='Skip tautomers in ligand preparation.'
    )
    parser.add_argument(
        '--timeout', type=int, default=60,
        help='Per-ligand timeout in seconds for external tools (scrub, mk_prepare_ligand).'
    )

    args = parser.parse_args()
    prepare_ligands(args)
