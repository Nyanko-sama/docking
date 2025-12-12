import argparse
import subprocess
import platform
import shutil
import sys
from pathlib import Path
from utils import locate_file

def find_clustalo():
    """Find clustalo executable, handling both Windows and Linux."""
    if platform.system() == 'Windows':
        try:
            clustal = locate_file(Path.cwd().parent, 'clustalo.exe', 'clustalo.exe')
            return str(clustal)
        except FileNotFoundError:
            raise FileNotFoundError(
                "clustalo.exe not found. Please install ClustalO or place clustalo.exe "
                "in the parent directory of the source folder."
            )
    else:
        # On Linux/Unix, try to find clustalo in PATH first
        clustal_path = shutil.which('clustalo')
        if clustal_path:
            print(f"Using clustalo from PATH: {clustal_path}")
            return clustal_path
        
        # Fall back to searching in parent directory
        try:
            clustal = locate_file(Path.cwd().parent, 'clustalo', 'clustalo')
            return str(clustal)
        except FileNotFoundError:
            # Not found - provide helpful error message
            raise FileNotFoundError(
                "clustalo not found in PATH or project directory.\n"
                "Please install it using one of the following methods:\n"
                "  1. conda install -c bioconda clustalo\n"
                "  2. Or add it to your environment.yaml file under dependencies:\n"
                "     - clustalo\n"
                "     And add 'bioconda' to channels:\n"
                "     - bioconda\n"
                "  3. Or install via your system package manager (e.g., apt-get install clustalo)"
            )

# Initialize clustal_cmd at module level
try:
    clustal_cmd = find_clustalo()
except FileNotFoundError as e:
    clustal_cmd = None
    print(f"WARNING: {e}", file=sys.stderr)

def align_sequences(seqs_path, output_path):
    if clustal_cmd is None:
        raise RuntimeError(
            "clustalo is not available. Please install it first.\n"
            "Run: conda install -c bioconda clustalo"
        )
    subprocess.run(f'{clustal_cmd} -i {seqs_path} -o {output_path} --auto', shell=True, check=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_path', type=str, default='../data/pdb_sequences.fasta', help='Input FASTA file containing sequences to be aligned')
    parser.add_argument('-o', '--output_path', type=str, default='../data/aligned_sequences.fasta', help='Output folder path')


    args = parser.parse_args()
    align_sequences(args.input_path, args.output_path)