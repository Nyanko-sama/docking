import argparse
import subprocess
from pathlib import Path
from utils import locate_file

clustal = locate_file(Path.cwd().parent, 'clustalo.exe', 'clustalo.exe')

def align_sequences(seqs_path, output_path):
    subprocess.run(f'{clustal} -i {seqs_path} -o {output_path} --auto', shell=True, check=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_path', type=str, default='../data/pdb_sequences.fasta', help='Input FASTA file containing sequences to be aligned')
    parser.add_argument('-o', '--output_path', type=str, default='../data/aligned_sequences.fasta', help='Output folder path')


    args = parser.parse_args()
    align_sequences(args.input_path, args.output_path)