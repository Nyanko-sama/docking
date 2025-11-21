import argparse
import pandas as pd
from glob import glob
from Bio.PDB import PDBParser
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def extract_seqs(pdb_folder, output_folder):
    parser = PDBParser()

    residue_translation = {
        "ALA": "A",  # Alanine
        "ARG": "R",  # Arginine
        "ASN": "N",  # Asparagine
        "ASP": "D",  # Aspartic Acid
        "CYS": "C",  # Cysteine
        "GLN": "Q",  # Glutamine
        "GLU": "E",  # Glutamic Acid
        "GLY": "G",  # Glycine
        "HIS": "H",  # Histidine
        "ILE": "I",  # Isoleucine
        "LEU": "L",  # Leucine
        "LYS": "K",  # Lysine
        "MET": "M",  # Methionine
        "PHE": "F",  # Phenylalanine
        "PRO": "P",  # Proline
        "SER": "S",  # Serine
        "THR": "T",  # Threonine
        "TRP": "W",  # Tryptophan
        "TYR": "Y",  # Tyrosine
        "VAL": "V",  # Valine
        # Note: X is often used for unknown or non-standard residues
    }
    records = []
    fasta_recs = []
    for pdb_path in glob(f'{pdb_folder}/*.pdb', recursive=True):
        path = Path(pdb_path)
        struct = parser.get_structure(path.stem, pdb_path)
        tuple_rec = (path.stem, ''.join([residue_translation[r.get_resname()] for r in struct.get_residues()]))
        fasta_rec = SeqRecord(Seq(tuple_rec[1]), tuple_rec[0], description='')
        records.append(tuple_rec)
        fasta_recs.append(fasta_rec)

    pd.DataFrame.from_records(records, columns=['uniprot_id', 'sequence']).to_csv(f'{output_folder}/pdb_sequences.csv', index=False)
    SeqIO.write(fasta_recs, f'{output_folder}/pdb_sequences.fasta', format='fasta')
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_folder', '-i', type=str, default='../data/pdbs', help='Extracts sequences from all PDB files in this folder')
    parser.add_argument('--output_folder', '-o', type=str, default='../data', help='Folder where the output files will be written')

    args = parser.parse_args()
    extract_seqs(args.input_folder, args.output_folder)