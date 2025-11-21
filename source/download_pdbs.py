import pandas as pd
import requests
import argparse
import os

def load_prots(prots_path):
    return pd.read_csv(prots_path, sep='\t')['uniprot_id']

def load_key(key_path):
    with open(key_path, 'r') as f:
        return f.read()
    
def download_pdbs(ids, key):
    for prot_id in ids:
        URL = f"https://alphafold.ebi.ac.uk/api/prediction/{prot_id}?"
        params = {'key' : key}
        if not os.path.exists('../data/pdbs'):
            os.mkdir('../data/pdbs/')

        r = requests.get(URL, params)
        if r.status_code == 200:
            data = r.json()[0]
            pdb_url = data['pdbUrl']
            
            r = requests.get(pdb_url)
            path = f'../data/pdbs/{prot_id}.pdb'
            with open(path, 'wb') as f:
                f.write(r.content)

            print(f'{prot_id} saved to {path}')

        else:
            print(f'Protein {prot_id} not found in the AlphaFold database, skipped. Received response {r.status_code}')
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--prots_path', default='../data/prots.csv', help='Path to a protein tsv-like file. Should contain a column name "uniprot_id", and the values should be Uniprot IDs of the proteins to be downloaded from AFDB.')
    parser.add_argument('--key_path', default='../data/key.txt', help='Location of the security key file for AFDB. Text file should contain ONLY the key.')
    args = parser.parse_args()
    download_pdbs(load_prots(args.prots_path), load_key(args.key_path))