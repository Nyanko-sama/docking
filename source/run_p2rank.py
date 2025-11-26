import os
import argparse
import glob
import subprocess

def create_ds_file(pdbs_path):
    paths = glob.glob(f'{pdbs_path}/*.pdb')
    with open('../data/pdbs.ds', 'w') as f:
        f.write('\n'.join(paths))

def run_p2rank(p2rank_path, ds_path, outputd_dir):
    if not os.path.exists(outputd_dir):
        os.mkdir(outputd_dir)

    subprocess.run([
        os.path.join(p2rank_path, 'prank'),
        'predict', ds_path, 
        '-o', outputd_dir
    ])
    #os.system(f'{p2rank_path}\\prank predict {ds_path} -o {outputd_dir}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--p2rank_path', default='..\\p2rank', help='p2rank root directory')
    parser.add_argument('-o', '--output_path', default='../data/p2rank_output', help='Output directory path')
    parser.add_argument('-i', '--pdbs_path', default='../data/pdbs', help='Path to the directory where .pdb files are stored')

    args = parser.parse_args()
    create_ds_file(args.pdbs_path)
    run_p2rank(args.p2rank_path, '../data/pdbs.ds', args.output_path)
    