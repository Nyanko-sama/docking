# PrankDock - Automated pocket prediction and docking tool



## Setup

### 1. Environmnent setup
- Download and install the Anaconda package manager, and create a new virtual environment for prankdock using the `.yml` file.    
```
conda env create -f environment.yml
```
- Activate the environment using
```
conda activate dock
```
- Download the [AutodockVina 4](https://github.com/ccsb-scripps/AutoDock-Vina/releases) binary and place it somewhere in the prankdock directory
-  Install [P2rank](https://github.com/rdk/p2rank) and [molscrub](https://github.com/forlilab/molscrub). Defaults behaviour of scripts in `source/` expects these tools to be located in the root directory of prankdock. If you install them somewhere else, do not forget to provide the correct path in the command line arguments of `run_p2rank.py` and `run_docking.py` 

### 2. PDB Downloads
If you wish to download PDBs via `download_pdbs.py`, you need to obtain an API key for AlphaFold Database. Create a `.txt` file containing the key at this path: `data/key.txt`.
You will need to install the 


## Usage
- NOTE: The scripts mentioned in this section expect to be run from the `source/` directory, so first run
```
cd source
```

- If you're unsure about anything, you can run any script with the `-h` or `--help` option to see argument explanations.
1. Download PDBs from AFDB via `source/download_pdbs.py`. 
    - This script requires a `.csv`containing the Uniprot IDs of requested receptors, inside a column named `uniprot_id`.
    - Default path is `data/prots.csv`
    - If you already have PDBs downloaded, place them into `data/pdbs`.
2. Predict pockets with P2rank using `source/run_p2rank.py`.
    -  By default, the script gathers the paths of all `.pdb` files in `data/pdbs`, runs P2rank, and deposits the results into `data/p2rank_output`
4. Prepare ligands using `source/prepare_ligands.py`.
    - This script requires a `.csv` file containing the following fields: `name` and `smiles`. `name` is the name of the ligand, and `smiles` is its smiles string.
    - Default location of this file is `data/ligands.csv`
    - If you wish to skip isomers, use the `--skip_tautomer` and `--skip_acidbase` arguments

5. Prepare receptors using `source/prepare_receptors.py`

    - This script handles receptor preparation, including addition of hydrogen atoms and docking pocket selection.
        - We estimated the residues close to the external part of the membrane by computing the MSA of the proteins from `prots.csv` via `source/extract_sequences.py` and `align_sequences.py`. Then, we inspected the alignment and chose a few subsets of global MSA indices that roughly corresponded to residues near the relevant region in multiple proteins. 
        - Pocket selection for docking has a few modes, the default one is `close_all`. It works as follows:
            - First, selected indices from the global MSA are matched to residues in the given protein
            - Then, a centroid is computed from the selected residues, which is hopefully close to the external part of the membrane
            - For each pocket, it computes the distance between the closest atom in the pocket and the centroid
            - If the pocket is close enough (<= 20 A), it is considered for docking
        - use `python run_docking.py -h` to get information about the other modes
    - Pocket size is computed automatically according to the p2rank pocket prediction
    - Each protonated receptor-pocket pair is then passed to `mk_prepare_receptor.py`
    - Prepared receptors and vina configs are stored in `data/docking_files` by default
3. Dock ligands into predicted pockets with `source/run_docking.py`. 


    - After the receptors are prepared, it runs the AutoDock Vina 4 using the Vina forcefield on each combination of receptor and ligand.
    - Outputs are deposited in the `output/` folder
        - each pocket has its own subfolder (e.g. `output/output/A0A067XG43_p1` for pocket with rank 1 from the P2rank prediction) containing the `.pdbqt` files of docked ligands
