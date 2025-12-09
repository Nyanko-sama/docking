#!/bin/bash

#PBS -N docking
#PBS -l select=1:ncpus=2:mem=8gb:scratch_local=20gb:cluster=alfrid
#PBS -l walltime=6:00:00
#PBS -m ae

# Path to docking project root
HOMEDIR=/storage/praha1/home/nelia_k/docking

export TMPDIR=$SCRATCHDIR
export CACHEDIR=$SCRATCHDIR
export LOCALCACHEDIR=$SCRATCHDIR

cd "$HOMEDIR/prankdock"

# Load mambaforge and enable "mamba activate" in this non-interactive shell
module add mambaforge
eval "$(mamba shell hook -s bash)"
mamba activate dock

# Debug info to confirm env is correct (keep for now)
echo "=== PYTHON INFO ==="
which python
python -c "import sys; print('python exe:', sys.executable)"
python -c "import molscrub; print('molscrub module:', molscrub.__file__)"

cd source
echo "Running full pipeline with array index $PBS_ARRAY_INDEX"

# if needed:
python run_p2rank.py --p2rank_path ../p2rank_2.5.1

python prepare_ligands.py --timeout 240
python prepare_receptors.py
python run_docking.py -d ../data/splits_$PBS_ARRAY_INDEX
