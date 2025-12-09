#!/bin/bash

#PBS -N docking
#PBS -l select=1:ncpus=2:mem=8gb:scratch_local=20gb:cluster=alfrid
#PBS -l walltime=6:00:00
#PBS -m ae

# Path to docking project root - adjust if needed
# On server, the folder is called "docking" (not "prankdock")
# Try multiple possible locations
if [ -d "/storage/plzen1/home/nelia_k/docking" ]; then
    HOMEDIR=/storage/plzen1/home/nelia_k/docking
elif [ -d "/storage/praha1/home/nelia_k/docking" ]; then
    HOMEDIR=/storage/praha1/home/nelia_k/docking
elif [ -d "$HOME/docking" ]; then
    HOMEDIR="$HOME/docking"
else
    echo "ERROR: Cannot find docking directory." >&2
    echo "Checked:" >&2
    echo "  /storage/plzen1/home/nelia_k/docking" >&2
    echo "  /storage/praha1/home/nelia_k/docking" >&2
    echo "  $HOME/docking" >&2
    echo "Current HOME: $HOME" >&2
    echo "Please set HOMEDIR manually in the script." >&2
    exit 1
fi

echo "Using HOMEDIR: $HOMEDIR"

export TMPDIR=$SCRATCHDIR
export CACHEDIR=$SCRATCHDIR
export LOCALCACHEDIR=$SCRATCHDIR

cd "$HOMEDIR" || {
    echo "ERROR: Failed to cd to $HOMEDIR" >&2
    exit 1
}

# Load mambaforge and enable "mamba activate" in this non-interactive shell
module add mambaforge

# Initialize mamba/conda for bash
if command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
    conda activate dock
elif command -v mamba >/dev/null 2>&1; then
    eval "$(mamba shell hook -s bash)"
    mamba activate dock
else
    echo "ERROR: Neither conda nor mamba found after loading module" >&2
    exit 1
fi

# Verify activation worked
if ! command -v python >/dev/null 2>&1; then
    echo "ERROR: Python not found after environment activation" >&2
    exit 1
fi

# Debug info to confirm env is correct
echo "=== PYTHON INFO ==="
which python
python -c "import sys; print('python exe:', sys.executable)"

# Install molscrub if not already installed
echo "=== INSTALLING MOLSCRUB ==="
if python -c "import molscrub" 2>/dev/null; then
    echo "molscrub already installed"
    python -c "import molscrub; print('molscrub module:', molscrub.__file__)"
else
    echo "Installing molscrub from $HOMEDIR/molscrub"
    pip install -e "$HOMEDIR/molscrub" || {
        echo "ERROR: Failed to install molscrub" >&2
        exit 1
    }
    python -c "import molscrub; print('molscrub module:', molscrub.__file__)"
fi

cd source || {
    echo "ERROR: Failed to cd to source directory" >&2
    exit 1
}

echo "Running full pipeline with array index $PBS_ARRAY_INDEX"
echo "Current directory: $(pwd)"

# Run pipeline steps with error checking
set -e  # Exit on any error

# if needed:
python run_p2rank.py --p2rank_path ../p2rank_2.5.1

python prepare_ligands.py --timeout 240

python prepare_receptors.py

python run_docking.py -d ../data/splits_$PBS_ARRAY_INDEX
