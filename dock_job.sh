#!/bin/bash

#PBS -N docking
#PBS -l select=1:ncpus=2:mem=8gb:scratch_local=20gb:cluster=alfrid
#PBS -l walltime=48:00:00
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

# Install RDKit via conda-forge (required dependency for molscrub)
echo "=== INSTALLING RDKIT ==="
if python -c "import rdkit" 2>/dev/null; then
    echo "RDKit already installed"
    python -c "import rdkit; print('RDKit version:', rdkit.__version__)"
else
    echo "Installing RDKit from conda-forge..."
    # Use mamba if available (faster), otherwise conda
    if command -v mamba >/dev/null 2>&1; then
        mamba install -y rdkit -c conda-forge || {
            echo "ERROR: Failed to install RDKit via mamba" >&2
            exit 1
        }
    else
        conda install -y rdkit -c conda-forge || {
            echo "ERROR: Failed to install RDKit via conda" >&2
            exit 1
        }
    fi
    echo "RDKit installed successfully"
fi

# Install molscrub if not already installed
echo "=== INSTALLING MOLSCRUB ==="
# Always add molscrub to PYTHONPATH to ensure subprocess calls can find it
export PYTHONPATH="$HOMEDIR/molscrub:$PYTHONPATH"

# Check if molscrub is properly installed (can import Scrub class)
if python -c "from molscrub import Scrub" 2>/dev/null; then
    echo "molscrub already installed and working"
    python -c "import molscrub; print('molscrub module:', molscrub.__file__ if hasattr(molscrub, '__file__') else 'namespace package')"
else
    # Check if molscrub directory exists
    if [ ! -d "$HOMEDIR/molscrub" ]; then
        echo "ERROR: molscrub directory not found at $HOMEDIR/molscrub" >&2
        exit 1
    fi
    
    echo "Installing molscrub from $HOMEDIR/molscrub"
    # Install in editable mode (rdkit already installed via conda, pip will install other deps like rich, numpy)
    pip install -e "$HOMEDIR/molscrub" || {
        echo "ERROR: Failed to install molscrub" >&2
        exit 1
    }
    # Verify installation worked
    echo "Verifying molscrub installation..."
    
    if python -c "from molscrub import Scrub" 2>/dev/null; then
        echo "molscrub installed successfully"
        python -c "import molscrub; print('molscrub module:', molscrub.__file__ if hasattr(molscrub, '__file__') else 'namespace package')"
    else
        echo "WARNING: Direct import failed, but PYTHONPATH set as fallback" >&2
        # Try importing with explicit path
        python -c "import sys; sys.path.insert(0, '$HOMEDIR/molscrub'); from molscrub import Scrub; print('molscrub works with explicit path')" 2>&1 || {
            echo "ERROR: molscrub still not importable" >&2
            python -c "from molscrub import Scrub" 2>&1 | head -10 || true
            echo "WARNING: Continuing - subprocess calls may fail" >&2
        }
    fi
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

python prepare_receptors.py

# Use splits directory if PBS_ARRAY_INDEX is set, otherwise use default docking_files
# Speed options: -e for exhaustiveness (lower=faster, default=32), -n for num_modes (lower=faster, default=9)
# Example for faster docking: python run_docking.py -e 16 -n 5
if [ -n "$PBS_ARRAY_INDEX" ] && [ -d "../data/splits_$PBS_ARRAY_INDEX" ]; then
    python run_docking.py -d ../data/splits_$PBS_ARRAY_INDEX -e 16 -n 3
else
    echo "PBS_ARRAY_INDEX not set or splits directory doesn't exist, using default docking_files directory"
    python run_docking.py -d ../data/docking_files
fi
