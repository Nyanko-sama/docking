#!/bin/bash

#PBS -N docking
# Alfrid lebo ztadial sa mi dal spustit Mamba(conda) environment, z nejakoho dovodu to nefungovalo na elmo clusteri
#PBS -l select=1:ncpus=2:mem=4gb:scratch_local=10gb:cluster=alfrid
#PBS -l walltime=6:00:00
#PBS -m ae

# Path to home directory
HOMEDIR=/storage/praha1/home/nexuso1/

# Path to singularity image we want to use
#IMAGE=/cvmfs/singularity.metacentrum.cz/NGC/PyTorch\:24.10-py3.SIF

# Probably not necessary
export TMPDIR=$SCRATCHDIR
export CACHEDIR=$SCRATCHDIR
export LOCALCACHEDIR=$SCRATCHDIR

cd $HOMEDIR
cd prankdock
module add mambaforge
mamba env create -f environment.yaml
mamba activate dock
cd source
echo "Running full pipeline with array index $PBS_ARRAY_INDEX"
python run_p2rank.py
python prepare_ligands.py --timeout 240
python prepare_receptors.py
python run_docking.py -d ../data/splits_$PBS_ARRAY_INDEX