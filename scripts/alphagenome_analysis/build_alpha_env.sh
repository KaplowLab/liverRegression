#!/bin/bash
#SBATCH --job-name=BuildAlphaEnv
#SBATCH -p
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00         
#SBATCH --output

# Initialize and create env
source ~/miniconda3/bin/activate
conda create -n alphagenome_env python=3.11 -y
conda activate alphagenome_env

# Configure channels
conda config --env --add channels conda-forge
conda config --env --set channel_priority strict

# Explicitly install before alphagenome install
conda install -y \
    gcc_linux-64=12 \
    gxx_linux-64=12 \
    binutils \
    scipy \
    openblas \
    h5py \
    pyarrow \
    biopython \
    pandas

pip install alphagenome

# Check proper installation
echo "VERIFICATION CHECK:"
python -c "import alphagenome; print('alphagenome version', alphagenome.__version__)"
