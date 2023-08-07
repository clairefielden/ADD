#!/bin/bash
# Go to https://repo.anaconda.com/archive/ to choose the appropriate conda installation
CONDA_VERSION="Miniconda3-py37_4.10.3-Linux-x86_64.sh"

# Color variables
red='\033[0;31m'
magenta='\033[0;35m'
# Clear the color after that
clear='\033[0m'

# Download file if not downloaded already
if ! test -f "$CONDA_VERSION"; then
  wget -O miniconda.sh https://repo.anaconda.com/miniconda/$CONDA_VERSION  
fi
# Give execute permission
chmod +x miniconda.sh
#install
bash ./miniconda.sh -b -f -p /usr/local/
#remove
rm miniconda.sh

#configure micromamba
conda config --add channels conda-forge
conda install -y mamba
mamba update -qy --all
mamba clean -qafy

# Verify Installation
which conda
conda --version