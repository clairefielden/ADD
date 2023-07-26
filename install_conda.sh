#!/bin/bash
# Go to https://repo.anaconda.com/archive/ to choose the appropriate conda installation
CONDA_VERSION="Anaconda3-2023.07-1-Linux-x86_64.sh"

# MacOSX version: Anaconda3-2023.07-1-MacOSX-x86_64.sh

# Color variables
red='\033[0;31m'
magenta='\033[0;35m'
# Clear the color after that
clear='\033[0m'

# Download file if not downloaded already
if ! test -f "$CONDA_VERSION"; then
  wget https://repo.continuum.io/archive/$CONDA_VERSION  
fi
# Give execute permission
chmod +x $CONDA_VERSION
# Install
bash ./$CONDA_VERSION -b -f -p /usr/local/

# Verify Installation
which conda
conda --version

echo -e "${red}Note:!${clear}"
echo -e "You may want to set the defaults of your conda environment. You can do that using:"
echo -e "${magenta}conda install --channel defaults conda python=3.10 --yes${clear}"
echo -e "${magenta}conda update --channel defaults --all --yes${clear}"