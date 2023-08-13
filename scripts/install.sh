#!/bin/bash

OUTPUT_DIR="CL_Results"

# Script for preparing Reinvent and Dockstream
red='\033[0;31m'
magenta='\033[0;35m'
# Clear the color after that
clear='\033[0m'

echo -e "${red}Note:${clear} Ensure that conda has been installed correctly."
# Prepare Reinvent environment
conda shell init
cd Reinvent && conda env create -f reinvent.yml
cd ..
# Prepare Dockstream Environments
cd DockStream && conda env create -f environment.yml
cd DockStream && conda env create -f environment_full.yml
cd ..
#Prepare RLMOL Environment
cd rlmol && conda env create -f env.yml
cd ..
# Init bash
conda init bash
#install ccdc conda channel
#conda install --channel=https://conda.ccdc.cam.ac.uk csd-python-api


# Create output directory
if ! test -d "$OUTPUT_DIR"; then
  mkdir $OUTPUT_DIR 
fi

echo -e "${red}Note:${clear} Your shell needs to be ${magenta}REBOOTED${clear}."