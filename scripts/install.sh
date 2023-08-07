#!/bin/bash

OUTPUT_DIR="CL_Results"

# Script for preparing Reinvent and Dockstream
red='\033[0;31m'
magenta='\033[0;35m'
# Clear the color after that
clear='\033[0m'

echo -e "${red}Note:${clear} Ensure that conda has been installed correctly."
# Prepare Reinvent environment
mamba shell init
cd Reinvent && mamba env create -f reinvent.yml
cd ..
# Prepare Dockstream Environment
cd DockStream && mamba env create -f environment.yml
cd ..
#Prepare RLMOL Environment
cd rlmol && mamba env create -f env.yml
# Init bash
mamba init bash

# Create output directory
if ! test -d "$OUTPUT_DIR"; then
  mkdir $OUTPUT_DIR 
fi

echo -e "${red}Note:${clear} Your shell needs to be ${magenta}REBOOTED${clear}."