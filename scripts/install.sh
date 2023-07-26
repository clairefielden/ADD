#!/bin/bash

OUTPUT_DIR="Dockstream_CL"

# Script for preparing Reinvent and Dockstream
red='\033[0;31m'
magenta='\033[0;35m'
# Clear the color after that
clear='\033[0m'

echo -e "${red}Note:${clear} Ensure that conda has been installed correctly."

# Prepare Reinvent environment
cd Reinvent && conda env create -f reinvent.yml
cd ..
# Prepare Dockstream Environment
cd DockStream && conda env create -f environment.yml
cd ..
# Init bash
conda init bash

# Create output directory
if ! test -d "$OUTPUT_DIR"; then
  mkdir $OUTPUT_DIR 
fi

echo -e "${red}Note:${clear} Your shell needs to be ${magenta}REBOOTED${clear}."