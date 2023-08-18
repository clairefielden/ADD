#!/bin/sh
#PBS -n explore_high
#PBS -l select=1
#PBS -P CSCI1142
#PBS -q serial
#PBS -l walltime=40:00:00
#PBS -o /mnt/lustre/cfielden/fldcla001/ADD/results/explore_high.out
#PBS -m abe
#PBS -M fldcla001@myuct.ac.za

module load chpc/chem/anaconda3-2021.11
cd /scratch/fldcla001/ADD/Reinvent
mkdir results
conda activate reinvent.v3.2
python ./Reinvent/input.py ./src/exploreHigh/exploreH_a.py
python ./Reinvent/input.py ./src/exploreHigh/exploreH_b.py
python ./Reinvent/input.py ./src/exploreHigh/exploreH_c.py
python ./Reinvent/input.py ./src/exploreHigh/exploreH_d.py
python ./Reinvent/input.py ./src/exploreHigh/exploreH_e.py
python ./Reinvent/input.py ./src/exploreHigh/exploreH_f.py