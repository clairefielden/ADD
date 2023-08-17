#!/bin/sh
#SBATCH --account compsci
#SBATCH --partition=ada
#SBATCH --time=10:00:00
#SBATCH --nodes=1 --ntasks=4
#SBATCH --job-name="Explore_High"
#SBATCH --mail-user=fldcla001@myuct.ac.za
#SBATCH --mail-type=ALL

module load python/anaconda-python-3.7
cd /scratch/fldcla001/ADD/Reinvent
conda env create -f reinvent.yml
cd ..
cd DockStream
conda env create -f environment.yml
cd ..
mkdir results
python ./Reinvent/input.py ./src/exploreHigh/exploreH_a.py
python ./Reinvent/input.py ./src/exploreHigh/exploreH_b.py
python ./Reinvent/input.py ./src/exploreHigh/exploreH_c.py
python ./Reinvent/input.py ./src/exploreHigh/exploreH_d.py
python ./Reinvent/input.py ./src/exploreHigh/exploreH_e.py
python ./Reinvent/input.py ./src/exploreHigh/exploreH_f.py