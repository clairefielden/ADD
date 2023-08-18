#!/bin/sh
#SBATCH --account compsci
#SBATCH --partition=a100
#SBATCH --gres=gpu:a100-1g-5gb:1
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
python ./Reinvent/input.py ./src/exploitHigh/exploitH_a.py
python ./Reinvent/input.py ./src/exploitHigh/exploitH_b.py
python ./Reinvent/input.py ./src/exploitHigh/exploitH_c.py
python ./Reinvent/input.py ./src/exploitHigh/exploitH_d.py
python ./Reinvent/input.py ./src/exploitHigh/exploitH_e.py
python ./Reinvent/input.py ./src/exploitHigh/exploitH_f.py