# Improving Simulated Molecular Docking using Curriculum Learning
> This work improves upon current Reinforcement Learning approaches 
towards the design of potential antimalarial drugs.
> A demonstration of how each aspect of this codebase was used for this research is provided [_here_](https://github.com/clairefielden/ADD/blob/41485ea9503113c026543628ce8b35f815ef0b5e/Curriculum_Learning_Dockstream_Demo.ipynb)

## Table of Contents
* [Contributions](#contributions)
* [Requirements](#platforms)
* [Experiments](#experiment)
* [Algorithms](#algorithm)
* [Chemistry](#chem)
* [Models](#model)
* [Results and Documentation](#results)
* [Future Work](#future)
* [Acknowledgements](#acknowledgements)
* [Contact](#contact)


## Contributions
In the context of RL-based hit-to-lead optimization, an RL agent cannot prioritize docking scores in isolation. Docking scores are returned from docking simulators, such as [Glide](https://www.schrodinger.com/products/glide), as a prediction of binding affinity. The goal is to maximize the binding affinity between the ligand generated by the RL agent and the PfPI4K receptor, in the case of designing antimalarials.

As docking simulations are computationally expensive in terms of time and resources, Curriculum Learning (CL) is investigated as a reward shaping mechanism. This work will:
1. Evaluate `Curriculum Learning` as a reward shaping mechanism to improve the productivity of an RL agent.
2. Provide an optimal curriculum for the generation of  an antagonist with high `binding affinity` and `synthesizability`.
3. Employ the [Schrödinger](https://www.schrodinger.com/) software to perform High-Throughput Virtual Screening for the PfPI4K target receptor.
4. Add an additional scoring function to [REINVENT's](https://github.com/MolecularAI/ReinventCommunity) framework in the form of Ertl and Schuffenhauer’s SA Score.

## Requirements
1. A license to the [Schrödinger Software Suite](https://www.schrodinger.com/)
2. Access to a GPU (for these experiments, the computational resources required were 3 GPUs and 30 CPU cores)
3. Generate the receptor grid using [Maestro](https://www.schrodinger.com/products/maestro)
   - Download a [PDB file](https://github.com/clairefielden/ADD/blob/d269dacc6855e29b3f2793c4c2fc2a7f2b1a090d/chem/4d0l.pdb) containing the interaction between the target receptor and a sample ligand
   - Use the Protein Preparation Wizard to create a grid around the target receptor's binding cavity
   - Download the receptor-grid file as a [zip file](https://github.com/clairefielden/ADD/blob/d269dacc6855e29b3f2793c4c2fc2a7f2b1a090d/models/pfpi4k_grid.zip)

4. Compile the configuration files
   - The configuration files for these experiments are provided in the `experiments` subdirectory.
   - For more information on the possibilities of docking backends and their `json` files: [DockStream](https://github.com/MolecularAI/DockStreamCommunity/blob/master/notebooks/demo_Glide.ipynb)
   - For more information on the possibilities for RL-based drug design and their `json` files: [REINVENT3.2](https://github.com/MolecularAI/ReinventCommunity/blob/master/notebooks/Reinforcement_Learning_Demo.ipynb)


## Experiments
### Installation
1. Ensure you are using a GPU in your Runtime Evironment
2. Clone this project's [GitHub Repository](https://github.com/clairefielden/ADD): <br>
```bash
$ git clone https://clairefielden/ADD.git
$ cd ADD
```
3. Install [conda](https://docs.anaconda.com/free/anaconda/install/index.html): <br>
```bash
$ module load chpc/python/anaconda/3-2021.11
```
4. Install [Schrödinger](https://www.schrodinger.com/): <br>
```bash
$ module load chpc/schrodinger/2019-4
```
5. Create the environments
```bash
$ conda env create -f ./ADD/DockStream/environment.yml
$ conda env create -f ./ADD/Reinvent/reinvent.yml
```

### Usage

#### All experiments were run on the [(CHPC)](https://www.chpc.ac.za/) with the scripts in the `experiments` directory:
1. `Agent1`: Exploitative and Naiive
2. `Agent2`: Exploitative and Sophisticated
3. `Agent3`: Explorative and Naiive
4. `Agent4`: Explorative and Sophisticated
5. `Benchmark`: The baseline RL experiment used as a point of comparison

#### Each experiment consists of the following files:
1. `agent*.job`: The job file used to submit the experiment to the cluster
2. `curriculum.json`: The JSON file used to instantiate the Curriculum, Production and Reinforcement Learning Strategy using [REINVENT](https://github.com/MolecularAI/ReinventCommunity)
3. `production.json`: The JSON file used to input SMILES generated by the agent into [DockStream](https://github.com/MolecularAI/DockStreamCommunity/blob/master/notebooks/demo_Glide.ipynb), which provides access to the [Glide](https://www.schrodinger.com/products/glide) and [LigPrep](https://www.schrodinger.com/products/ligprep) backends

#### Each job script will navigate into the appropriate directory and perform the following actions:

1. Activate the environments: <br>
```bash
$ source /apps/chpc/anaconda/3-2021.11/etc/profile.d/conda.sh
$ source activate /home/<name>/.conda/envs/reinvent.v3.2
```
2. Input the `curriculum.json` file into the [REINVENT](https://github.com/MolecularAI/ReinventCommunity) framework:
```
$ cd /mnt/lustre/users/<src_directory>
$ python ./ADD/Reinvent/input.py ./ADD/<experiments>/<experiment>/<CL_config_file>
```
3. During the *Production Phase*, the `production.json` file is input to [DockStream](https://github.com/MolecularAI/DockStreamCommunity/blob/master/notebooks/demo_Glide.ipynb) through the following commands:
```
$ conda activate DockStream
$ python ./ADD/DockStream/docker.py -conf production.json
```

## Algorithms
### jaccard_distance.py and fpscores.pkl.gz
In order to measure the synthesizability of the ligands generated by the agents in this research, the algorithm for SA Score was added as a [REINVENT](https://github.com/MolecularAI/ReinventCommunity) framework.  The equation is based on work done by [Peter Ertl and Ansgar Schuffenhauer](http://www.jcheminf.com/content/1/1/8) and is publicly available on RDKit's [GitHub Repository](https://github.com/rdkit/rdkit/blob/master/Contrib/SA_Score/sascorer.py).
To add SA Score to the *Curriculum Learning* framework as a scoring function:
  1. Copy the [jaccard_distance.py](https://github.com/clairefielden/ADD/blob/478e333503d9500fa8c1869d88dc01d9659998c2/algorithms/jaccard_distance.py) file in `algorithms`
  2. Navigate to `/home/<user_name>/.conda/envs/reinvent.v3.2/lib/python3.7/site-packages/reinvent_scoring/scoring/score_components/standard/jaccard_distance.py` and replace it with the copied file
  3. Ensure [fpscores.pkl.gz](https://github.com/clairefielden/ADD/blob/478e333503d9500fa8c1869d88dc01d9659998c2/algorithms/fpscores.pkl.gz) is in `/mnt/lustre/users/<src_directory>/ADD/algorithms` - this is the file used to calculate the presence of synthesizable fragments in the generated compound
### productivity.py
This is the algorithm used to measure the `productivity` of the agent based on the results from the experiments. To run this algorithm:
1. Download the `avg_nll_agent.csv` and `docking_score.csv` files generated by TensorBoard after the experiment has completed
2. Place them in the `results` file under the relevant experiment
3. Change the relevant file paths. The equations are based on the following variables:
    - `average convergence rate`: The rate at which the `nll` of the agent improves/degenerates
    - `average docking score`: The average score returned by [Glide](https://www.schrodinger.com/products/glide) throughout the *Production Phase*

### average_sa.py
This is the algorithm used to measure the average `Multi-Parameter Optimization` abilities of each agent. To run this algorithm:
1. After the experiment has completed, navigate to it `results` file under the relevant experiment
2. Change the file paths to point to the `ligand_scores.csv` file saved during the *Production Phase*
3. The equations are based on the following variables:
    - `GlideScore`: The GlideScore for each ligand
    - `smiles`: The SMILES String of the generated ligand
    - `SA Score`: The SA Score measured, as per [RDKit](https://github.com/rdkit/rdkit/blob/master/Contrib/SA_Score/sascorer.py)

## Chemistry

### 4d0l.pdb
The PDB file used to generate the receptor grid using [Maestro](https://www.schrodinger.com/products/maestro)

### best_ligand.tex
The `chemfig` of the best-scoring ligand in this work

### pre_existing_ligands
The SMILES Strings of the pre-existing antimalarials used for the first `Curriculum Objective`, which guided the agent towards favourable regions of chemical space using Tanimoto Similarity. This CSV file was input into [DockStream](https://github.com/MolecularAI/DockStreamCommunity/blob/master/notebooks/demo_Glide.ipynb) to obtain the `ligands_csv_input_embedded.sdf` file, as a result of performing molecular docking simulations with [Glide](https://www.schrodinger.com/products/glide)

## Models
The `models` directory contains the PfPI4K docking grid `zip` file and `.sh` file that were generated using [Maestro](https://www.schrodinger.com/products/maestro). They are used by the `production.json` files during molecular docking simulations.
 
### random.prior.new
A pretrained RNN making use of 3 stacked LSTM layers. It has been fully trained on the ChemBL dataset to possess a large vocabulary of SMILES Strings. This is the model used as the `prior`, which is pretrained in the case of *Curriculum Learning*. For the benchmarking algorithm, this is the agent used for a standard RL *Production Phase*.

## Results and Documentation
#### The experiments in the `experiments` subdirectory resulted in the following directory:
1. `Agent1`: Exploitative and Naiive
2. `Agent2`: Exploitative and Sophisticated
3. `Agent3`: Explorative and Naiive
4. `Agent4`: Explorative and Sophisticated
5. `Benchmark`: The baseline RL experiment used as a point of comparison

#### Each folder consists of the following:
1. `results`: 
   - `*.model`: The model of the agent when training was termined (at the end of its *Production Phase*
   - `ligand_scores.csv`: The scores for the ligands in its last *Production* episode
   - `ligands_docked.sdf`: The conformations of the ligands, including their coordinates relative to the target receptor and grid
2. `*_nll.csv`: The NLL loss of the agent throughout its *Production Phase*
3. `*_ds.csv`: The docking scores throughout the agent's *Production Phase*

#### You can access my [research paper](https://github.com/clairefielden/ADD/blob/4523faa1112a4d54c010ec1030715ef9dcf8cdfe/Documentation/FLDCLA001%20-%20ADD%20Final%20Paper.pdf) and [PDB validation](https://github.com/clairefielden/ADD/blob/4523faa1112a4d54c010ec1030715ef9dcf8cdfe/Documentation/4d0l_full_validation.pdf) of *Phosphatidylinositol 4-kinase III beta-PIK93 in a complex with Rab11a-GTP gammaS* in the `Documentation` folder

## Future Work
In this work, the *Production Phase* consisted of 100 episodes, and thus the true benefit of Diversity Oriented Synthesis was not exhibited. For exploration to improve standard RL hit-to-lead optimization, sophisticated
Curriculum Progression Criteria should be made use of to avoid mode collapse. Moreover, a longer Curriculum Phase is required to accommodate the combination of sophistication and exploration. This shows promise in proposing novel regions of chemical space. Such propositions would be worth the Curriculum Phase training time, which is computationally inexpensive, due to expected improvement in Production Phase convergence time and ligand quality.

## Acknowledgements
- This project was inspired by Rob Maccallum's work, RLMOL
- Thank you to my supervisor :)
- Many thanks to the CSIR’s Centre for High Performance Computing [(CHPC)](https://www.chpc.ac.za/), who provided computational resources and access to the Schrödinger Software Suite


## Contact
Created by [@clairefielden](FLDCLA001@myuct.ac.za) - feel free to contact me!


## License
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
