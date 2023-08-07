# load dependencies
import os
import re
import json
import tempfile

# --------- change these path variables as required
reinvent_dir = os.path.expanduser("./Reinvent")
reinvent_env = os.path.expanduser("/usr/local/envs/reinvent.v3.2")
output_dir = os.path.expanduser("./AutoDock_Vina_demo")

try:
    os.mkdir(output_dir)
except FileExistsError:
    pass

# DockStream variables
dockstream_dir = os.path.expanduser("./DockStream")
dockstream_env = os.path.expanduser("/usr/local/envs/DockStream")
# generate the path to the DockStream entry points
docker_path = os.path.join(dockstream_dir, "docker.py")
vina_binary_location = "./autodock_vina_1_1_2_linux_x86"
target_preparator = dockstream_dir + "/target_preparator.py"

# generate the paths to the files shipped with this implementation
apo_1UYD_path = "./DockStreamCommunity/data/1UYD/1UYD_apo.pdb"
reference_ligand_path = "./DockStreamCommunity/data/1UYD/PU8.pdb"
smiles_path = "./DockStreamCommunity/data/1UYD/ligands_smiles.txt"

# generate output paths for the configuration file, the "fixed" PDB file and the "Gold" receptor
target_prep_path = output_dir + "/ADV_target_prep.json"
fixed_pdb_path = output_dir + "/ADV_fixed_target.pdb"
adv_receptor_path = output_dir + "/ADV_receptor.pdbqt"
log_file_target_prep = output_dir + "/ADV_target_prep.log"
log_file_docking = output_dir + "/ADV_docking.log"

# generate output paths for the configuration file, embedded ligands, the docked ligands and the scores
docking_path = output_dir + "/ADV_docking.json"
ligands_conformers_path = output_dir + "/ADV_embedded_ligands.sdf"
ligands_docked_path = output_dir + "/ADV_ligands_docked.sdf"
ligands_scores_path = output_dir + "/ADV_scores.csv"

# execute this in a command-line environment after replacing the parameters
#!{dockstream_env}/bin/python {target_preparator} -conf {target_prep_path}
#python3 ./DockStream/target_preparator.py -conf ./AutoDock_Vina_demo/ADV_target_prep.json
#!head -n 25 {adv_receptor_path}
#head -n 25 ./AutoDock_Vina_demo/ADV_receptor.pdbqt


# execute this in a command-line environment after replacing the parameters
#!{dockstream_env}/bin/python {target_preparator} -conf {target_prep_path}
#python3 ./DockStream/target_preparator.py -conf ./AutoDock_Vina_demo/ADV_target_prep.json
#!head -n 25 {adv_receptor_path}
#head -n 25 ./AutoDock_Vina_demo/ADV_receptor.pdbqt

# load the smiles (just for illustrative purposes)
# here, 15 moleucles will be used
with open(smiles_path, 'r') as f:
    smiles = [smile.strip() for smile in f.readlines()]
print(smiles)

#!cat {log_file_target_prep}

# specify the embedding and docking JSON file as a dictionary and write it out
ed_dict = {
  "docking": {
    "header": {                                         # general settings
      "logging": {                                      # logging settings (e.g. which file to write to)
        "logfile": log_file_docking
      }
    },
    "ligand_preparation": {                             # the ligand preparation part, defines how to build the pool
      "embedding_pools": [
        {
          "pool_id": "Corina_pool",                     # here, we only have one pool
          "type": "Corina",
          "parameters": {
            "prefix_execution": "module load corina"    # only required, if a module needs to be loaded to execute "Corina"
          },
          "input": {
            "standardize_smiles": False,
            "type": "smi",
            "input_path": smiles_path
          },
          "output": {                                   # the conformers can be written to a file, but "output" is
                                                        # not required as the ligands are forwarded internally
            "conformer_path": ligands_conformers_path, 
            "format": "sdf"
          }
        }
      ]
    },
    "docking_runs": [
    {
      "backend": "AutoDockVina",
      "run_id": "AutoDockVina",
      "input_pools": ["Corina_pool"],
      "parameters": {
        "binary_location": vina_binary_location,        # absolute path to the folder, where the "vina" binary
                                                        # can be found
        "parallelization": {
          "number_cores": 4
        },
        "seed": 42,                                     # use this "seed" to generate reproducible results; if
                                                        # varied, slightly different results will be produced
        "receptor_pdbqt_path": [adv_receptor_path],     # paths to the receptor files
        "number_poses": 2,                              # number of poses to be generated
        "search_space": {                               # search space (cavity definition); see text
          "--center_x": 3.3,
          "--center_y": 11.5,
          "--center_z": 24.8,
          "--size_x": 15,
          "--size_y": 10,
          "--size_z": 10
        }
      },
      "output": {
        "poses": { "poses_path": ligands_docked_path },
        "scores": { "scores_path": ligands_scores_path }
}}]}}

with open(docking_path, 'w') as f:
    json.dump(ed_dict, f, indent=2)

# print out path to generated JSON
print(docking_path)

# execute this in a command-line environment after replacing the parameters
#!{dockstream_env}/bin/python {docker} -conf {docking_path} -print_scores


