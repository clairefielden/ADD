# load dependencies
import os
import re
import json
import tempfile

# --------- change these path variables as required
abs_path_to_file = "/home/cfielden/lustre/fldcla001/ADD"

reinvent_dir = os.path.join(abs_path_to_file,"Reinvent")
reinvent_env = os.path.expanduser("/home/cfielden/.conda/envs/reinvent.v3.2")
output_dir = os.path.join(abs_path_to_file,"ADV_Results")

# if required, generate a folder to store the results
try:
    os.mkdir(output_dir)
except FileExistsError:
    pass

# DockStream variables
dockstream_dir = os.path.join(abs_path_to_file,"DockStream")
dockstream_env = os.path.expanduser("/home/cfielden/.conda/envs/DockStream/bin/python")
# generate the path to the DockStream entry points
docker = os.path.join(dockstream_dir, "docker.py")
target_preparator = dockstream_dir + "/target_preparator.py"

# Vina docking variables
vina_binary_location = os.path.join(abs_path_to_file,"autodock_vina_1_1_2_linux_x86/bin")

docking_configuration_path = os.path.join(output_dir, "ADV_dock_2.json")

# generate the paths to the **PfPI4K** receptor
adv_receptor_path = os.path.join(abs_path_to_file,"models/pfpi4k-no-ligand.pdbqt")
log_file_target_prep = output_dir + "/ADV_target_prep.log"
log_file_docking = output_dir + "/ADV_docking.log"

# generate output paths for the configuration file, embedded ligands, the docked ligands and the scores
docking_path = output_dir + "/ADV_dock_2.json"
ligands_conformers_path = output_dir + "/ADV_embedded_ligands_2.sdf"
ligands_docked_path = output_dir + "/ADV_ligands_docked_2.sdf"
ligands_scores_path = output_dir + "/ADV_scores_2.csv"

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
          "pool_id": "RDkit_pool",                     # here, we only have one pool
          "type": "RDkit",
          "parameters": {
            "removeHs": False,
            "coordinate_generation": {
                  "method": "UFF",
                  "maximum_iterations": 450
              }
          },
          "input": {
            "standardize_smiles": False,
          },
          "output": {                                   # the conformers can be written to a file, but "output" is not required as the ligands are forwarded internally
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
      "input_pools": ["RDkit_pool"],
      "parameters": {
        "binary_location": vina_binary_location,        # absolute path to the folder, where the "vina" binary
                                                        # can be found
        "parallelization": {
          "number_cores": 4
        },
        "seed": 42,                                     # use this "seed" to generate reproducible results; if
                                                        # varied, slightly different results will be produced
        "receptor_pdbqt_path": [adv_receptor_path],     # paths to the receptor files
        "number_poses": 20,                              # number of poses to be generated
        "search_space": {                               # search space (cavity definition); see text
          "--center_x": 6.980,
          "--center_y": 335.600,
          "--center_z": 8.060,
          "--size_x": 40,
          "--size_y": 40,
          "--size_z": 40
        }
      },
      "output": {
        "poses": { "poses_path": ligands_docked_path, "overwrite": False },
        "scores": { "scores_path": ligands_scores_path, "overwrite": False }
}}]}}

with open(docking_path, 'w') as f:
    json.dump(ed_dict, f, indent=2)

# print out path to generated JSON
print(docking_path)

# initialize the dictionary
configuration = {
    "version": 3,                                 # we are going to use REINVENT's newest release
    "run_type": "curriculum_learning",            # other run types: "sampling", "reinforcement_learning",
    #                  "transfer_learning",
    #                  "scoring" and "create_model"
    "model_type": "default"
}

# add block to specify whether to run locally or not and
# where to store the results and logging
configuration["logging"] = {
    "sender": "http://0.0.0.1",            # only relevant if "recipient" is set to "remote"
    "recipient": "local",                  # either to local logging or use a remote REST-interface
    "logging_frequency": 1000,              # log every x-th steps
    "logging_path": os.path.join(output_dir, "progress.log"), # load this folder in tensorboard
    "result_folder": os.path.join(output_dir, "results_2"),     # will hold the compounds (SMILES) and summaries
    "job_name": "Explorative Behaviour with High Threshold",         # set an arbitrary job name for identification
    "job_id": "ErHT"                       # only relevant if "recipient" is set to a specific REST endpoint
}

# add the "parameters" block
configuration["parameters"] = {}

# First add the paths to the Prior, Agent, and set the curriculum type to automated
configuration["parameters"]["prior"] = "./models/random.prior.new"
configuration["parameters"]["agent"] = "./models/random.prior.new"
configuration["parameters"]["curriculum_type"] = "automated"

# set up the Curriculum Strategy
configuration["parameters"]["curriculum_strategy"] = {
    "name": "user_defined",         # denotes that the order of Curriculum Objectives is defined by the user
    "max_num_iterations": 20000,     # denotes the total number of epochs to spend in the Curriculum Phase
    # if by the end of the total epochs the last Curriculum Objective is not
    # satisfied (based on the agent achieving a score >= threshold), the run stops
    "batch_size": 128,              # specifies how many molecules are generated per epoch
    "learning_rate": 0.0001,        # sets how strongly the agent is influenced by each epoch
    "sigma": 128,                   # used to calculate the "augmented likelihood", see publication
    "learning_strategy": {
        "name": "dap_single_query",
        "parameters": {
            "sigma": 120
        }
    },
    "diversity_filter": {
        "name": "IdenticalMurckoScaffold",         # other options are: "IdenticalTopologicalScaffold",
        #                    "IdenticalMurckoScaffold", and "ScaffoldSimilarity"

        "bucket_size": 25,          # the bin size; penalization will start once this is exceeded
        "minscore": 0.4,            # the minimum total score to be considered for binning
        "minsimilarity": 0.4        # the minimum similarity to be placed into the same bin
    },
    "inception": {
        "smiles": [],               # fill in a list of SMILES here that can be used (or leave empty)
        "memory_size": 100,         # sets how many molecules are to be remembered
        "sample_size": 10           # how many are to be sampled each epoch from the memory for experience replay
    },
    # Curriculum Objectives are all the scoring functions that are to be sequentially activated
    "curriculum_objectives": [{
        # 1st scoring function/curriculum objective is to obtain a high tanimoto similarity to the existing Malaria drugs
        #napthyridine = 'C1=CC2=C(C=CC=N2)N=C1'
        #aminopyridine = 'C1=CC=NC(=C1)N'
        #imidazopyridazine = 'C1=CC2=NC=CN2N=C1'
        "scoring_function": {
            "name": "custom_product",     # this is our default one (alternative: "custom_sum")
            "parallel": False,
            "parameters": [{
                "component_type": "tanimoto_similarity",
                "name": "Tanimoto similarity",         # arbitrary name for the component
                "weight": 1,                           # the weight of the component (default: 1)
                "specific_parameters": {
                    "smiles": ["C1=CC2=C(C=CC=N2)N=C1", "C1=CC=NC(=C1)N","C1=CC2=NC=CN2N=C1"], # a list of SMILES can be provided
                }
            }]             # the weight of the component (default: 1)
        },
        "score_threshold": 0.75            # agent must achieve an average score of this before
        # progressing to the next Curriculum Objective
    },
        # 2nd scoring function/curriculum objective is to obtain high QED
        {
            "scoring_function": {
                # add component: calculate the QED drug-likeness score (using RDkit)
                "name":"custom_product",
                "parallel":False,
                "parameters":[{
                    "component_type": "qed_score",
                    "name": "QED Score",# arbitrary name for the component
                    "weight": 1 # the weight of the component (default: 1)
                }]
            },
            "score_threshold": 0.75
        },
        # 3rd scoring function: Matching Substructure
        {
            "scoring_function": {
            "name": "custom_product",
            "parallel": False,
            "parameters": [{
                "component_type": "matching_substructure",
                "name": "Full Substructure",
                "specific_parameters": {
                    "smiles": ["C1=CC2=C(C=CC=N2)N=C1", "C1=CC=NC(=C1)N","C1=CC2=NC=CN2N=C1"]
                },
                "weight": 1}]
            },
        "score_threshold": 0.75
        }
    ]
}

# set up the Production Strategy
configuration["parameters"]["production_strategy"] = {
    "name": "standard",
    "retain_inception": True,       # option to retain the inception from the Curriculum Phase
    # retain it here since the last Curriculum Objective is the same as
    # Production Objective. Previous top compounds will be relevant

    "number_of_steps": 5000,         # number of epochs to run the Production Phase
    "batch_size": 128,              # specifies how many molecules are generated per epoch
    "learning_rate": 0.0001,        # sets how strongly the agent is influenced by each epoch
    "sigma": 128,                   # used to calculate the "augmented likelihood", see publication
    "learning_strategy": {
        "name": "dap_single_query",
        "parameters": {
            "sigma": 120
        }
    },
    "diversity_filter": {
        "name": "NoFilter",         # other options are: "IdenticalTopologicalScaffold",
        #                    "IdenticalMurckoScaffold"", and "ScaffoldSimilarity"

        "bucket_size": 25,          # the bin size; penalization will start once this is exceeded
        "minscore": 0.4,            # the minimum total score to be considered for binning
        "minsimilarity": 0.4        # the minimum similarity to be placed into the same bin
    },
    "inception": {
        "smiles": [],               # fill in a list of SMILES here that can be used (or leave empty)
        "memory_size": 100,         # sets how many molecules are to be remembered
        "sample_size": 10           # how many are to be sampled each epoch from the memory for experience replay
    },
    # the Production Objective contains the final scoring function to be activated
    # here, it is the same scoring function as the last Curriculum Objective
    # as we want to continue sampling the target substructure
    "scoring_function": {
        "name": "custom_product",
        "parallel": False,
        "parameters": [{
            "component_type": "dockstream",
            "name": "Autodock Vina Docking",
            "weight": 1,
            "specific_parameters": {
                "transformation": {
                    "transformation_type": "reverse_sigmoid",         # lower/more negative Autodock scores are better - use reverse
                    # sigmoid transformation
                    "low": -12,
                    "high": -8,
                    "k": 0.25
                },
                "configuration_path": docking_configuration_path,
                "docker_script_path": docker,
                "environment_path": dockstream_env
            },
            "weight": 1}]
    }
}

configuration_JSON_path = os.path.join(output_dir, "ADV_2.json")
with open(configuration_JSON_path, 'w') as f:
    json.dump(configuration, f, indent=4)
