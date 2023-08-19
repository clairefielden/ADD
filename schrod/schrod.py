# load dependencies

import os
import re
import json
import tempfile

# --------- change these path variables as required
abs_path_to_file = "/home/cfielden/lustre/fldcla001/ADD"

reinvent_dir = os.path.join(abs_path_to_file,"Reinvent")
reinvent_env = os.path.expanduser("/home/cfielden/.conda/envs/reinvent.v3.2")
output_dir = os.path.join(abs_path_to_file,"schrod/results")
log_file_docking = output_dir+"/schrod_docking.log"

try:
    os.makedirs(output_dir)
except FileExistsError:
    pass

# Glide docking variables
# DockStream variables
dockstream_dir = os.path.join(abs_path_to_file,"DockStream")
dockstream_env = os.path.join("/home/cfielden/.conda/envs/DockStream/bin/python")
# generate the path to the DockStream entry points
docker_path = os.path.join(dockstream_dir, "docker.py")
grid_file_path = os.path.join(abs_path_to_file,"models/pfpi4k_grid.zip")
output_ligands_docked_poses_path = os.path.join(abs_path_to_file,"schrod/results/docked_poses")
output_ligands_docking_scores_path = os.path.join(abs_path_to_file,"schrod/results/docking_scores")

try:
    os.mkdir(output_ligands_docked_poses_path)
except FileExistsError:
    pass

try:
    os.mkdir(output_ligands_docking_scores_path)
except FileExistsError:
    pass

docking_configuration_path = os.path.join(abs_path_to_file, "schrod/glide_conf.json")

# specify the embedding and docking JSON file as a dictionary and write it out
ed_dict = {
  "docking": {
    "header": {                                   # general settings
      "logging": {
          "logfile": log_file_docking
      }
    },
    "ligand_preparation": {                       # the ligand preparation part, defines how to build the pool
      "embedding_pools": [
        {
          "pool_id": "Ligprep_pool",
          "type": "Ligprep",
          "parameters": {
            "prefix_execution": "module load schrodinger/2019-4",
            "parallelization": {
                "number_cores": 4
            },
            "use_epik": {
              "target_pH": 7.0,
              "pH_tolerance": 2.0
            },
            "force_field": "OPLS3e"
          },
          "input": {
            "standardize_smiles": False,
            "input_path": "console"                              # expected input is a text file with smiles
          }
        }
      ]
    },
    "docking_runs": [
        {
          "backend": "Glide",
          "run_id": "Glide_run",
        "input_pools": ["Ligprep_pool"],
        "parameters": {
          "prefix_execution": "module load schrodinger/2019-4", # will be executed before a program call
          "parallelization": {                              # if present, the number of cores to be used
                                                            # can be specified
            "number_cores": 4
          },
          "glide_flags": {                                  # all all command-line flags for Glide here 
            "-HOST": "localhost"
          },
          "glide_keywords": {                               # add all keywords for the "input.in" file here
            "AMIDE_MODE": "trans",
            "EXPANDED_SAMPLING": "True",
            "GRIDFILE": grid_file_path,
            "NENHANCED_SAMPLING": "2",
            "POSE_OUTTYPE": "ligandlib_sd",
            "POSES_PER_LIG": "3",
            "POSTDOCK_NPOSE": "15",
            "POSTDOCKSTRAIN": "True",
            "PRECISION": "HTVS",
            "REWARD_INTRA_HBONDS": "True"
          }
        },
        "output": {
          "poses": { "poses_path": output_ligands_docked_poses_path},
          "scores": { "scores_path": output_ligands_docking_scores_path}
        }
      }]}}

with open(docking_configuration_path, 'w') as f:
    json.dump(ed_dict, f, indent=2)

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
    "logging_frequency": 100,              # log every x-th steps
    "logging_path": os.path.join(output_dir, "progress.log"), # load this folder in tensorboard
    "result_folder": os.path.join(output_dir, "results"),     # will hold the compounds (SMILES) and summaries
    "job_name": "Schrodinger Demo",         # set an arbitrary job name for identification
    "job_id": "Demo"                       # only relevant if "recipient" is set to a specific REST endpoint
}

# add the "parameters" block
configuration["parameters"] = {}

# First add the paths to the Prior, Agent, and set the curriculum type to automated
configuration["parameters"]["prior"] = os.path.join(abs_path_to_file,"models/random.prior.new")
configuration["parameters"]["agent"] = os.path.join(abs_path_to_file,"models/random.prior.new")
configuration["parameters"]["curriculum_type"] = "automated"

# set up the Curriculum Strategy
configuration["parameters"]["curriculum_strategy"] = {
    "name": "user_defined",         # denotes that the order of Curriculum Objectives is defined by the user
    "max_num_iterations": 1500,     # denotes the total number of epochs to spend in the Curriculum Phase
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
        "name": "NoFilter",         # other options are: "IdenticalTopologicalScaffold",
        #                    "IdenticalMurckoScaffold", "NoFilter" and "ScaffoldSimilarity"

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
        # 1st scoring function below
        "scoring_function": {
            "name": "custom_product",     # this is our default one (alternative: "custom_sum")
            "parallel": False,
            "parameters": [{
                "component_type": "qed_score",     # enforce the match to a given substructure
                "name": "QED Score",     # arbitrary name for the component
                "weight": 1}]             # the weight of the component (default: 1)
        },
        "score_threshold": 0.25            # agent must achieve an average score of this before
        # progressing to the next Curriculum Objective
    },
        # 2nd scoring function below
        {
            "scoring_function": {
                "name": "custom_product",
                "parallel": False,
                "parameters": [{
                    "component_type": "tanimoto_similarity",
                    "name": "Tanimoto Similarity",
                    "specific_parameters": {
                        "smiles": [
                            "C1=CC2=C(C=CC=N2)N=C1", "C1=CC=NC(=C1)N","C1=CC2=NC=CN2N=C1"
                        ]
                    },
                    "weight": 1}]
            },
            "score_threshold": 0.5
        },
        # 3rd scoring function: Matching Substructure
        {
            "scoring_function": {
                # add component: calculate the synthetic accessibility (using RDkit)
                "name":"custom_product",
                "parallel":True,
                "parameters":[{
                    "component_type": "jaccard_distance",
                    "name": "SA Score",# arbitrary name for the component
                    "weight": 1 # the weight of the component (default: 1)
                }]
            },
            "score_threshold": 5
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
           "name": "custom_product",                  # this is our default one (alternative: "custom_sum")
            "parallel": False,                         # sets whether components are to be executed
                                               # in parallel; note, that python uses "False" / "True"
                                               # but the JSON "false" / "true"

            # the "parameters" list holds the individual components
            "parameters": [

            # add component: use 
            {
              "component_type": "dockstream",                           # use DockStream as a Scoring Function component      
              "name": "Glide LigPrep Docking",                          # arbitrary name
              "weight": 1,
              "specific_parameters": {
                  "transformation": {
                      "transformation_type": "reverse_sigmoid",         # lower Glide scores are better - use reverse
                                                              # sigmoid transformation
                      "low": -11,
                      "high": -5,
                      "k": 0.25
                  },
        "configuration_path": docking_configuration_path,
        "docker_script_path": docker_path,
        "environment_path": dockstream_env
        }
    }]
    }
}

configuration_JSON_path = os.path.join(abs_path_to_file, "schrod/schrod_conf.json")
with open(configuration_JSON_path, 'w') as f:
    json.dump(configuration, f, indent=4)
