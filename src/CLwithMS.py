# load dependencies
import os
import re
import json
import tempfile

# --------- change these path variables as required
reinvent_dir = os.path.expanduser("./Reinvent")
reinvent_env = os.path.expanduser("/usr/local/envs/reinvent.v3.2")
output_dir = os.path.expanduser("./Dockstream_CL")

# DockStream variables
dockstream_dir = os.path.expanduser("./DockStream")
dockstream_env = os.path.expanduser("/usr/local/envs/DockStream")
# generate the path to the DockStream entry points
docker_path = os.path.join(dockstream_dir, "docker.py")

# Glide docking variables
grid_file_path = os.path.expanduser("./DockStream/1UYD_grid.zip")
output_ligands_docked_poses_path = os.path.expanduser("./Dockstream_CL/docked_poses")
output_ligands_docking_scores_path = os.path.expanduser("./Dockstream_CL/docking_scores")

try:
    os.mkdir(output_ligands_docked_poses_path)
except FileExistsError:
    pass

try:
    os.mkdir(output_ligands_docking_scores_path)
except FileExistsError:
    pass

docking_configuration_path = os.path.join(output_dir, "Glide_DockStream_Conf.json")

# specify the embedding and docking JSON file as a dictionary and write it out
ed_dict = {
    "docking": {
        "header": {                                   # general settings
            "environment": {
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
                            "number_cores": 2
                        },
                        "use_epik": {
                            "target_pH": 7.0,
                            "pH_tolerance": 2.0
                        },
                        "force_field": "OPLS3e"
                    },
                    "input": {
                        "standardize_smiles": False,
                        "type": "console"                     # input type "console" when using DockStream with REINVENT
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
                    "prefix_execution": "module load schrodinger/2020-4", # will be executed before a program call
                    "parallelization": {                                  # if present, the number of cores to be used
                        # can be specified
                        "number_cores": 2
                    },
                    "glide_flags": {                                  # all all command-line flags for Glide here
                        "-HOST": "localhost"
                    },
                    "glide_keywords": {                               # add all keywords for the "input.in" file here
                        # this is the minimum keywords that needs to be
                        # specified and represents a simple `Glide`
                        # docking configuration

                        "GRIDFILE": grid_file_path,
                        "POSE_OUTTYPE": "ligandlib_sd",
                        "PRECISION": "HTVS"
                    }
                },
                "output": {
                    "poses": { "poses_path": os.path.join(output_ligands_docked_poses_path, "docked_poses.sdf")},
                    "scores": { "scores_path": os.path.join(output_ligands_docking_scores_path, "docking_scores.csv")}
                }
            }]}}

with open(docking_configuration_path, 'w') as f:
    json.dump(ed_dict, f, indent=2)

# if required, generate a folder to store the results
try:
    os.mkdir(output_dir)
except FileExistsError:
    pass

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
    "job_name": "Automated Curriculum Learning Demo",         # set an arbitrary job name for identification
    "job_id": "Demo"                       # only relevant if "recipient" is set to a specific REST endpoint
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
        "score_threshold": 0.5            # agent must achieve an average score of this before
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
            "score_threshold": 0.5
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
        "score_threshold": 0.5
        }
    ]
}

# set up the Curriculum Strategy
configuration["parameters"]["production_strategy"] = {
    "name": "standard",
    "retain_inception": True,       # option to retain the inception from the Curriculum Phase
    # retain it here since the last Curriculum Objective is the same as
    # Production Objective. Previous top compounds will be relevant

    "number_of_steps": 3500,         # number of epochs to run the Production Phase
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
            "name": "Glide LigPrep Docking",
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
            },
            "weight": 1}]
    }
}

#configuration["parameters"]["scoring_function"] = scoring_function

configuration_JSON_path = os.path.join(output_dir, "CLwithDocking_config.json")
with open(configuration_JSON_path, 'w') as f:
    json.dump(configuration, f, indent=4)