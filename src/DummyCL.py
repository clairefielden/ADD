# load dependencies
import os
import re
import json
import tempfile

reinvent_dir = os.path.expanduser("./Reinvent")
reinvent_env = os.path.expanduser("/usr/local/envs/reinvent.v3.2")
output_dir = os.path.expanduser("./Dockstream_CL")

try:
    os.mkdir(output_dir)
except FileExistsError:
    pass

configuration = {
    "version": 3,                                 # we are going to use REINVENT's newest release
    "run_type": "curriculum_learning",            # other run types: "sampling", "reinforcement_learning",
    #                  "transfer_learning",
    #                  "scoring" and "create_model"
    "model_type": "default"
}

configuration["logging"] = {
    "sender": "http://0.0.0.1",            # only relevant if "recipient" is set to "remote"
    "recipient": "local",                  # either to local logging or use a remote REST-interface
    "logging_frequency": 100,              # log every x-th steps
    "logging_path": os.path.join(output_dir, "progress.log"), # load this folder in tensorboard
    "result_folder": os.path.join(output_dir, "results"),     # will hold the compounds (SMILES) and summaries
    "job_name": "Automated Curriculum Learning Demo",         # set an arbitrary job name for identification
    "job_id": "Demo"                       # only relevant if "recipient" is set to a specific REST endpoint
}

configuration["parameters"] = {}

# First add the paths to the Prior, Agent, and set the curriculum type to automated
configuration["parameters"]["prior"] = "./models/random.prior.new"
configuration["parameters"]["agent"] = "./models/random.prior.new"
configuration["parameters"]["curriculum_type"] = "automated"

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
        # 1st scoring function below
        "scoring_function": {
            "name": "custom_product",     # this is our default one (alternative: "custom_sum")
            "parallel": False,
            "parameters": [{
                "component_type": "matching_substructure",     # enforce the match to a given substructure
                "name": "Pyrimidine",     # arbitrary name for the component
                "specific_parameters": {
                    "smiles": [
                        "[c]1[c][c]n[c]n1"     # a match with this substructure is required
                    ]
                },
                "weight": 1}]             # the weight of the component (default: 1)
        },
        "score_threshold": 0.8            # agent must achieve an average score of this before
        # progressing to the next Curriculum Objective
    },
        # 2nd scoring function below
        {
            "scoring_function": {
                "name": "custom_product",
                "parallel": False,
                "parameters": [{
                    "component_type": "matching_substructure",
                    "name": "H-Bonding Ring",
                    "specific_parameters": {
                        "smiles": [
                            "[c]1[c][c]nc(n1)[N]"
                        ]
                    },
                    "weight": 1}]
            },
            "score_threshold": 0.8
        },
        # 3rd scoring function below
        {
            "scoring_function": {
                "name": "custom_product",
                "parallel": False,
                "parameters": [{
                    "component_type": "matching_substructure",
                    "name": "H-Bonding Ring with Phenyl",
                    "specific_parameters": {
                        "smiles": [
                            "[c]1[c][c]c([c][c]1)[N]c2n[c][c][c]n2"
                        ]
                    },
                    "weight": 1}]
            },
            "score_threshold": 0.8
        },
        # 4th scoring function below
        {
            "scoring_function": {
                "name": "custom_product",
                "parallel": False,
                "parameters": [{
                    "component_type": "matching_substructure",
                    "name": "Double Ring",
                    "specific_parameters": {
                        "smiles": [
                            "[c]1[c][c]c([c][c]1)[N]c2n[c]c3c(n2)-[c][c][C][C]3"
                        ]
                    },
                    "weight": 1}]
            },
            "score_threshold": 0.8
        },
        # 5th scoring function below
        {
            "scoring_function": {
                "name": "custom_product",
                "parallel": False,
                "parameters": [{
                    "component_type": "matching_substructure",
                    "name": "Triple Ring",
                    "specific_parameters": {
                        "smiles": [
                            "[c]1[c][c]c([c][c]1)[N]c2n[c]c3c(n2)-c4c([c]n[n]4)[C][C]3"
                        ]
                    },
                    "weight": 1}]
            },
            "score_threshold": 0.8
        },
        # 6th scoring function below
        {
            "scoring_function": {
                "name": "custom_product",
                "parallel": False,
                "parameters": [{
                    "component_type": "matching_substructure",
                    "name": "Full Substructure",
                    "specific_parameters": {
                        "smiles": [
                            "[*]NC(=O)c1nn([*])c2c1CCc3cnc(Nc4ccccc4)nc23"
                        ]
                    },
                    "weight": 1}]
            },
            "score_threshold": 0.8
        },
    ]
}

# set up the Curriculum Strategy
configuration["parameters"]["production_strategy"] = {
    "name": "standard",
    "retain_inception": True,       # option to retain the inception from the Curriculum Phase
    # retain it here since the last Curriculum Objective is the same as
    # Production Objective. Previous top compounds will be relevant

    "number_of_steps": 100,         # number of epochs to run the Production Phase
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
            "component_type": "matching_substructure",
            "name": "Full Substructure",
            "specific_parameters": {
                "smiles": [
                    "[*]NC(=O)c1nn([*])c2c1CCc3cnc(Nc4ccccc4)nc23"
                ]
            },
            "weight": 1}]
    }
}

configuration_JSON_path = os.path.join(output_dir, "AutoCL_config.json")
with open(configuration_JSON_path, 'w') as f:
    json.dump(configuration, f, indent=4)
