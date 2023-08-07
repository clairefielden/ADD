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

# specify the target preparation JSON file as a dictionary and write it out
tp_dict = {
  "target_preparation":
  {
    "header": {                                   # general settings
      "logging": {                                # logging settings (e.g. which file to write to)
        "logfile": log_file_target_prep
      }
    },
    "input_path": apo_1UYD_path,                  # this should be an absolute path
    "fixer": {                                    # based on "PDBFixer"; tries to fix common problems with PDB files
      "enabled": True,
      "standardize": True,                        # enables standardization of residues
      "remove_heterogens": True,                  # remove hetero-entries
      "fix_missing_heavy_atoms": True,            # if possible, fix missing heavy atoms
      "fix_missing_hydrogens": True,              # add hydrogens, which are usually not present in PDB files
      "fix_missing_loops": False,                 # add missing loops; CAUTION: the result is usually not sufficient
      "add_water_box": False,                     # if you want to put the receptor into a box of water molecules
      "fixed_pdb_path": fixed_pdb_path            # if specified and not "None", the fixed PDB file will be stored here
    },
    "runs": [                                     # "runs" holds a list of backend runs; at least one is required
      {
        "backend": "AutoDockVina",                # one of the backends supported ("AutoDockVina", "OpenEye", ...)
        "output": {
          "receptor_path": adv_receptor_path      # the generated receptor file will be saved to this location
        },
        "parameters": {
          "pH": 7.4,                              # sets the protonation states (NOT used in Vina)
          "extract_box": {                        # in order to extract the coordinates of the pocket (see text)
            "reference_ligand_path": reference_ligand_path,   # path to the reference ligand
            "reference_ligand_format": "PDB"                  # format of the reference ligand
          }
}}]}}

with open(target_prep_path, 'w') as f:
    json.dump(tp_dict, f, indent="    ")
