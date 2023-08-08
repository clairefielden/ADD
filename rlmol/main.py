from src import agent, environment, docking
import argparse
import pathlib
import logging
import yaml
import csv
import pathlib

"""logger = logging.getLogger(__name__)
logging.basicConfig(format='%(levelname)s : %(module)s : %(asctime)s : %(message)s',
                    datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.DEBUG)
logger.info('Started')

parser = argparse.ArgumentParser(
    prog='rlmol',
    description='Uses RL to explore drug space'
)
parser.add_argument('-d', '--docker', help='the docking program to use',
                    choices=['vina', 'gpu', 'cpu'], default=42)
parser.add_argument('-r', '--receptor', help='path to the receptors pdb file')
args = parser.parse_args()

receptor_path = pathlib.Path(args.receptor)
if not receptor_path.exists():
    print("Receptor pdb file does not exist")
    raise SystemExit(1)

with open('etc/config.yml', 'r') as file:
    config = yaml.safe_load(file)

if args.docker == 'vina':
     v = docking.get_vina_instance()
     v.set_receptor(receptor_path)
     e = environment.initialise_vina(v)
     a = agent.initialise(e)
     a.train()

target_receptor = 'pfpi4k-no-ligand'
docking.prepare_receptor(target_receptor)
grid_params = docking.read_grid_file('pfpi4k-docking-box.txt')

napthyridine = 'C1=CC2=C(C=CC=N2)N=C1'
aminopyridine = 'C1=CC=NC(=C1)N'
imidazopyridazine = 'C1=CC2=NC=CN2N=C1'
"""
#e = environment.MoleculeEnvironment(target=target_receptor,
                #grid=grid_params,
                #nruns=10,
                #starting_mol=napthyridine)


files = [f for f in pathlib.Path().glob("./CL_Scaffolds/*.csv")]
print(len(files))

my_file =files[len(files)-1]
scaffolds = []

with open(my_file, 'r') as file:
     csvreader = csv.reader(file)
     for row in csvreader:
        scaffolds.append(row[1])
     

for scaffold in range(0,10):
        with open('./config.yml') as f:
            doc = yaml.safe_load(f)

        doc['training']['starting_molecule'] = scaffolds.pop()
        print(doc['training']['starting_molecule'])

        with open('./config.yml', 'w') as f:
            yaml.safe_dump(doc, f)

        with open('./config.yml', 'r') as file:
            config = yaml.safe_load(file)

        e = environment.AutodockEnvironment(config)

        a = agent.Agent(e)

        a.train()

