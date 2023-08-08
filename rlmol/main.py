from src import agent, environment, docking
import argparse
import pathlib
import logging
import yaml
import csv
import pathlib

files = [f for f in pathlib.Path().glob("./CL_Results/results/*.csv")]
print(len(files))

my_file =files[len(files)-1]
scaffolds = []

with open(my_file, 'r') as file:
     csvreader = csv.reader(file)
     for row in csvreader:
        scaffolds.append(row[1])

for scaffold in range(1,11):
        with open('./rlmol/config.yml') as f:
            doc = yaml.safe_load(f)

        doc['training']['starting_molecule'] = scaffolds[scaffold]
        print(doc['training']['starting_molecule'])

        with open('./rlmol/config.yml', 'w') as f:
            yaml.safe_dump(doc, f)

        with open('./rlmol/config.yml', 'r') as file:
            config = yaml.safe_load(file)

        e = environment.AutodockEnvironment(config)

        a = agent.Agent(e)

        a.train()
