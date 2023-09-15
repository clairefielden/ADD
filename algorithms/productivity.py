#
# calculation of productivity as described in:
#Average convergence rate of evolutionary algorithms in continuous optimization,
#Yu Chen and Jun He
# https://doi.org/10.1016/j.ins.2020.12.076
# https://www.sciencedirect.com/science/article/pii/S0020025520312421},

import pandas as pd

#the name of the file downloaded from tensorboard containing loss values
convergence_file = "ADD/results/benchmark/benchmark_nll.csv"
data = pd.read_csv(convergence_file, nrows=101)
x = data['Step'].tolist()
y = data['Value'].tolist()

total = 0
for i in y:
    total = total+i

acr = total/100

# the name of the file downloaded from tensorboard containing docking scores
convergence_file = "ADD/results/benchmark/benchmark_ds.csv"
data = pd.read_csv(convergence_file, nrows=101)
x = data['Step'].tolist()
y = data['Value'].tolist()

total = 0
for i in y:
    total = total + i

ds = total / 100

productivity = ds*acr
print(productivity)