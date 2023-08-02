.ONESHELL:
SHELL=/bin/zsh
PYTHON:=python3
PIP:=pip3
ENV_NAME:=rlmol
ENV_FILE:=env.yml
RUN:=main
RECEPTOR_PREP:=src.docker
RECEPTOR_PREP_PYTHON:=~/micromamba/envs/$(ENV_NAME)/bin/python2.7
ACTIVATE=source $($ MAMBA_ROOT_PREFIX)/etc/profile.d/micromamba.sh ; micromamba activate
docker=vina
receptor=""

.PHONY: list install uninstall create remove receptor test run clean

# Setup
default: list
install:
	curl micro.mamba.pm/install.sh | bash
	micromamba shell init
uninstall:
	rm -rf $(MAMBA_ROOT_PREFIX)
list:
	$(ACTIVATE) $(ENV_NAME) ; micromamba list
create: 
	micromamba create -f ${ENV_FILE} -y
remove:
	micromamba remove --name $(ENV_NAME) --all -y
	rm -rf $(MAMBA_ROOT_PREFIX)/envs/$(ENV_NAME)

# Run
receptor:
	$(ACTIVATE) $(ENV_NAME) ; $(RECEPTOR_PREP_PYTHON) -m $(RECEPTOR_PREP)
test:
	$(ACTIVATE) $(ENV_NAME) ; pytest test
run:
	$(ACTIVATE) $(ENV_NAME) ; $(PYTHON) $(RUN).py -d $(docker) -r $(receptor)
clean:
	rm -f *mol *mol2 *pdbqt *dlg *dpf summary_of_results* *map* *gpf *gpf_log 
test_docking:
	$(ACTIVATE) $(ENV_NAME) ; pytest test/test_docking.py
test_chemistry:
	$(ACTIVATE) $(ENV_NAME) ; pytest test/test_chemistry.py
test_environment:
	$(ACTIVATE) $(ENV_NAME) ; pytest test/test_environment.py

