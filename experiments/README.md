IMPROVING SIMULATED MOLECULAR DOCKING USING CURRICULUM LEARNING
=================================================================================================================

Installation
-------------

1. Install [Conda](https://conda.io/projects/conda/en/latest/index.html)
2. Clone this Git repository
   
        $ git clone https://clairefielden/ADD.git

2. Install [Schrodinger](https://www.schrodinger.com/)
   
        $ module load chpc/schrodinger/2019-4

3. Create the required environments:
**This will take a long time if you do not have access to multiple cores/GPUs**
   
        $ conda env create -f ./ADD/DockStream/environment.yml
        $ conda env create -f ./ADD/Reinvent/reinvent.yml

4. Further instructions on how to run these experiments are provided under (#usage)
   
        $ conda activate reinvent.v3.2

5. Run the curriculum learning configuration script:
   
        $ python ./ADD/Reinvent/input.py ./ADD/experiments/Agent4/curriculum.json


1. Ensure you are using a GPU in your Runtime Evironment
2. Clone the github repo with all the necessary files: <br>
`git clone https://clairefielden/ADD.git` <br>
`cd ADD`
3. Install [conda](https://docs.anaconda.com/free/anaconda/install/index.html): <br>
`module load chpc/python/anaconda/3-2021.11`
4. Install [Schrodinger](https://www.schrodinger.com/): <br>
`module load chpc/schrodinger/2019-4`
5. Activate the environments: <br>
`source /apps/chpc/anaconda/3-2021.11/etc/profile.d/conda.sh`<br>
`source activate /home/<name>/.conda/envs/reinvent.v3.2`
5. Compile the configuration files: <br>
- Compiling json files for docking: [DockStream](https://github.com/MolecularAI/DockStreamCommunity/blob/master/notebooks/demo_Glide.ipynb)
- Compiling json files for reinforcement learning: [REINVENT3.2](https://github.com/MolecularAI/ReinventCommunity/blob/master/notebooks/Reinforcement_Learning_Demo.ipynb)
- The configuration files for these experiments are provided in the [<experiments>]() subdirectory of this repository.

     
System Requirements
-------------------

* Python 3.7
* Cuda-enabled GPU
* `REINVENT` has been tested on Linux


Tutorials / `jupyter` notebooks
-----
There is another repository containing useful `jupyter` notebooks related to `REINVENT` 
called [ReinventCommunity](https://github.com/MolecularAI/ReinventCommunity). Note, that it uses a
different `conda` environment to execute, so you have to set up a separate environment.


Usage
-----

For concrete examples, you can check out the Jupyter notebook examples in the ReinventCommunity repo.
Running each example will result in a template file.There are templates for many running modes. 
Each running mode can be executed by `python input.py some_running_mode.json` after activating the environment.
    
   Templates can be manually edited before using. The only thing that needs modification for a standard run are the file 
   and folder paths. Most running modes produce logs that can be monitored by `tensorboard`.


Tests 
-----
The REINVENT project uses the `unittest` framework for its tests; before you run them you first have to create a 
configuration, which the tests will use. In the project directory, create a `config.json` file in the `configs/` directory; 
you can use the example config (`example.config.json`) as a base. The simplest way is to make a copy of the `example.config.json`
and name it `config.json`. At this point, `REINVENT` can be executed. If you want to further run the unit tests, relevant paths 
will need to be specified in your 'config.json' file, e.g, if testing reinforcement learning, the corresponding unit tests 
will require a prior: specify the path to the prior in the "PRIOR_PATH" field.

Important: Make sure that you set `MAIN_TEST_PATH` to a non-existent directory; it is where temporary files will be
written during the tests; if it is set to an existing directory, that directory will be removed once the tests have finished.

Some tests require a proprietary OpenEye license; you have to set up a few things to make the tests read your
license.  The simple way is to just set the `OE_LICENSE` environment variable to the path of the file containing the
license.  If you just want to set the license in the `reinvent_scoring` Conda environment, it is a bit more complicated,
but you only have to do it once.

```
(reinvent-scoring) $ cd $CONDA_PREFIX
$ mkdir -p etc/conda/activate.d
$ mkdir -p etc/conda/deactivate.d
```

Put the following in `etc/conda/activate.d/env_vars.sh`.

```
#!/bin/sh
export OE_LICENSE='</path/to/your/oe_license/file>'
```

And put the following in `etc/conda/deactivate.d/env_vars.sh`.

```
#!/bin/sh
unset OE_LICENSE
```

Once you have created the files, deactivate and re-activate the environment, and `echo $OE_LICENSE` should output the
path to the license file.
Once you have a configuration and your license can be read, you can run the tests.

```
$ python main_test.py
```
Liscense 
-----

MIT © [Open Sauced](LICENSE)