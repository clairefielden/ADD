> **How to run this notebook (terminal)?**
1. Ensure you are using a GPU in your Runtime Evironment on Google Colab Pro
2. Clone the github repo with all the necessary files: <br>
`git clone https://clairefielden/ADD.git` <br>
`cd ADD`
3. Install Miniconda
4. Activate the environments: <br>
`chmod +x ./scripts/install.sh` <br>
`make` <br>
`exit` <br>
`activate reinvent.v3.2`
5. Compile the configuration files: <br>
`cd ADD`
`python3 ./src/name_of_configuration_file.json`
6. Train the agent and perform docking simulation <br>
`python3 ./Reinvent/input.py ./Dockstream_CL/name_of_configuration_file.json`
7. Put results in a tar file and download it
`tar cvzf name_of_file.tar.gz Dockstream_CL`<br>
`from google.colab import files`<br>
`files.download('./ADD/name_of_file.tar.gz')`<br>
8. Start a tensorboard session on a terminal on your personal computer
`tar -xvzf name_of_file.tar.gz`<br>
`cd Dockstream_CL`<br>
`tensorboard --logdir progress.log`