#!/bin/sh
#This shell is called by "Run_dock_simulations.sh" at the run docker phase.
#It builds the neurodamus configuration, runs the simulations, and prepare the environment to plot the simulation results.
cd /mnt/mydata/
cp $NEURODAMUS_MODS_DIR/* MainOlfactoryBulb/mod/
build_neurodamus.sh MainOlfactoryBulb/mod/
export HOC_LIBRARY_PATH=/mnt/mydata/MainOlfactoryBulb/mod:$HOC_LIBRARY_PATH
./x86_64/special -mpi -python $NEURODAMUS_PYTHON/init.py --configFile=simulations/run_no_dynamics/simulation_MC_config.json
apt-get update
yes | apt-get install vim
pip install bluepysnap

