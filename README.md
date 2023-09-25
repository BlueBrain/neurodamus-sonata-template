# neurodamus-sonata-template
Template with minimal SONATA circuit with 1 neuron morphology and emodel to be tested with neurodamus docker image


Hippocampus use case:

1. clone this repo
2. go to to_docker/manipulate_nodes directory

3. to undertand the contents of nodes_template.h5 run:

`h5ls -rd nodes_correct.h5`

which will give the tree of hdf5 file with the values we want to obtain from nodes_to_change.h5. We want to modify
model_template and morphology for node population 'popA/0' 

```
/nodes/popA/0/model_template Dataset {1}
    Data:
         "hoc:CA1_pyr_cACpyr_mpg141017_a1_2_idC_2019032814340"
...
/nodes/popA/0/morphology Dataset {1}
    Data:
         "dend-050921AM2_axon-mpg141017_a1-2_idC"
```

4. Then we can specify which morphology and emodel to overwrite using the script. First lets copy the template.

```
cp nodes_to_change.h5 nodes_modified.h5
```

Now we will use update_nodes.py to modify nodes_modified.h5's contents

The first argument is to give node file to edit, then we give one by one the dataset field to change and the value we want it to be.

`python update_nodes.py {nodes_file} {hdf5_path_to_dataset} {new_value} {hdf5_path_to_dataset2} {new_value2}`

which translates into:

`python update_nodes.py nodes_modified.h5 /nodes/popA/0/model_template hoc:CA1_pyr_cACpyr_mpg141017_a1_2_idC_2019032814340 /nodes/popA/0/morphology dend-050921AM2_axon-mpg141017_a1-2_idC`

Use this nodes.h5 file in circuit directory. Better to copy it there.

`cp nodes_modified.h5 ../circuit/nodes.h5`

or you can make a symbolic link in order not to duplicate files

5. Now we follow the instructions in neurodamus to pull its docker image, compile our mod files and run the simulation

`docker pull bluebrain/neurodamus`

will give us the docker image

you can check if it is there with `docker image ls`

6. Then we mount our files into docker which we put before in to_docker/

`docker run --rm -it -v full/path/to/to_docker:/mnt/mydata bluebrain/neurodamus`

This will instantiate the docker and your terminal will change to dockers. You can always exit the image with `exit` command.

7. Compile your mod files into neurodamus. Here we use hippocampus mod files and also common mod files for all circuits which I put under to_docker/ for you.

```
cd /mnt/mydata/
cp $NEURODAMUS_MODS_DIR/* hippocampus/mod/
build_neurodamus.sh hippocampus/mod/
```

this will create a new folder x86_64/ under mnt/mydata. Then we run:

`export HOC_LIBRARY_PATH=/mnt/mydata/hippocampus/mod:$HOC_LIBRARY_PATH`

8. Once it compiles the mod files and we export its path to HOC_LIBRARY_PATH, we can run the simulation from mnt/mydata/ directory:

`/x86_64/special -mpi -python $NEURODAMUS_PYTHON/init.py --configFile=simulations/run_no_dynamics/simulation_config.json`


Then results will be under simulations/run_no_dynamics


Useful Links:
1. Neurodamus Docker: https://github.com/BlueBrain/neurodamus/tree/main/docker
2. BBP mod files : https://bbpgitlab.epfl.ch/hpc/sim/models

