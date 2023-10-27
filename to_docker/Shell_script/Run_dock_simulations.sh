#!/bin/sh
#This script setups the configurations for the the simulation of the Mitral, middle Tufted, and Granule cells
#before running a docker to run the simulation and obtain a plot.
#Each cell type is ran indepedently in different dockers


#Mitral:
printf "Setting the simulation for a Single Mitral Cell"
cd /gpfs/bbp.cscs.ch/project/proj157/scratch/home/shared/SingleCellTesting/neurodamus-sonata-template/to_docker/manipulate_nodes
cp nodes_to_change.h5 nodes_modified.h5
python update_nodes.py nodes_modified.h5 /nodes/popA/0/model_template hoc:MOBmi_MC_bAC /nodes/popA/0/morphology MC0000

cp nodes_modified.h5 ../circuit/nodes.h5
printf "Starting Docker to run the simulation of a Single Mitral Cell"
docker run --name sing-cell-test --rm -v /gpfs/bbp.cscs.ch/project/proj157/scratch/home/shared/SingleCellTesting/neurodamus-sonata-template/to_docker:/mnt/mydata bluebrain/neurodamus -c './mnt/mydata/Shell_script/Run_in_docker_MC.sh && cd /mnt/mydata/simulations/ && python plot_sonata_soma_trace.py 0'
printf "\nMitral Cell is done, plot can be found in \n /gpfs/bbp.cscs.ch/project/proj157/scratch/home/shared/SingleCellTesting/neurodamus-sonata-template/to_docker/simulations/MC_soma_trace.png\n\n=============\n"


#Middle Tufted:
printf "Setting the simulation for a Single middle Tufted Cell"
python update_nodes.py nodes_modified.h5 /nodes/popA/0/model_template hoc:MOBopl_mTC_bSTUT /nodes/popA/0/morphology mTC0635

cp nodes_modified.h5 ../circuit/nodes.h5
printf "Starting Docker to run the simulation of a Single middle Tufted Cell"
docker run --name sing-cell-test --rm -v /gpfs/bbp.cscs.ch/project/proj157/scratch/home/shared/SingleCellTesting/neurodamus-sonata-template/to_docker:/mnt/mydata bluebrain/neurodamus -c './mnt/mydata/Shell_script/Run_in_docker_mTC.sh && cd /mnt/mydata/simulations/ && python plot_sonata_soma_trace.py 1'
printf "\nmiddle Tufted Cell is done, plot can be found in \n /gpfs/bbp.cscs.ch/project/proj157/scratch/home/shared/SingleCellTesting/neurodamus-sonata-template/to_docker/simulations/mTC_soma_trace.png\n\n=============\n"


#Granule:
printf "Setting the simulation for a Single Granule Cell"
python update_nodes.py nodes_modified.h5 /nodes/popA/0/model_template hoc:MOBgr_GC_dNAC /nodes/popA/0/morphology MOBgr_GC000

cp nodes_modified.h5 ../circuit/nodes.h5
printf "Starting Docker to run the simulation of a Single middle Tufted Cell"
docker run --name sing-cell-test --rm -v /gpfs/bbp.cscs.ch/project/proj157/scratch/home/shared/SingleCellTesting/neurodamus-sonata-template/to_docker:/mnt/mydata bluebrain/neurodamus -c './mnt/mydata/Shell_script/Run_in_docker_GC.sh && cd /mnt/mydata/simulations/ && python plot_sonata_soma_trace.py 2'
printf "\nGranule Cell is done, plot can be found in \n /gpfs/bbp.cscs.ch/project/proj157/scratch/home/shared/SingleCellTesting/neurodamus-sonata-template/to_docker/simulations/GC_soma_trace.png\n\n=============\n"

printf "All plots can be found in \n /gpfs/bbp.cscs.ch/project/proj157/scratch/home/shared/SingleCellTesting/neurodamus-sonata-template/to_docker/simulations/\n"
