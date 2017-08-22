#!/bin/bash
# run simple Fluent job in parallel on 1 node

# set default resource requirements for job (4 processors on 1 node for 4 hours)
# - these can be overridden on the qsub command line
#PBS -l nodes=1:ppn=4
#PBS -l walltime=04:00:00

#Change to directory from which job was submitted
cd ~/aq-tngapms/

# set number of processors to run on (using list of node names in file $PBS_NODEFILE)
nprocs=`wc -l $PBS_NODEFILE | awk '{ print $1 }'`
#nprocs=4

# load fluent module so that we find the fluent command
module load R

# Run default version of Fluent in 3d mode in parallel over $nprocs processors
# Fluent commands are in file elbow.jou, output messages go to output_file
#fluent 3d -t$nprocs -g -i elbow.jou > output_file
R --vanilla < model/spacetime/spacetime_epa_ca_ozone_v4.R &> logs/st.log
