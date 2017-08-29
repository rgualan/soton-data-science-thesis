#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l walltime=6:00:00
cd ~/aq-tngapms/model/spTimer/spAir_no2_best_model/
module load R
R --vanilla < no2_best_model.R &> simulation.log

