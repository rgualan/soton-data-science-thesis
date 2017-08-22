#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=12:00:00
cd ~/aq-tngapms/
module load R
R --vanilla < preprocessing/impute_ozone.R &> logs/ozone_imputation.log

