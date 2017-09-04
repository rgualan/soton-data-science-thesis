#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=06:00:00
cd ~/aq-tngapms/
module load R
R --vanilla < preprocessing/feature_importance.R &> logs/importance.log

