#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=12:00:00
cd ~/aq-tngapms/
module load R
R --vanilla < model/ml/svr_gs.R &> logs/svr_gs.log

