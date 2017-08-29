#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l walltime=24:00:00
cd ~/aq-tngapms/
module load R
R --vanilla < model/RF/rf_epa_ca.R &> logs/rf.log

