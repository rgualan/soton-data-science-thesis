#!/bin/bash
#PBS -l nodes=1:ppn=11
#PBS -l walltime=6:00:00
cd ~/aq-tngapms/
module load R
R --vanilla < model/RF/rf_epa_ca_parallel.R &> logs/rfp.log

