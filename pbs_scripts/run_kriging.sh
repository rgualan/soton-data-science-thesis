#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l walltime=06:00:00
cd ~/aq-tngapms/
module load R
R --vanilla < model/spatial/kriging.R &> logs/sp.log

