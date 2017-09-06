#!/bin/bash
#PBS -l nodes=1:ppn=6
#PBS -l walltime=12:00:00
cd ~/aq-tngapms/
module load R
R --vanilla < model/spacetime/spacetime_epa_ca_ozone.R &> logs/st.log

