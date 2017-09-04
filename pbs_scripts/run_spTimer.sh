#!/bin/bash
#PBS -l nodes=1:ppn=11
#PBS -l walltime=12:00:00
cd ~/aq-tngapms/
module load R
R --vanilla < model/spTimer/sptimer_epa_parallel_v2.R &> logs/spTimer.p.test.log

