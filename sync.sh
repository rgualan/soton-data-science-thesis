#!/bin/bash

rsync -av --info=flist model preprocessing util data pbs_scripts rmgs1u16@lyceum2.soton.ac.uk:/home/rmgs1u16/aq-tngapms/ 

#--delete
# The delete option also deletes the outputs of the pbs jobs

