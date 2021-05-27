#!/bin/sh

#$ -N VEP_C2_eu  #name
#$ -cwd  #current dir
#$ -V
#$ -l h_rt=24:00:00   #runtime
#$ -l h_vmem=8G


# Module load
module load igmm/apps/vep/97

vep --database -i vep_input.txt -o VEP_C2_eu.txt
