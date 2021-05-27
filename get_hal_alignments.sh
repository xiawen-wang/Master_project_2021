# 10 jobs

#!/bin/sh

#$ -N COVID_hal  #name
#$ -cwd  #current dir
#$ -V
#$ -l h_rt=24:00:00   #runtime
#$ -l h_vmem=4G


# Module load
module load igmm/apps/R/3.3.0 igmm/apps/BEDTools/2.25.0 igmm/libs/ncurses/6.0 igmm/apps/samtools/1.3 igmm/apps/bcftools/1.3 igmm/apps/vcftools/0.1.13 igmm/libs/ensembl_api/86 igmm/apps/last/847 igmm/compilers/gcc/5.5.0

echo python get_hal_alignments.py $SGE_TASK_ID

python get_hal_alignments.py $SGE_TASK_ID
