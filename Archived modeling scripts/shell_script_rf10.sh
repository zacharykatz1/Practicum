#!/bin/bash
#$ -cwd -l mem=40g,time=10:00:00 -S /bin/bash -N JOBrf10 -j y -t 1-400 -M zacharykatz1@gmail.com -m aes

module load R/4.2.2

# R=/nfs/apps/R/4.2.2/bin/R
export R_LIBS_USER=/ifs/scratch/msph/software/R/library422:/ifs/home/msph/biostat/zak2132/R_LIB:/ifs/scratch/msph/software/R/library:$R_LIBS_USER

currind=$SGE_TASK_ID

Rscript rf_models_scenario10.R $currind
~
