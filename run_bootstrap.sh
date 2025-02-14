#!/bin/bash

PARTITION=$1
SETTING=$2
NSEEDS=$3

module purge
module load R/4.4.0

export PATH=/apps/R/4.4.0/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin
 
sbatch --array=1-$NSEEDS \
  --partition=$PARTITION \
  --job-name=tb_hiv \
  --output=/projects/dbenkes/allison/protectr/scratch/${SETTING}_boot_${SLURM_ARRAY_TASK_ID}_%J.out \
  --export=SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID,SETTING=$SETTING \
  --wrap "/apps/R/4.4.0/bin/Rscript run_bootstrap.R"
