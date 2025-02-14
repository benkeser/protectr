#!/bin/bash

PARTITION=$1
SETTING=$2

module purge
module load R/4.4.0

export PATH=/apps/R/4.4.0/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin
 
sbatch 
  --array=1
  --partition=$PARTITION \
  --cpus-per-task=1 \
  --mem-per-cpu=10G \
  --job-name=tb_hiv \
  --output=/projects/dbenkes/allison/protectr/scratch/${SETTING}_boot_%J.out \
  --export=SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID,SETTING=$SETTING \
  --wrap "/apps/R/4.4.0/bin/Rscript run_bootstrap.R"
