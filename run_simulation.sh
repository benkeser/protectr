#!/bin/bash

PARTITION=$1
SETTING=$2

module purge
module load R/4.4.0

export PATH=/apps/R/4.4.0/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin
 
sbatch --partition=$PARTITION \
  --nodes=1 \
  --ntasks-per-node=1 \
  --cpus-per-task=4 \
  --mem-per-cpu=46G \
  --job-name=tb_hiv \
  --output=/projects/dbenkes/allison/protectr/scratch/${SETTING}_%J.out \
  --export=SETTING=$SETTING \
  --wrap "/apps/R/4.4.0/bin/Rscript run_simulation.R"
