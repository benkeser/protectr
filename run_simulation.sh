#!/bin/bash

PARTITION=$1
SETTING=$2

module purge
module load R/4.4.0

export PATH=/apps/R/4.4.0/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin
 
sbatch --partition=$PARTITION \
  --nodes=1 \
  --ntasks-per-node=1 \
  --cpus-per-task=32 \
  --mem-per-cpu=6G \
  --job-name=tb_%a \
  --export=SETTING=$SETTING \
  --wrap "/apps/R/4.4.0/bin/Rscript run_simulation.R"