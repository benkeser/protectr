#!/bin/bash

PARTITION=$1
SETTING=$2

module purge
module load R/4.4.0

export PATH=/apps/R/4.4.0/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin
 
sbatch --partition=$PARTITION \
  --nodes=1 \
  --ntasks-per-node=1 \
  --cpus-per-task=31 \
  --mem-per-cpu=6G \
  --job-name=tb_.%a \
  --output=/projects/dbenkes/allison/protectr/scratch/%a_%J.out \
  --export=SETTING=$SETTING \
  --wrap "/apps/R/4.4.0/bin/Rscript run_simulation.R"
