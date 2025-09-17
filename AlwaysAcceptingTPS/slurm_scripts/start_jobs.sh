#!/bin/bash

for i in {1..3}
do
  sbatch --job-name="TPS_$i" submit_4_jobs.slurm
done

