#!/bin/bash
#SBATCH --job-name=TPS180-1          # Job name (takes the name from sbatch --job-name)
#SBATCH --output=logs/tps.out       # Log file (named after the job)
#SBATCH --partition=boost_usr_prod # Partition
##SBATCH --time=1-00:00:00            # Maximum runtime
#SBATCH --qos=boost_qos_lprod  # quality of service
#SBATCH -N 1                   # 1 node
#SBATCH --time 4-00:00:00      # format: D-HH:MM:SS
#SBATCH --gres=gpu:1               # Request 1 GPU
#SBATCH --ntasks-per-node=1    # 1 tasks out of 32
#SBATCH --mem=128000                  # Memory per node
#SBATCH --mail-type=ALL            # Email notifications

module load cuda/12.1
module load hdf5
source /leonardo/home/userexternal/mhaeupl0/clathrate_env/bin/activate

python Clathrate_TPS.py
