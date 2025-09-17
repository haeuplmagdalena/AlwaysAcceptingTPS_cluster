#!/bin/bash
#SBATCH --job-name=a_sim_%x          # Job name (takes the name from sbatch --job-name)
#SBATCH --output=logs/%x.out       # Log file (named after the job)
#SBATCH --partition=boost_usr_prod # Partition
#SBATCH --time=23:59:59            # Maximum runtime
#SBATCH --gres=gpu:1               # Request 1 GPU
#SBATCH --mem=32G                  # Memory per node
#SBATCH --mail-type=ALL            # Email notifications

module load cuda/12.1
module load hdf5
source /leonardo/home/userexternal/mhaeupl0/clathrate_env/bin/activate

# Arguments passed from the Python script
INDEX=$1
SIM_NUM=$2

echo "Running simulation $SIM_NUM for index $INDEX"
python production_run_commitor.py --idx_mcg $INDEX --sim_num $SIM_NUM  --N 100000000 --P 500 --T 260
