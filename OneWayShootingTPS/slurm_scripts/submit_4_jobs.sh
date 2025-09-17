#!/bin/bash
##SBATCH --job-name=175-2
#SBATCH --account=L-AUT_017
#SBATCH --output=logs/tps_%j.out
#SBATCH --error=logs/tps_%j.err
#SBATCH --partition=boost_usr_prod
#SBATCH --qos=boost_qos_lprod #also change time
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --gres=gpu:4
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=8
#SBATCH --time=4-00:00:00
##SBATCH --mem=128000
#SBATCH --mail-type=ALL

module load cuda/12.1
module load hdf5
source ../../clathrate_env/bin/activate

# Start 4 jobs in parallel using srun with exclusive CPU/GPU binding
srun --gpus-per-task=1 python Clathrate_TPS.py & 

srun --gpus-per-task=1 python Clathrate_TPS.py &

srun --gpus-per-task=1 python Clathrate_TPS.py &

srun --gpus-per-task=1 python Clathrate_TPS.py & 

wait
