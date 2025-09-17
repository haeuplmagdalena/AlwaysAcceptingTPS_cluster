#!/bin/bash
#SBATCH --job-name=dcd2gro
#SBATCH --output=logs/tps_%j.out
#SBATCH --error=logs/tps_%j.err
#SBATCH --partition=boost_usr_prod
#SBATCH --account=L-AUT_017
#SBATCH --qos=boost_qos_dbg
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --gres=gpu:4
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=8
##SBATCH --time=0-02:00:00
##SBATCH --mem=128000
#SBATCH --mail-type=ALL

module load cuda/12.1
module load hdf5
source ../../clathrate_env/bin/activate

python scripts/DCD2GRO.py
