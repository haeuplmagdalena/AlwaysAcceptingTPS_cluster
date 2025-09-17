#!/bin/bash
#SBATCH --account=L-AUT_Meier
#SBATCH --job-name=grade
#SBATCH --output=logs/tps_%j.out
#SBATCH --error=logs/tps_%j.err
#SBATCH --partition=boost_usr_prod
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --gres=gpu:4
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=8
#SBATCH --time=0-12:00:00
##SBATCH --mem=128000
#SBATCH --mail-type=ALL

module load cuda/12.1
module load hdf5
source ../../clathrate_env/bin/activate

python start_simulation.py --N 200000 --N_out 1 --input_file /leonardo_scratch/large/userexternal/mhaeupl0/InitialTrajectories/Crystalline/cv_10/N_out_1_N_2e5/frame_cv_10.dcd
