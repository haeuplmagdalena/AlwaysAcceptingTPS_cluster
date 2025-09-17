#!/bin/bash
#SBATCH --account=L-AUT_017
#SBATCH --job-name=grade_AA_01
#SBATCH --output=logs/grade_%j.out
#SBATCH --error=logs/grade_%j.err
#SBATCH --partition=boost_usr_prod
##SBATCH --qos=boost_qos_bprod
#SBATCH --nodes=1
##SBATCH --exclusive
##SBATCH --gres=gpu:4
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
##SBATCH --mem=128000
#SBATCH --mail-type=ALL

module load cuda/12.1
module load hdf5
source ../../../clathrate_env/bin/activate

python -u ../scripts/GradeAnalysis.py
