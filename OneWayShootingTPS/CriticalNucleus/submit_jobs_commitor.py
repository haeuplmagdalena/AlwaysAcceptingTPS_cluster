import subprocess
import numpy as np
from scipy.stats import norm
import mdtraj as md
import os
import shutil
import random


optimal_mcg_values = np.array([145, 155, 160, 165, 175])

base_input_path = '/leonardo_work/L-AUT_Coretti/CODE/ClathrateTPS'
traj_list_file = 'CrystallineTrajectoriesList.txt'
traj_list = traj_list = np.loadtxt(traj_list_file, dtype=str)


gro_folder = "../gro_files"

topology_file = os.path.join(gro_folder, "conf.gro")

os.makedirs("logs", exist_ok=True)

#frameidx = np.unique(og_file_idx[closest_indices])
#shooting_point = md.load_frame(xtc_file, index=frameidx[0], top=topology_file)


folder = "$FAST/CriticalNucleus/Crystalline"



for optimal_mcg in optimal_mcg_values:#closest_indices:

    print(optimal_mcg)

    random_numbers = random.sample(range(len(traj_list)), 14)


    for sim_num in range(len(random_numbers)):

        chosen_trajectory = traj_list[random_numbers[sim_num]]
        chosen_trajectory_folder = os.path.dirname(chosen_trajectory)

        print(chosen_trajectory)

        traj_index = chosen_trajectory.split('_')[-1].split('.')[0]


        ##move all files into this folder

        mcg_file = os.path.join(base_input_path, chosen_trajectory_folder, f"cv_{traj_index}.txt")
        
        mcg_frames, mcg_values = np.loadtxt(mcg_file, skiprows = 1, unpack = True)
        closest_indices = np.abs(mcg_values[:, None] - optimal_mcg).argmin(axis=0)
        mcg_value = int(mcg_values[closest_indices])

        shooting_point_xyz = md.load_frame(os.path.join(base_input_path, chosen_trajectory), index=closest_indices, top=topology_file)
        box_vectors = shooting_point_xyz.unitcell_vectors

        
        sp_dir = os.path.expandvars(f"{folder}/shooting_point_mcg_{optimal_mcg}")
        os.makedirs(sp_dir, exist_ok = True)
        print(f"created directory: {sp_dir}")
        base_path = os.path.expandvars(f"{folder}/shooting_point_mcg_{optimal_mcg}/sim_nr_{sim_num}")
        os.makedirs(base_path, exist_ok = True)

        with open(os.path.join(base_path, f'log.txt'), "a") as file:  # Open file inside the function
            row = f"Starting at MCG {mcg_value} from the trajectory loaded from path {os.path.join(base_input_path, chosen_trajectory)}"
            file.write(row + "\n")
            
        shutil.copy(os.path.join(base_input_path, chosen_trajectory), base_path)
        
        
        shooting_point_path = os.path.join(base_path, "shooting_point.xyz")
        shooting_point_xyz.save_xyz(shooting_point_path)
        np.save(os.path.join(base_path, "box_vectors.npy"), box_vectors)

        shutil.copyfile(os.path.expandvars("$HOME/gro_files/conf.gro"), os.path.expandvars(f"{folder}/shooting_point_mcg_{optimal_mcg}/sim_nr_{sim_num}/conf.gro"))
        shutil.copyfile(os.path.expandvars("$HOME/gro_files/topol.top"), os.path.expandvars(f"{folder}/shooting_point_mcg_{optimal_mcg}/sim_nr_{sim_num}/topol.top"))
        
        
        job_name = f"sim_{mcg_value}_{sim_num}"
        subprocess.run([
            "sbatch",
            "--job-name", job_name,
            "--gres=gpu:1",  # Request 1 GPU per job
            "--output", f"logs/{job_name}.out",  # Log file
            "--partition", "boost_usr_prod",  # Specify the developer/test queue (this depends on the actual name of the developer queue on the system)
            "--time", "23:59:59",  # Maximum runtime
            "--mem", "32G",  # Memory per node
            "--mail-type", "ALL",  # Email notifications
            "run_simulation.sh",  # SLURM script
            str(mcg_value),  # Pass index as argument
            str(sim_num)  # Pass simulation number
        ])
