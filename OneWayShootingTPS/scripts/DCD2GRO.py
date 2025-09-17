import numpy as np
import mdtraj as md
import os


burst = 1e5

folder = 'data_mu175_chain1_batch_4'
savefolder = 'cages_mu175_chain1_batch4'

gro_folder = 'gro_files'

os.makedirs(savefolder, exist_ok = True)

for run_folder in os.listdir(folder):
    run_folder_path = os.path.join(folder, run_folder)

    # Skip non-directories and hidden folders
    if not os.path.isdir(run_folder_path) or run_folder.startswith('.'):
        continue

    # Process trajectory files
    for succ_file in os.listdir(run_folder_path):
        if succ_file.startswith('traj_') and not succ_file.endswith('_0.dcd'):

            idx = succ_file.split('_')[1].split('.')[0]
            succ_file_path = os.path.join(run_folder_path, succ_file)
            # Input files

            mcg_file = os.path.join(run_folder_path, f'cv_{idx}.txt')

            mcg_frames, mcg_values = np.loadtxt(mcg_file, skiprows=1, unpack=True)


            cutoff_mask = mcg_values >= 300
            if np.any(cutoff_mask):
                first_cutoff_idx = np.argmax(cutoff_mask)
                B_cutoff = int( mcg_frames[first_cutoff_idx] / burst )
            else:
                B_cutoff = None  # or -1 or np.nan if you prefer a sentinel value

            print(B_cutoff)
            topology_file = os.path.join(gro_folder, "conf.gro")


            # Load the trajectory with a stride of n (e.g., every 10th frame)
            n = 1  # Change this to the desired interval
            traj = md.load(succ_file_path, top=topology_file, stride=n)
            original_frame_indices = np.arange(traj.n_frames) * n

            traj.time = np.arange(0, traj.n_frames) * burst * n

            if len(mcg_frames) != len(traj):
                print(error)
                continue

            traj = traj[: B_cutoff + 1]

            stride = 10

            strided_traj = traj[::stride]

            # --- Rename 'HOH' residues to 'SOL' ---
            for residue in strided_traj.topology.residues:
                if residue.name == 'HOH':
                    residue.name = 'SOL'

            print(f"Trajectory loaded using every {n}th step")
            savesubfolder = f'{savefolder}/{run_folder}'
            os.makedirs(savesubfolder, exist_ok = True)

            strided_traj.save(f'{savesubfolder}/traj_{idx}_stride_{stride}.gro')

