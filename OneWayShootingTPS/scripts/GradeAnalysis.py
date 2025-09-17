import numpy as np
import mdtraj as md
import os
import subprocess
import fnmatch

SCRIPT="/leonardo_work/L-AUT_Coretti/AAAClathrates/scripts/GRADE_updated/GRADE"

savefolder = 'cages_mu175_chain1_batch4'
stride = 10

gro_folder = 'gro_files'

for run_folder in os.listdir(savefolder):
    run_folder_path = os.path.join(savefolder,run_folder)

    # Skip non-directories and hidden folders
    if not os.path.isdir(run_folder_path) or run_folder.startswith('.'):
        continue
    
    # Process trajectory files
    for gro_file in os.listdir(run_folder_path):
        pattern = f'traj_*_stride_{stride}*.gro'
        print(f'checking for pattern {pattern}')
        if fnmatch.fnmatch(gro_file, pattern):
            print(f"Matched: {gro_file}")
            gro_file_path = os.path.join(run_folder_path, gro_file)
            print(f'Run on {gro_file_path}')

            command = [SCRIPT, "-i", gro_file_path]
            
            try:
                result = subprocess.run(
                    command,
                    check=True,          # Raises error if command fails
                    capture_output=True, # Captures stdout/stderr
                    text=True           # Returns output as string (not bytes)
                )
                print(f"Output for {gro_file}:")
                print(result.stdout)
            except subprocess.CalledProcessError as e:
                print(f"Error processing {gro_file}:")
                print(e.stderr)

