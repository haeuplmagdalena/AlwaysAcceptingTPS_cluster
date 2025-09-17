import numpy as np
import os
import matplotlib.pyplot as plt
import mcg_gro_nl_torch
import mdtraj as md
import argparse

parser = argparse.ArgumentParser(description='folder inputs')
parser.add_argument('--input_file', dest='input_file', required = True)
parser.add_argument('--n_skip', dest='n_skip', required = True, type = int)
parser.add_argument('--output_folder', dest='output_folder')
parser.add_argument('--gro_folder', dest='gro_folder', default='/leonardo/home/userexternal/mhaeupl0/gro_files')
args = parser.parse_args()

frames2us =  2e-9

if args.output_folder is None:
    args.output_folder = os.path.dirname(os.path.abspath(args.input_file))

# Input files
xtc_file = args.input_file
topology_file = os.path.join(args.gro_folder, "conf.gro")

# Load the trajectory with a stride of n (e.g., every 10th frame)
n = args.n_skip  
traj = md.load(xtc_file, stride=n)
original_frame_indices = np.arange(traj.n_frames) * n

# Save the reduced trajectory as a .gro file
# Output file
# gro_file = "trajectory.gro"
# traj.save_gro(gro_file)

print(f"Trajectory loaded using every {n}th step")


outfile_path = os.path.join(args.output_folder, f"mcg_results_skip_{n}")

times_ps, CAR_COM_frames, WATER_COM_frames, box_lengths = mcg_gro_nl_torch.COM_calculation_mdtraj(traj)

MCG_OP_results = []

with open(outfile_path, 'w') as outFile:
    # Write the header
    outFile.write(f"{'OG frame idx':<15}{'Frame':<10}{'Time (ps)':<20}{'MCG_OP':<15}\n")

    for idx, time_ps in enumerate(times_ps):
        CAR_COM = CAR_COM_frames[idx]
        WATER_COM = WATER_COM_frames[idx]
        box_length = box_lengths[idx]
        time_us = time_ps * 1e-6

        frame = round(time_us / frames2us)

        # Assuming `MCG` takes these arguments and computes the order parameter
        MCG_OP, _, monomer_coords = mcg_gro_nl_torch.MCG_optimized(CAR_COM, WATER_COM, box_length, order=3, guest_cutoff=0.9)
        #np.savetxt( '../mcg_tests/monomers.xyz', monomer_coords)
        print(f'The largest cluster size in the system at frame {frame} and time {round(time_ps, 1)}ps is: {MCG_OP}')
        MCG_OP_results.append(MCG_OP)
        outFile.write(f"{original_frame_indices[idx]:<15}{frame:<10}{round(time_ps, 10):<20}{MCG_OP:<15.6f}\n")
