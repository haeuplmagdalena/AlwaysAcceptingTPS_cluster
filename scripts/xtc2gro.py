import mdtraj as md
import argparse
import os

parser = argparse.ArgumentParser(description='folder inputs')
parser.add_argument("--input_file", type=str, required=True)
parser.add_argument("--n_skip")
args = parser.parse_args()

xtc_file = args.input_file
folder = os.path.dirname(os.path.abspath(xtc_file))
topology_file = os.path.join(folder, "conf.gro")

# Load the trajectory with a stride of n (e.g., every 10th frame)
n = args.n_skip  # Change this to the desired interval
traj = md.load(xtc_file, top=topology_file, stride=n)

for residue in traj.topology.residues:
    if residue.name == "HOH":
        residue.name = "SOL"

# Save the reduced trajectory as a .gro file
# Output file

gro_file = os.path.join(folder, f"trajectory_skip_{n}.gro")
traj.save_gro(gro_file)

print(f"Trajectory loaded using every {n}th step")
