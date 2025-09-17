import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
#from analysis_scripts import *
import argparse
import os
import scipy.constants as cs
import mcg_gro as mcg_gro
import pandas as pd
from matplotlib_config import colors
from scipy.fft import fft, fftfreq
import numpy as np

step2us = 2e-9

plt.rcParams.update({
    # 'font.family': 'sans-serif',
    # 'font.sans-serif': ['Arial', 'Helvetica', 'sans-serif'],
    'font.family': 'serif',  # Use a LaTeX-compatible serif font
    'font.serif': ['Computer Modern Roman'],  # Match LaTeX's default font
    'figure.figsize': (6.5, 4.5),
    'figure.dpi': 100,
    'text.usetex': True,
    'font.size': 18,
    'axes.labelsize': 20,
    'axes.titlesize': 16,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 16,
    'mathtext.default':'regular',
})


folder = "/leonardo_scratch/fast/L-AUT_Coretti/CommittorAnalysisData250"


data = {}

for sp_folder in os.listdir(folder):

    sp_nr = int(sp_folder.split("_")[-1])
    data[sp_nr] = {}
    sp_folder_path = os.path.join(folder, sp_folder)

    if os.path.isdir(sp_folder_path):
        for sim_nr_folder in os.listdir(sp_folder_path):

            sim_nr_folder_path = os.path.join(sp_folder_path, sim_nr_folder)
            if os.path.isdir(sim_nr_folder_path):
                output_file = os.path.join(sim_nr_folder_path, 'energy.csv')
                mcg_file = os.path.join(sim_nr_folder_path, 'mcg_results_skip_50')

                step, epot, ekin, etot, temp, vol = np.loadtxt(output_file, skiprows=1, unpack=True, delimiter=",")
                _, frame, time_ps, mcg_value = np.loadtxt(mcg_file, skiprows=1, unpack=True)

                sim_nr = int(sim_nr_folder.split("_")[-1])


                data[sp_nr][sim_nr] = {
                    'step': step,
                    'epot': epot,
                    'ekin': ekin,
                    'etot': etot,
                    'temp': temp,
                    'vol': vol,
                    'mcg_step': frame,
                    'mcg': mcg_value
                }

sp_list = list(sorted(data.keys()))
sim_nr_list = list(sorted(data[sp_list[0]].keys()))

plt.rcParams.update({
    # 'font.family': 'sans-serif',
    # 'font.sans-serif': ['Arial', 'Helvetica', 'sans-serif'],
    'font.family': 'serif',  # Use a LaTeX-compatible serif font
    'font.serif': ['Computer Modern Roman'],  # Match LaTeX's default font
    'figure.figsize': (6.5, 4.5),
    'figure.dpi': 100,
    'text.usetex': True,
    'font.size': 25,
    'axes.labelsize': 25,
    'axes.titlesize': 25,
    'xtick.labelsize': 25,
    'ytick.labelsize': 25,
    'legend.fontsize': 25,
    'mathtext.default':'regular',
})


# Determine the number of rows and columns for the grid
n_species = len(sp_list)
n_cols = 4  # Number of columns in the grid (adjust as needed)
n_rows = (n_species + n_cols - 1) // n_cols  # Calculate the number of rows

# Create a figure with subplots in a grid layout
fig, axs = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=(5 * n_cols, 5 * n_rows))

cmap = plt.cm.plasma  # You can change this to any other colormap
colors_plt = cmap(np.linspace(0, 1, len(sim_nr_list)))  # Generate a color cycle from the colormap

# Set the color cycle globally using rcParams
plt.rcParams.update({'axes.prop_cycle': plt.cycler('color', colors_plt)})


# Flatten the axs array for easy iteration (if n_rows > 1 or n_cols > 1)
if n_rows > 1 or n_cols > 1:
    axs = axs.flatten()
else:
    axs = [axs]  # Wrap in a list if there's only one subplot

# Iterate over each species and create a subplot
for i, sp in enumerate(sp_list):
    ax = axs[i]  # Get the current subplot axis
    for sim_nr in sim_nr_list:
        step = data[sp][sim_nr]['mcg_step']
        mcg = data[sp][sim_nr]['mcg']
        ax.plot(step*step2us, mcg)

    # Add labels, title, and legend to the subplot
    ax.set_xlabel(r"Time ($\mu$s)")
    ax.set_ylabel('MCG')
    ax.set_title(f'Start MCG {sp}')
    ax.axhline(20, label = 'Liquid State')
    ax.axhline(300, label = 'Solid State', color = colors['lightblue'])
    ax.legend(fontsize=15)
    ax.grid(True)

# Hide empty subplots (if any)
for j in range(i + 1, n_rows * n_cols):
    axs[j].axis('off')  # Turn off unused subplots

# Adjust layout to prevent overlap
plt.tight_layout()

# Show the plot

plt.savefig("CA_trajectories.png", bbox_inches='tight', dpi=300)  # Increase DPI as needed
#plt.show()


def state_definition(last_mcg):
    if last_mcg <= 20:
        return 'A'
    elif last_mcg >= 300:
        return 'B'


state_ratio_arr = np.zeros(len(sp_list))

for idx, sp in enumerate(sp_list):
    state_definitions = []
    for sim_nr in data[sp].keys():
        step = data[sp][sim_nr]['mcg_step']
        mcg = data[sp][sim_nr]['mcg']

        state = state_definition(mcg[-1])
        state_definitions.append(state)
    count_state_a = state_definitions.count("A")
    count_state_b = state_definitions.count("B")
    state_ratio = count_state_b / (count_state_a + count_state_b)
    state_ratio_arr[idx] = state_ratio


plt.scatter(sp_list, state_ratio_arr)
plt.xlabel("MCG")
plt.ylabel(r"$C_B / ( C_A + C_B)$")
plt.tight_layout()
plt.savefig("state_ratio.png")
