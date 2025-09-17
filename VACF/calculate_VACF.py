import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import h5py
import os

import time
from tqdm import tqdm

# %%
structure = 'Amorphous'

# %%
mcg = 300

# %%
topology = '../ClathrateTPS/gro_files/conf.gro'

# %%

# %%
def autocorrelation(data, maxlag=None, stride=1):
    """
    Compute the autocorrelation function (VACF) with a stride over the lag.
    
    Parameters:
        data (numpy.ndarray): Velocity array of shape (n_frames, n_atoms, 3).
        maxlag (int, optional): Maximum lag to compute. Defaults to total frames.
        stride (int): Step size between computed lags (must be >=1). Defaults to 1.
    
    Returns:
        numpy.ndarray: Autocorrelation function over time at specified lags.
    """
    n_frames = data.shape[0]
    
    # Set maxlag to n_frames if not specified or too large
    if maxlag is None or maxlag > n_frames:
        print(f'Setting maxlag to {n_frames}!')
        maxlag = n_frames
    
    # Validate stride
    stride = max(1, stride)  # Ensure stride is at least 1
    lags = range(0, maxlag, stride)  # Lags to compute (0, stride, 2*stride, ...)
    n_lags = len(lags)
    vacf = np.zeros(n_lags)
    
    print(f'Computing {n_lags} lags (stride={stride}) up to maxlag {maxlag}')
    
    # Compute autocorrelation only at selected lags
    for idx, t in enumerate(lags):
        print(f"\r{idx + 1}/{n_lags} (lag={t})", end="", flush=True)
        
        # Slice data for current lag
        v_t = data[:n_frames - t]  # Initial frames: [0, n_frames - t - 1]
        v_t_dt = data[t:]          # Shifted frames: [t, n_frames - 1]
        
        # Dot product and averaging
        dot_product = np.einsum('ijk,ijk->ij', v_t, v_t_dt)  # Sum over coordinates
        particle_mean = np.mean(dot_product, axis=1)          # Average over atoms
        vacf[idx] = np.mean(particle_mean)                   # Average over time
    
    print()  # Newline after progress
    vacf /= vacf[0]  # Normalize by t=0 value
    return vacf


# %% [markdown]
# ## Short timescale

# %%
stride = 10

burst = 1

#end = 100000

# %%
file_path = f'/leonardo_scratch/large/userexternal/mhaeupl0/InitialTrajectories/{structure}/cv_{mcg}/N_out_1_N_10e5/trajectory.h5'

# %%
with h5py.File(file_path, 'r') as f:
    print(f["coordinates"].shape) 

# %%
with h5py.File(file_path, "r") as f:
    total_frames = f["coordinates"].shape[0]
    n_atoms = f["coordinates"].shape[1]
    
    # Calculate indices you're going to load
    indices = np.arange(0, total_frames, stride)
    n_frames = len(indices)

    # Preallocate arrays
    coordinates = np.empty((n_frames, n_atoms, 3), dtype=np.float32)
    velocities = np.empty((n_frames, n_atoms, 3), dtype=np.float32)

    for i, idx in enumerate(tqdm(indices, desc="Loading trajectory")):
        coordinates[i] = f["coordinates"][idx]
        velocities[i] = f["velocities"][idx]

# %%
lag_stride = 10

# %%
acf_vel = autocorrelation(velocities, maxlag = None, stride = lag_stride)
acf_pos = autocorrelation(coordinates, maxlag = None, stride = lag_stride)


np.save('acf_vel', acf_vel)
np.save('acf_pos', acf_pos)
