import h5py
import numpy as np
from tqdm import tqdm

structure = 'Amorphous'
mcg = 300

file_path = f'/leonardo_scratch/large/userexternal/mhaeupl0/InitialTrajectories/{structure}/cv_{mcg}/N_out_1_N_10e5/trajectory.h5'


input_file = file_path
output_file = "output_short.hdf5"
start_frame = 0
end_frame = 200000  # inclusive
chunk_size = 1000   # Adjust based on memory (1000-10000)
position_path = "coordinates"  # Path to position dataset
velocity_path = "velocities"   # Path to velocity dataset

with h5py.File(input_file, "r") as f_in, h5py.File(output_file, "w") as f_out:
    # Copy root attributes
    for key, value in f_in.attrs.items():
        f_out.attrs[key] = value

    # Copy non-dataset objects (groups, metadata)
    def copy_items(name, obj):
        if isinstance(obj, h5py.Group):
            # Create group and copy attributes
            group = f_out.create_group(name)
            for key, value in obj.attrs.items():
                group.attrs[key] = value
        elif isinstance(obj, h5py.Dataset) and name not in [position_path, velocity_path]:
            # Copy non-trajectory datasets
            f_in.copy(name, f_out, name=name)

    f_in.visititems(copy_items)

    # Process trajectory datasets
    for path in [position_path, velocity_path]:
        dset_in = f_in[path]
        n_frames = end_frame - start_frame + 1
        total_chunks = int(np.ceil(n_frames / chunk_size))
        
        # Create output dataset with same properties
        dset_out = f_out.create_dataset(
            path,
            shape=(n_frames,) + dset_in.shape[1:],
            dtype=dset_in.dtype,
            chunks=dset_in.chunks,
            compression=dset_in.compression,
            compression_opts=dset_in.compression_opts
        )
        
        # Initialize progress bar
        pbar = tqdm(
            total=n_frames,
            desc=f"Copying {path.split('/')[-1]}",
            unit="frame",
            dynamic_ncols=True
        )
        
        # Copy data in chunks
        for i in range(0, n_frames, chunk_size):
            chunk_end = min(i + chunk_size, n_frames)
            chunk = dset_in[start_frame + i : start_frame + chunk_end]
            dset_out[i:chunk_end] = chunk
            pbar.update(len(chunk))  # Update progress bar
            
        pbar.close()
        
        # Copy dataset attributes
        for key, value in dset_in.attrs.items():
            dset_out.attrs[key] = value

print(f"\nSaved frames {start_frame}-{end_frame} to {output_file}")
