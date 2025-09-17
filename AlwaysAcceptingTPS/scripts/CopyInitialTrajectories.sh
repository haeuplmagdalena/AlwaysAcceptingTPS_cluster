#!/bin/bash

# Define the data folder FIRST
DATA_DIR="data_mu175_chain0_batch_3"

# Create the destination directory using the base name of DATA_DIR
DEST_DIR="InitialTrajectories/$(basename "$DATA_DIR")_last"
mkdir -p "$DEST_DIR"

# Loop through each run_* directory inside DATA_DIR
for run_dir in "$DATA_DIR"/run_*; do
    if [ -d "$run_dir" ]; then
        run_name=$(basename "$run_dir")
        mkdir -p "$DEST_DIR/$run_name"

        # Extract the maximum traj number
        max_num=$(ls "$run_dir"/traj_*.dcd 2>/dev/null | sed -E 's/.*traj_([0-9]+)\.dcd/\1/' | sort -n | tail -n 1)
        if [ -n "$max_num" ]; then
            target_file="$run_dir/traj_${max_num}.dcd"
            # Copy if the file exists
            if [ -f "$target_file" ]; then
                cp "$target_file" "$DEST_DIR/$run_name"/
            else
                echo "Warning: $target_file does not exist in $run_dir"
            fi
        fi
    fi
done
