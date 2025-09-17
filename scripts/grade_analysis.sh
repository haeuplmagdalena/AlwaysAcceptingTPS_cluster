#!/bin/bash

source /leonardo/home/userexternal/mhaeupl0/clathrate_env/bin/activate

# Define the base directory containing the first set of folders
BASE_DIR="/leonardo_scratch/fast/L-AUT_Coretti/CommittorAnalysisData"

# Define the Python script path
SCRIPT="/leonardo_work/L-AUT_Coretti/CODE/GRADE_added_rings/GRADE"

# Loop over each parent folder
for parent_folder in "$BASE_DIR"/*; do
    if [[ -d "$parent_folder" ]]; then  # Check if it's a directory
	folder_name=$(basename "$parent_folder")
        folder_number=$(echo "$folder_name" | grep -oE '[0-9]+' | head -1)

        # Check if the number is greater than 80
        if [[ "$folder_number" -gt 176 ]]; then
		echo "Entering parent folder: $parent_folder"

		# Loop over each subfolder
		for subfolder in "$parent_folder"/*; do
		    if [[ -d "$subfolder" ]]; then  # Check if it's a directory
			echo "Processing subfolder: $subfolder"

			# Find all .gro files in the subfolder
			for gro_file in "$subfolder"/*500.gro; do
			    if [[ -f "$gro_file" ]]; then  # Check if file exists
				echo "Running analysis on: $gro_file"

				# Run the script
				"$SCRIPT" \
				    -i "$gro_file" \
				    -f4 yes

				echo "Finished processing $gro_file"
			    fi
			done
		    fi
		    
		done
	fi
    fi
done

echo "All analyses completed!"

