#!/bin/bash

#SBATCH --account=p1068_tsd
#SBATCH --job-name=rename_files_dryrun
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=10GB

# Base directory where the sample folders are located
base_directory="/PATH/CAMISIM_output/some_directory/"

cd "$base_directory" # Ensure you are using the variable with a $ sign

for file in *sample_*; do
  # Extract the part before '_sample_' and the number after it
  prefix=$(echo $file | sed -e 's/_sample_[0-9]*//')
  number=$(echo $file | grep -o -E '_sample_([0-9]+)' | sed -e 's/_sample_//')

  # Add 500 to the number
  new_number=$((10#$number + 1100)) # The 10# is to ensure the number is treated as decimal

  # Construct the new filename
  new_file="${prefix}_sample_${new_number}"

  # Print old and new filenames instead of renaming
  echo "Would rename '$file' to '$new_file'"
done
