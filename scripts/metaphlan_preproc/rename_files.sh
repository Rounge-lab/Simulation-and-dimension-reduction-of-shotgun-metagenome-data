#!/bin/bash

#SBATCH --account=p1068_tsd
#SBATCH --job-name=rename_files
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=10GB

# Base directory where the sample folders with reads are located
base_directory="/PATH/CAMISIM_output/some_directory/"

cd "$base_directory"

for file in *sample_*; do
  # Extract the part before '_sample_' and the number after it
  prefix=$(echo $file | sed -e 's/_sample_[0-9]*//')
  number=$(echo $file | grep -o -E '_sample_([0-9]+)' | sed -e 's/_sample_//')

  # Add given number to the sample id
  add_number=100
  new_number=$((10#$number + $add_number)) # The 10# is to ensure the number is treated as decimal

  # Construct the new filename
  new_file="${prefix}_sample_${new_number}"

  # Rename the file
  mv "$file" "$new_file"
done
