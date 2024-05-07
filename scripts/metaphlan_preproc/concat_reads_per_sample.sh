#!/bin/bash

#SBATCH --account=p1068_tsd
#SBATCH --job-name=concat_fastq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0-04:00:00
#SBATCH --mem-per-cpu=10GB

# Base directory where the sample folders with reads from CAMISIM are located
base_directory="/PATH/CAMISIM_output/some_directory/"

# Output directory for the concatenated files
output_directory="/PATH/metaphlan_input/fastq/some_directory"

# Create the output directory if it doesn't exist
mkdir -p "$output_directory"

# Loop through each sample directory matching the pattern
for sample_dir in "$base_directory"/*_sample_*; do
    if [ -d "$sample_dir/reads" ]; then
        echo "Processing: $sample_dir"

        # Extract the sample name from the folder name
        folder_name=$(basename "$sample_dir")
        sample_name=$(echo $folder_name | grep -o 'sample_[0-9]*')

        # Change to reads directory
        cd "$sample_dir/reads"

        # Concatenate and compress .01.fq.gz files into sample_x_R1.fastq.gz in the output directory
        zcat *01.fq.gz | gzip > "$output_directory/${sample_name}_R1.fastq.gz"
        echo "Created ${sample_name}_R1.fastq.gz in the output directory"

        # Concatenate and compress .02.fq.gz files into sample_x_R2.fastq.gz in the output directory
        zcat *02.fq.gz | gzip > "$output_directory/${sample_name}_R2.fastq.gz"
        echo "Created ${sample_name}_R2.fastq.gz in the output directory"

    else
        echo "Reads directory not found in $sample_dir"
    fi
done

echo "Processing complete."
