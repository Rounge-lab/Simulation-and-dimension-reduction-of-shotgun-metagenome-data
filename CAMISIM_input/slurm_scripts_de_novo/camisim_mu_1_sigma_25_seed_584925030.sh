#!/bin/bash

#SBATCH --account=p1068_tsd
#SBATCH --job-name=CAMISIM
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20:00:00
#SBATCH --mem-per-cpu=10GB
#SBATCH --output=output_%j.txt
#SBATCH --error=error_%j.txt
#SBATCH --cpus-per-task=8

module --quiet purge  # Reset the modules to the system default

module load matplotlib/3.4.3-foss-2021b
module load Biopython/1.79-foss-2021b
module load biom-format/2.1.12-foss-2021b
module load scikit-learn/1.0.1-foss-2021b
module load SAMtools/1.14-GCC-11.2.0
module load Perl/5.34.0-GCCcore-11.2.0

cd /PATH/CAMISIM/

python metagenomesimulation.py ./config_files/config_mu_1_sigma_25_seed_58492503.ini