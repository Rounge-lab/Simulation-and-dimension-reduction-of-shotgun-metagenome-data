# This code can be used to write the abundance profiles to file, and find the full paths of the abundance.
# These paths are needed in the configuration file of the CAMISIM pipeline.

# Imports
library(tidyverse)

# Source function file
source("/PATH/scripts/functions/all_functions.R")

# -----------------------------------------------------------------------------
# Preparing CRCbiome data to make profiles
# -----------------------------------------------------------------------------

# Reading data by sample data. Collecting only baseline data
data_by_sample.tbl <- readRDS("/PATH/datasets/data_by_FIT_sample/data_by_sample.Rds") %>% 
  select(sample_id, Prøvetype, Total_Bases_QC_ATLAS) %>% 
  filter(Prøvetype == "Baseline") %>% 
  filter(Total_Bases_QC_ATLAS > 1e9)

# Read median coverage table
mag_path <- "/PATH/datasets/MAGs/raw/counts/median_coverage_genomes.tsv"

# Reading tsv file
median_coverage_genomes <- read_delim(mag_path, delim = "\t") %>%
  rename(sample_id = `...1`) %>% 
  mutate(sample_id = gsub("-", "_", sample_id, fixed = TRUE)) %>% 
  filter(sample_id %in% data_by_sample.tbl$sample_id) %>%  # Only include baseline samples
  mutate(sample_id = paste("sample", seq_len(n()), sep = "_")) # New sample ids due to this being synthetic data

# -----------------------------------------------------------------------------
# Processing data
# -----------------------------------------------------------------------------

# Making the 20 MAGs 5x abundance profiles
df.long_20_5x <- abundance_profile(coverage_df = median_coverage_genomes, 
                                   rand_mags_n = 20, 
                                   seed = 123, 
                                   adjustment_factor = 5.0)

# Making the 50 MAGs 10x abundance profiles
df.long_50_10x <- abundance_profile(coverage_df = median_coverage_genomes, 
                                    rand_mags_n = 50, 
                                    seed = 123, 
                                    adjustment_factor = 10.0)


# Defining output path for abundance profiles (needs to exist)
out_dir <- "/PATH/CAMISIM_input/some_directory/"

# Writing to file
writing <- write_profiles_to_file(df.long_50_10x, out_dir) # Remember to change the abundance profile data

# Extracting the paths of the first 100 samples
writing$group2_tsv_files[1:100]

# Use this code to create a string of paths separated by comma (as required by CAMISIM)
paste(writing$group2_tsv_files[1:10], collapse = ",")
