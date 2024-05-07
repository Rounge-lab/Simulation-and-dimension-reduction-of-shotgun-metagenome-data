# Imports
library(tidyverse)
library(vegan)

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

# Defining path to data
mag_path <- "/PATH/datasets/MAGs/raw/counts/median_coverage_genomes.tsv" # median number of reads per position in each genome, for each sample

# Read median coverage table
median_coverage_genomes <- read_delim(mag_path, delim = "\t") %>%
  rename(sample_id = `...1`) %>% 
  mutate(sample_id = gsub("-", "_", sample_id, fixed = TRUE)) %>% 
  filter(sample_id %in% data_by_sample.tbl$sample_id) %>%  # Only include baseline samples
  mutate(sample_id = paste("sample", seq_len(n()), sep = "_")) # New sample ids due to this being synthetic data

# -----------------------------------------------------------------------------
# Preparing profiles with different adjustments
# -----------------------------------------------------------------------------

# 10 MAGs 2x
df.wide_10_2x <- abundance_profile(coverage_df = median_coverage_genomes, 
                                   rand_mags_n = 10, 
                                   seed = 123, 
                                   adjustment_factor = 2.0) %>% 
  pivot_wider(names_from = MAG, 
              values_from = rel_abundance, 
              id_cols = c(sample_id, group))

# 20 MAGs 5x
df.wide_20_5x <- abundance_profile(coverage_df = median_coverage_genomes, 
                                   rand_mags_n = 20, 
                                   seed = 123, 
                                   adjustment_factor = 5.0) %>% 
  pivot_wider(names_from = MAG, 
              values_from = rel_abundance, 
              id_cols = c(sample_id, group))

# 50 MAGs 10x
df.wide_50_10x <- abundance_profile(coverage_df = median_coverage_genomes, 
                                    rand_mags_n = 50, 
                                    seed = 123, 
                                    adjustment_factor = 10.0) %>% 
  pivot_wider(names_from = MAG, 
              values_from = rel_abundance, 
              id_cols = c(sample_id, group))

# 100 MAGs 20x
df.wide_100_20x <- abundance_profile(coverage_df = median_coverage_genomes, 
                                     rand_mags_n = 100, 
                                     seed = 123, 
                                     adjustment_factor = 20.0) %>% 
  pivot_wider(names_from = MAG, 
              values_from = rel_abundance, 
              id_cols = c(sample_id, group))

# -----------------------------------------------------------------------------
# Preparing data with real variables
# -----------------------------------------------------------------------------

# Variable = "kjonn" (Male or Female)
data_by_sample.tbl <- readRDS("/PATH/datasets/data_by_FIT_sample/data_by_sample.Rds") %>% 
  select(sample_id, deltaker_id, Prøvetype, Total_Bases_QC_ATLAS) %>% 
  filter(Prøvetype == "Baseline") %>% 
  filter(Total_Bases_QC_ATLAS > 1e9)

screening_proc <- readRDS("/PATH/CRCbiome/datasets/screening/241122_screening_proc.rds") %>% 
  select(deltaker_id, kjonn)

df <- left_join(data_by_sample.tbl, screening_proc, by = "deltaker_id") %>% 
  select(sample_id, kjonn)

df.wide_crcbiome_kjonn <- read_delim(mag_path, delim = "\t") %>%
  rename(sample_id = `...1`) %>% 
  mutate(sample_id = gsub("-", "_", sample_id, fixed = TRUE)) %>% 
  filter(sample_id %in% data_by_sample.tbl$sample_id) %>%  # Only include baseline samples
  left_join(df, by = "sample_id") %>% 
  select(sample_id, kjonn, everything()) %>% 
  # Convert data to long format
  pivot_longer(
    cols = -c(sample_id, kjonn), # exclude the id and group from the gathering
    names_to = "MAG",
    values_to = "median_coverage") %>% 
  group_by(sample_id) %>%
  mutate(sum_measurements = sum(median_coverage, na.rm = TRUE)) %>% # Calculate the sum of measurements for each sample
  mutate(rel_abundance = (median_coverage/ sum_measurements)) %>% # Calculate the relative abundance by total reads and sum of measurements
  ungroup() %>% 
  select(-c(sum_measurements)) %>% 
  pivot_wider(names_from = MAG, 
              values_from = rel_abundance, 
              id_cols = c(sample_id, kjonn))

# Variable = "senter" (Moss or Bærum)
data_by_sample.tbl <- readRDS("/PATH/datasets/data_by_FIT_sample/data_by_sample.Rds") %>% 
  select(sample_id, deltaker_id, Prøvetype, Total_Bases_QC_ATLAS) %>% 
  filter(Prøvetype == "Baseline") %>% 
  filter(Total_Bases_QC_ATLAS > 1e9)

screening_proc <- readRDS("/PATH/CRCbiome/datasets/screening/241122_screening_proc.rds") %>% 
  select(deltaker_id, senter)

df <- left_join(data_by_sample.tbl, screening_proc, by = "deltaker_id") %>% 
  select(sample_id, senter)

df.wide_crcbiome_senter <- read_delim(mag_path, delim = "\t") %>%
  rename(sample_id = `...1`) %>% 
  mutate(sample_id = gsub("-", "_", sample_id, fixed = TRUE)) %>% 
  filter(sample_id %in% data_by_sample.tbl$sample_id) %>%  # Only include baseline samples
  left_join(df, by = "sample_id") %>% 
  select(sample_id, senter, everything()) %>% 
  # Convert data to long format
  pivot_longer(
    cols = -c(sample_id, senter), # exclude the id and group from the gathering
    names_to = "MAG",
    values_to = "median_coverage") %>% 
  group_by(sample_id) %>%
  mutate(sum_measurements = sum(median_coverage, na.rm = TRUE)) %>% # Calculate the sum of measurements for each sample
  mutate(rel_abundance = (median_coverage/ sum_measurements)) %>% # Calculate the relative abundance by total reads and sum of measurements
  ungroup() %>% 
  select(-c(sum_measurements)) %>% 
  pivot_wider(names_from = MAG, 
              values_from = rel_abundance, 
              id_cols = c(sample_id, senter))

# -----------------------------------------------------------------------------
# Preparing MetaPhlAn data from simulations
# -----------------------------------------------------------------------------

# Paths to MetaPhlAn data
path_20_5x <- "/PATH/metaphlan_output/20_5x/tables/profiles.tsv"
path_50_10x <- "/PATH/metaphlan_output/50_10x/tables/profiles.tsv"

# Transforming values to relative abundance values and transforming to correct format
metaphlan_20_5x.wide <- process_metaphlan_table_rel_abundance_filtering(path_20_5x) %>% 
  process_metaphlan_table_pcoa_profile_based_sim() %>% 
  mutate(group = as.factor(if_else(as.numeric(str_extract(sample_id, "\\d+")) >= 600, 2, 1))) %>% # All sample ids above 600 are group 2
  select(sample_id, group, everything())

metaphlan_50_10x.wide <- process_metaphlan_table_rel_abundance_filtering(path_50_10x) %>% 
  process_metaphlan_table_pcoa_profile_based_sim() %>% 
  mutate(group = as.factor(if_else(as.numeric(str_extract(sample_id, "\\d+")) >= 600, 2, 1))) %>% # All sample ids above 600 are group 2
  select(sample_id, group, everything())

# -----------------------------------------------------------------------------
# Permanova analysis
# -----------------------------------------------------------------------------
permanova_10_2x <- adonis2(formula = df.wide_10_2x[, -c(1:2)] ~ group,
                             data = df.wide_10_2x,
                             method = "bray",
                             permutations = 1000)

permanova_20_5x <- adonis2(formula = df.wide_20_5x[, -c(1:2)] ~ group,
                             data = df.wide_20_5x,
                             method = "bray",
                             permutations = 1000)

permanova_50_10x <- adonis2(formula = df.wide_50_10x[, -c(1:2)] ~ group,
                             data = df.wide_50_10x,
                             method = "bray",
                             permutations = 1000)

permanova_100_20x <- adonis2(formula = df.wide_100_20x[, -c(1:2)] ~ group,
                             data = df.wide_100_20x,
                             method = "bray",
                             permutations = 1000)

permanova_CRCbiome_kjonn <- adonis2(formula = df.wide_crcbiome_kjonn[, -c(1:2)] ~ kjonn,
                                    data = df.wide_crcbiome_kjonn,
                                    method = "bray",
                                    permutations = 1000)

permanova_CRCbiome_senter <- adonis2(formula = df.wide_crcbiome_senter[, -c(1:2)] ~ senter,
                                     data = df.wide_crcbiome_senter,
                                     method = "bray",
                                     permutations = 1000)

permanova_20_5x_metaphlan <- adonis2(formula = metaphlan_20_5x.wide[, -c(1:2)] ~ group,
                                     data = metaphlan_20_5x.wide,
                                     method = "bray",
                                     permutations = 1000)

permanova_50_10x_metaphlan <- adonis2(formula = metaphlan_50_10x.wide[, -c(1:2)] ~ group,
                                      data = metaphlan_50_10x.wide,
                                      method = "bray",
                                      permutations = 1000)
