# Imports
library(tidyverse)
library(vegan)
library(glue)
library(patchwork)

# Source function file
source("/PATH/scripts/functions/all_functions.R")

# -----------------------------------------------------------------------------
# Preparing CRCbiome data
# -----------------------------------------------------------------------------

# Defining path to CRCbiome MetaPhlAn data
metaphlan_path_crcbiome <- "/PATH/datasets/taxonomy/profiles.tsv"

# Processing and filtering CRCbiome MetaPhlAn data
metaphlan.crcbiome <- process_metaphlan_table_rel_abundance_filtering(metaphlan_path_crcbiome)

# Processing CRCbiome MetaPhlAn data for PCoA
metaphlan_crcbiome.tbl <- process_metaphlan_table_pcoa(metaphlan.crcbiome, "CRCbiome")

# Extracting a subset of 50 samples from the PCoA processed CRCbiome MetaPhlAn data
metaphlan_crcbiome_subset.tbl <- get_random_samples(metaphlan_crcbiome.tbl, 50)

# Extracting 10 random samples from CRCbiome MetaPhlAn data to use for "species per sample" and "samples per species"
set.seed(123)
num_columns_to_sample <- 10
sampled_cols <- sample(2:ncol(metaphlan.crcbiome), num_columns_to_sample)
final_cols <- c(1, sampled_cols)
metaphlan.crcbiome_subset_10 <- metaphlan.crcbiome[, final_cols]

# -----------------------------------------------------------------------------
# Collecting the data from mu simulations
# -----------------------------------------------------------------------------
# Defining paths to MetaPhlAn data
metaphlan_path_mu_1 <- "/PATH/metaphlan_output/mu_1_sigma_2_seed_632741178/tables/profiles.tsv"
metaphlan_path_mu_2 <- "/PATH/metaphlan_output/mu_2_sigma_2_seed_632741178/tables/profiles.tsv"
metaphlan_path_mu_3 <- "/PATH/metaphlan_output/mu_3_sigma_2_seed_632741178/tables/profiles.tsv"

# Processing and filtering MetaPhlAn data
metaphlan.mu_1 <- process_metaphlan_table_rel_abundance_filtering(metaphlan_path_mu_1)
metaphlan.mu_2 <- process_metaphlan_table_rel_abundance_filtering(metaphlan_path_mu_2)
metaphlan.mu_3 <- process_metaphlan_table_rel_abundance_filtering(metaphlan_path_mu_3)

# -----------------------------------------------------------------------------
# PCoA for the mu simulation data
# -----------------------------------------------------------------------------

# Performing PCoA on MetaPhlAn data
pcoa_mu <- perform_pcoa_and_plot_3_df(metaphlan_df_1 = metaphlan.mu_1,
                                      metaphlan_df_2 = metaphlan.mu_2,
                                      metaphlan_df_3 = metaphlan.mu_3,
                                      metaphlan_df_crc = metaphlan_crcbiome_subset.tbl,
                                      source_1 = "mu = 1",
                                      source_2 = "mu = 2",
                                      source_3 = "mu = 3",
                                      method = "bray",
                                      k = 2,
                                      remove_columns = c(1:2))

pcoa_plot_mu <- pcoa_mu$pcoa_plot + labs(color = "")

# -----------------------------------------------------------------------------
# Samples per species for the mu simulation data
# -----------------------------------------------------------------------------

# Calculating samples per species for each dataset, and creating histograms
samples_per_species_mu <- samples_per_species_histograms(metaphlan.mu_1, 
                                                         metaphlan.mu_2, 
                                                         metaphlan.mu_3, 
                                                         metaphlan.crcbiome_subset_10,
                                                         source_1 = "mu = 1", 
                                                         source_2 = "mu = 2", 
                                                         source_3 = "mu = 3")

histogram_grid_mu <- samples_per_species_mu$histogram_grid

wilcoxon_results_mu <- samples_per_species_mu$wilcoxon_res

# -----------------------------------------------------------------------------
# Species per sample for the mu simulation data
# -----------------------------------------------------------------------------

# Calculating species per sample for each dataset, and creating box plot
species_per_sample_mu <- species_per_sample_boxplot(metaphlan_df_1 = metaphlan.mu_1,
                                                    metaphlan_df_2 = metaphlan.mu_2,
                                                    metaphlan_df_3 = metaphlan.mu_3,
                                                    metaphlan_df_crc = metaphlan.crcbiome_subset_10,
                                                    source_1 = "mu = 1", 
                                                    source_2 = "mu = 2", 
                                                    source_3 = "mu = 3")
box_plot_mu <- species_per_sample_mu$box_plot

wilcoxon_results_mu <- species_per_sample_mu$wilcoxon_res

# -----------------------------------------------------------------------------
# Gather all plots from mu simulation in one grid
# -----------------------------------------------------------------------------

plot <- ((pcoa_plot_mu + box_plot_mu) / histogram_grid_mu) + 
  plot_annotation(tag_levels = 'A')

# -----------------------------------------------------------------------------
# Collecting the data from seed simulations
# -----------------------------------------------------------------------------

# Defining paths to MetaPhlAn data
metaphlan_path_seed_632741178 <- "/PATH/metaphlan_output/mu_1_sigma_25_seed_632741178/tables/profiles.tsv"
metaphlan_path_seed_584925030 <- "/PATH/metaphlan_output/mu_1_sigma_25_seed_584925030/tables/profiles.tsv"
metaphlan_path_seed_470267527 <- "/PATH/metaphlan_output/mu_1_sigma_25_seed_470267527/tables/profiles.tsv"

# Processing and filtering MetaPhlAn data
metaphlan.seed_632741178 <- process_metaphlan_table_rel_abundance_filtering(metaphlan_path_seed_632741178)
metaphlan.seed_584925030 <- process_metaphlan_table_rel_abundance_filtering(metaphlan_path_seed_584925030)
metaphlan.seed_470267527 <- process_metaphlan_table_rel_abundance_filtering(metaphlan_path_seed_470267527)

# -----------------------------------------------------------------------------
# PCoA for the seed simulation data
# -----------------------------------------------------------------------------

# Performing PCoA on MetaPhlAn data
pcoa_seed <- perform_pcoa_and_plot_3_df(metaphlan_df_1 = metaphlan.seed_632741178,
                                        metaphlan_df_2 = metaphlan.seed_584925030,
                                        metaphlan_df_3 = metaphlan.seed_470267527,
                                        metaphlan_df_crc = metaphlan_crcbiome_subset.tbl,
                                        source_1 = "seed = 632741178",
                                        source_2 = "seed = 584925030",
                                        source_3 = "seed = 470267527",
                                        method = "bray",
                                        k = 2,
                                        remove_columns = c(1:2))

pcoa_plot_seed <- pcoa_seed$pcoa_plot + labs(color = "")

# -----------------------------------------------------------------------------
# Samples per species for the seed simulation data
# -----------------------------------------------------------------------------

# Calculating samples per species for each dataset, and creating histograms
samples_per_species_seed <- samples_per_species_histograms(metaphlan_df_1 = metaphlan.seed_632741178,
                                                           metaphlan_df_2 = metaphlan.seed_584925030,
                                                           metaphlan_df_3 = metaphlan.seed_470267527,
                                                           metaphlan_df_crc = metaphlan.crcbiome_subset_10,
                                                           source_1 = "seed = 632741178",
                                                           source_2 = "seed = 584925030",
                                                           source_3 = "seed = 470267527")

histogram_grid_seed <- samples_per_species_seed$histogram_grid

wilcoxon_results_seed <- samples_per_species_seed$wilcoxon_res

# -----------------------------------------------------------------------------
# Species per sample for the seed simulation data
# -----------------------------------------------------------------------------

# Calculating species per sample for each dataset, and creating box plot
species_per_sample_seed <- species_per_sample_boxplot(metaphlan_df_1 = metaphlan.seed_632741178,
                                                      metaphlan_df_2 = metaphlan.seed_584925030,
                                                      metaphlan_df_3 = metaphlan.seed_470267527,
                                                      metaphlan_df_crc = metaphlan.crcbiome_subset_10,
                                                      source_1 = "seed = 632741178",
                                                      source_2 = "seed = 584925030",
                                                      source_3 = "seed = 470267527")

box_plot_seed <- species_per_sample_seed$box_plot + theme(axis.text.x = element_text(angle = 20, hjust = 1.1, vjust=1))

wilcoxon_results_seed <- species_per_sample_seed$wilcoxon_res

# -----------------------------------------------------------------------------
# Gather all plots from seed simulation in one grid
# -----------------------------------------------------------------------------

plot <- ((pcoa_plot_seed + box_plot_seed) / histogram_grid_seed) + 
  plot_annotation(tag_levels = 'A')

# -----------------------------------------------------------------------------
# Collecting the data from sigma simulations
# -----------------------------------------------------------------------------
# Defining paths to MetaPhlAn data
metaphlan_path_sigma_2 <- "/PATH/metaphlan_output/mu_1_sigma_2_seed_632741178/tables/profiles.tsv"
metaphlan_path_sigma_15 <- "/PATH/metaphlan_output/mu_1_sigma_15_seed_632741178/tables/profiles.tsv"
metaphlan_path_sigma_25 <- "/PATH/metaphlan_output/mu_1_sigma_25_seed_632741178/tables/profiles.tsv"

# Processing and filtering MetaPhlAn data
metaphlan.sigma_2 <- process_metaphlan_table_rel_abundance_filtering(metaphlan_path_sigma_2)
metaphlan.sigma_15 <- process_metaphlan_table_rel_abundance_filtering(metaphlan_path_sigma_15)
metaphlan.sigma_25 <- process_metaphlan_table_rel_abundance_filtering(metaphlan_path_sigma_25)

# -----------------------------------------------------------------------------
# PCoA for the sigma simulation data
# -----------------------------------------------------------------------------

# Performing PCoA on MetaPhlAn data
pcoa_sigma <- perform_pcoa_and_plot_3_df(metaphlan_df_1 = metaphlan.sigma_2,
                                         metaphlan_df_2 = metaphlan.sigma_15,
                                         metaphlan_df_3 = metaphlan.sigma_25,
                                         metaphlan_df_crc = metaphlan_crcbiome_subset.tbl,
                                         source_1 = "sigma = 2",
                                         source_2 = "sigma = 1.5",
                                         source_3 = "sigma = 2.5",
                                         method = "bray",
                                         k = 2,
                                         remove_columns = c(1:2))

pcoa_plot_sigma <- pcoa_sigma$pcoa_plot + labs(color = "")

# -----------------------------------------------------------------------------
# Samples per species for the sigma simulation data
# -----------------------------------------------------------------------------

# Calculating samples per species for each dataset, and creating histograms
samples_per_species_sigma <- samples_per_species_histograms(metaphlan_df_1 = metaphlan.sigma_2, 
                                                            metaphlan_df_2 = metaphlan.sigma_15, 
                                                            metaphlan_df_3 = metaphlan.sigma_25, 
                                                            metaphlan_df_crc = metaphlan.crcbiome_subset_10,
                                                            source_1 = "sigma = 2", 
                                                            source_2 = "sigma = 1.5", 
                                                            source_3 = "sigma = 2.5")

histogram_grid_sigma <- samples_per_species_sigma$histogram_grid

wilcoxon_results_sigma <- samples_per_species_sigma$wilcoxon_res

# -----------------------------------------------------------------------------
# Species per sample for the sigma simulation data
# -----------------------------------------------------------------------------

# Calculating species per sample for each dataset, and creating box plot
species_per_sample_sigma <- species_per_sample_boxplot(metaphlan_df_1 = metaphlan.sigma_2,
                                                       metaphlan_df_2 = metaphlan.sigma_15,
                                                       metaphlan_df_3 = metaphlan.sigma_25,
                                                       metaphlan_df_crc = metaphlan.crcbiome_subset_10,
                                                       source_1 = "sigma = 2",
                                                       source_2 = "sigma = 1.5",
                                                       source_3 = "sigma = 2.5")

box_plot_sigma <- species_per_sample_sigma$box_plot + theme(axis.text.x = element_text(angle = 40, hjust = 1.1, vjust=1))

wilcoxon_results_sigma <- species_per_sample_sigma$wilcoxon_res

# -----------------------------------------------------------------------------
# Gather all plots from sigma simulation in one grid
# -----------------------------------------------------------------------------

plot <- ((pcoa_plot_sigma + box_plot_sigma) / histogram_grid_sigma) + 
  plot_annotation(tag_levels = 'A')
