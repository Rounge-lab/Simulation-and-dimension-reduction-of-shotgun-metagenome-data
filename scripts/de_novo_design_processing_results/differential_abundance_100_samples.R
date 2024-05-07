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

# Extracting a subset of 100 samples from the PCoA processed CRCbiome MetaPhlAn data
metaphlan_crcbiome_subset.tbl <- get_random_samples(metaphlan_crcbiome.tbl, 100)

# Extracting 100 random samples from CRCbiome MetaPhlAn data to use for "species per sample" and "samples per species"
# (These are the same samples used in the PCoA, but in a different format)
set.seed(123)
num_columns_to_sample <- 100
sampled_cols <- sample(2:ncol(metaphlan.crcbiome), num_columns_to_sample)
final_cols <- c(1, sampled_cols)
metaphlan.crcbiome_subset_100 <- metaphlan.crcbiome[, final_cols]

# -----------------------------------------------------------------------------
# Collecting the data from simulations
# -----------------------------------------------------------------------------
# Defining paths to MetaPhlAn data
metaphlan_path_seed_632741178 <- "/PATH/metaphlan_output/mu_1_sigma_25_seed_632741178_100_samples/tables/profiles.tsv"
metaphlan_path_seed_584925030 <- "/PATH/metaphlan_output/mu_1_sigma_25_seed_584925030_100_samples/tables/profiles.tsv"
metaphlan_path_seed_470267527 <- "/PATH/metaphlan_output/mu_1_sigma_25_seed_470267527_100_samples/tables/profiles.tsv"

# Processing and filtering MetaPhlAn data
metaphlan.seed_632741178 <- process_metaphlan_table_rel_abundance_filtering(metaphlan_path_seed_632741178)
metaphlan.seed_584925030 <- process_metaphlan_table_rel_abundance_filtering(metaphlan_path_seed_584925030)
metaphlan.seed_470267527 <- process_metaphlan_table_rel_abundance_filtering(metaphlan_path_seed_470267527)

# -----------------------------------------------------------------------------
# PCoA 
# -----------------------------------------------------------------------------

# Performing PCoA on all the MetaPhlAn data
pcoa <- perform_pcoa_and_plot_3_df(metaphlan_df_1 = metaphlan.seed_632741178,
                                   metaphlan_df_2 = metaphlan.seed_584925030,
                                   metaphlan_df_3 = metaphlan.seed_470267527,
                                   metaphlan_df_crc = metaphlan_crcbiome_subset.tbl,
                                   source_1 = "seed = 632741178", 
                                   source_2 = "seed = 584925030", 
                                   source_3 = "seed = 470267527",
                                   method = "bray",
                                   k = 2,
                                   remove_columns = c(1:2))

pcoa_plot <- pcoa$pcoa_plot

# -----------------------------------------------------------------------------
# Samples per species
# -----------------------------------------------------------------------------

# Calculating samples per species for each dataset, and creating histograms
samples_per_species <- samples_per_species_histograms(metaphlan_df_1 = metaphlan.seed_632741178,
                                                      metaphlan_df_2 = metaphlan.seed_584925030,
                                                      metaphlan_df_3 = metaphlan.seed_470267527,
                                                      metaphlan_df_crc = metaphlan.crcbiome_subset_100,
                                                      source_1 = "seed = 632741178", 
                                                      source_2 = "seed = 584925030", 
                                                      source_3 = "seed = 470267527")

histogram_grid <- samples_per_species$histogram_grid

wilcoxon_results <- samples_per_species$wilcoxon_res

# -----------------------------------------------------------------------------
# Species per sample
# -----------------------------------------------------------------------------

# Calculating species per sample for each dataset, and creating box plot
species_per_sample <- species_per_sample_boxplot(metaphlan_df_1 = metaphlan.seed_632741178,
                                                 metaphlan_df_2 = metaphlan.seed_584925030,
                                                 metaphlan_df_3 = metaphlan.seed_470267527,
                                                 metaphlan_df_crc = metaphlan.crcbiome_subset_100,
                                                 source_1 = "seed = 632741178", 
                                                 source_2 = "seed = 584925030", 
                                                 source_3 = "seed = 470267527")

box_plot <- species_per_sample$box_plot + theme(axis.text.x = element_text(angle = 20, hjust = 1.1, vjust=1))

wilcoxon_results <- species_per_sample$wilcoxon_res

# -----------------------------------------------------------------------------
# Gathering all plots in one grid
# -----------------------------------------------------------------------------

plot <- ((pcoa_plot + box_plot) / histogram_grid)+ 
  plot_annotation(tag_levels = 'A')

# -----------------------------------------------------------------------------
# Species abundance and prevalence in each dataset
# -----------------------------------------------------------------------------

# Calculating abundances for all species in each dataset
all_abundances_a <- calculate_species_abundances(metaphlan.seed_632741178)$all_species_abundances
sorted_species_a <- calculate_species_abundances(metaphlan.seed_632741178)$sorted_species

all_abundances_b <- calculate_species_abundances(metaphlan.seed_584925030)$all_species_abundances
sorted_species_b <- calculate_species_abundances(metaphlan.seed_584925030)$sorted_species

all_abundances_c <- calculate_species_abundances(metaphlan.seed_470267527)$all_species_abundances
sorted_species_c <- calculate_species_abundances(metaphlan.seed_470267527)$sorted_species

# Plotting bar plot of abundances
abundance_plot <- plot_abundance(all_abundances_a = all_abundances_a,
                                 all_abundances_b = all_abundances_b,
                                 all_abundances_c = all_abundances_c,
                                 sorted_species_a = sorted_species_a,
                                 sorted_species_b = sorted_species_b,
                                 sorted_species_c = sorted_species_c,
                                 source_a = "seed = 632741178", 
                                 source_b = "seed = 584925030", 
                                 source_c = "seed = 470267527")

# Calculating prevalence for all species in each dataset
all_prevalences_a <- calculate_prevalence(metaphlan.seed_632741178)$taxa_prevalence
sorted_species_a <- calculate_prevalence(metaphlan.seed_632741178)$top_taxa_prevalence

all_prevalences_b <- calculate_prevalence(metaphlan.seed_584925030)$taxa_prevalence
sorted_species_b <- calculate_prevalence(metaphlan.seed_584925030)$top_taxa_prevalence

all_prevalences_c <- calculate_prevalence(metaphlan.seed_470267527)$taxa_prevalence
sorted_species_c <- calculate_prevalence(metaphlan.seed_470267527)$top_taxa_prevalence

# Plotting barplot of prevalences 
prevalence_plot <- plot_prevalence(all_prevalences_a = all_prevalences_a,
                                   all_prevalences_b = all_prevalences_b,
                                   all_prevalences_c = all_prevalences_c,
                                   sorted_species_a = sorted_species_a,
                                   sorted_species_b = sorted_species_b,
                                   sorted_species_c = sorted_species_c,
                                   source_a = "seed = 632741178",
                                   source_b = "seed = 584925030",
                                   source_c = "seed = 470267527")

# -----------------------------------------------------------------------------
# Gather abundance and prevalence plots in one grid
# -----------------------------------------------------------------------------

plot <- (abundance_plot + prevalence_plot) + 
  plot_annotation(tag_levels = 'A') + 
  plot_layout(guides = "collect",
              axes = "collect") &
  theme(legend.position='bottom')
