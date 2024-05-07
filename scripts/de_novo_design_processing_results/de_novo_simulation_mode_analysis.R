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

# Processing and filtering CRCbiome MetaPhlAn data and preparing data for PCoA
metaphlan_crcbiome.tbl <- process_metaphlan_table_rel_abundance_filtering(metaphlan_path_crcbiome) %>% 
  process_metaphlan_table_pcoa("CRCbiome (n = 50)")

# Extracting a subset of 50 samples from the PCoA processed CRCbiome MetaPhlAn data
metaphlan_crcbiome_subset.tbl <- get_random_samples(metaphlan_crcbiome.tbl, 50)

# -----------------------------------------------------------------------------
# Collecting differential abundance simulation data
# -----------------------------------------------------------------------------

# Defining path to MetaPhlAn data
metaphlan_path_diff <- "/PATH/metaphlan_output/pilot_10_diff/tables/profiles.tsv"

# Processing and filtering MetaPhlAn data and preparing data for PCoA
metaphlan_diff.tbl <- process_metaphlan_table_rel_abundance_filtering(metaphlan_path_diff) %>% 
  process_metaphlan_table_pcoa("Differential abundance mode (n = 10)")

# Joining simulation data and CRCbiome data and performing PCoA
pcoa_diff <- joining_data(metaphlan_crcbiome_subset.tbl, metaphlan_diff.tbl) %>% 
  perform_pcoa_and_plot(method = "bray", log_mu = 1, log_sigma = 2)

plot_diff <- pcoa_diff$pcoa_plot + 
  labs(title = "", color = "") +
  scale_x_continuous(limits = c(-0.3, 0.75)) + coord_fixed()

# -----------------------------------------------------------------------------
# Collecting replicates simulation data
# -----------------------------------------------------------------------------

# Defining path to MetaPhlAn data
metaphlan_path_rep <- "/PATH/metaphlan_output/pilot_10_rep/tables/profiles.tsv"

# Processing and filtering MetaPhlAn data and preparing data for PCoA
metaphlan_rep.tbl <- process_metaphlan_table_rel_abundance_filtering(metaphlan_path_rep) %>% 
  process_metaphlan_table_pcoa("Replicates mode (n = 10)")

# Joining simulation data and CRCbiome data and performing PCoA
pcoa_rep <- joining_data(metaphlan_crcbiome_subset.tbl, metaphlan_rep.tbl) %>% 
  perform_pcoa_and_plot(method = "bray", log_mu = 1, log_sigma = 2)

plot_rep <- pcoa_rep$pcoa_plot + 
  labs(title = "", color = "") +
  scale_x_continuous(limits = c(-0.3, 0.75)) + coord_fixed()

# -----------------------------------------------------------------------------
# Gather plots in one grid
# -----------------------------------------------------------------------------

plot <- (plot_diff | plot_rep) + 
  plot_annotation(tag_levels = "A") + 
  plot_layout(axes = "collect") &
  theme(legend.position = "bottom")
