# Imports
library(tidyverse)
library(vegan)
library(glue)
library(patchwork)
library(ggpubr)
library(gridExtra)
library(grid)

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
# Preparing profiles with different adjustments and performing PCoA on the data
# -----------------------------------------------------------------------------

# 10 MAGs 2x
df.long_10_2x <- abundance_profile(coverage_df = median_coverage_genomes, 
                                   rand_mags_n = 10, 
                                   seed = 123, 
                                   adjustment_factor = 2.0)

pcoa_10_2x <- profiles_pcoa_plot(df.long = df.long_10_2x)

# 10 MAGs 5x
df.long_10_5x <- abundance_profile(coverage_df = median_coverage_genomes, 
                                   rand_mags_n = 10, 
                                   seed = 123, 
                                   adjustment_factor = 5.0)

pcoa_10_5x <- profiles_pcoa_plot(df.long = df.long_10_5x)

# 10 MAGs 10x
df.long_10_10x <- abundance_profile(coverage_df = median_coverage_genomes, 
                                    rand_mags_n = 10, 
                                    seed = 123, 
                                    adjustment_factor = 10.0)

pcoa_10_10x <- profiles_pcoa_plot(df.long = df.long_10_10x)

# 10 MAGs 20x
df.long_10_20x <- abundance_profile(coverage_df = median_coverage_genomes, 
                                    rand_mags_n = 10, 
                                    seed = 123, 
                                    adjustment_factor = 20.0)

pcoa_10_20x <- profiles_pcoa_plot(df.long = df.long_10_20x)

# 10 MAGs 25x
df.long_10_25x <- abundance_profile(coverage_df = median_coverage_genomes, 
                                    rand_mags_n = 10, 
                                    seed = 123, 
                                    adjustment_factor = 25.0)

pcoa_10_25x <- profiles_pcoa_plot(df.long = df.long_10_25x)

# 20 MAGs 2x
df.long_20_2x <- abundance_profile(coverage_df = median_coverage_genomes, 
                                   rand_mags_n = 20, 
                                   seed = 123, 
                                   adjustment_factor = 2.0)

pcoa_20_2x <- profiles_pcoa_plot(df.long = df.long_20_2x)

# 20 MAGs 5x
df.long_20_5x <- abundance_profile(coverage_df = median_coverage_genomes, 
                                   rand_mags_n = 20, 
                                   seed = 123, 
                                   adjustment_factor = 5.0)

pcoa_20_5x <- profiles_pcoa_plot(df.long = df.long_20_5x)

# 20 MAGs 10x
df.long_20_10x <- abundance_profile(coverage_df = median_coverage_genomes, 
                                    rand_mags_n = 20, 
                                    seed = 123, 
                                    adjustment_factor = 10.0)

pcoa_20_10x <- profiles_pcoa_plot(df.long = df.long_20_10x)

# 20 MAGs 20x
df.long_20_20x <- abundance_profile(coverage_df = median_coverage_genomes, 
                                    rand_mags_n = 20, 
                                    seed = 123, 
                                    adjustment_factor = 20.0)

pcoa_20_20x <- profiles_pcoa_plot(df.long = df.long_20_20x)

# 20 MAGs 25x
df.long_20_25x <- abundance_profile(coverage_df = median_coverage_genomes, 
                                    rand_mags_n = 20, 
                                    seed = 123, 
                                    adjustment_factor = 25.0)

pcoa_20_25x <- profiles_pcoa_plot(df.long = df.long_20_25x)

# 50 MAGs 2x
df.long_50_2x <- abundance_profile(coverage_df = median_coverage_genomes, 
                                   rand_mags_n = 50, 
                                   seed = 123, 
                                   adjustment_factor = 2.0)

pcoa_50_2x <- profiles_pcoa_plot(df.long = df.long_50_2x)

# 50 MAGs 5x
df.long_50_5x <- abundance_profile(coverage_df = median_coverage_genomes, 
                                   rand_mags_n = 50, 
                                   seed = 123, 
                                   adjustment_factor = 5.0)

pcoa_50_5x <- profiles_pcoa_plot(df.long = df.long_50_5x)

# 50 MAGs 10x
df.long_50_10x <- abundance_profile(coverage_df = median_coverage_genomes, 
                                    rand_mags_n = 50, 
                                    seed = 123, 
                                    adjustment_factor = 10.0)

pcoa_50_10x <- profiles_pcoa_plot(df.long = df.long_50_10x)

# 50 MAGs 20x
df.long_50_20x <- abundance_profile(coverage_df = median_coverage_genomes, 
                                    rand_mags_n = 50, 
                                    seed = 123, 
                                    adjustment_factor = 20.0)

pcoa_50_20x <- profiles_pcoa_plot(df.long = df.long_50_20x)

# 50 MAGs 25x
df.long_50_25x <- abundance_profile(coverage_df = median_coverage_genomes, 
                                    rand_mags_n = 50, 
                                    seed = 123, 
                                    adjustment_factor = 25.0)

pcoa_50_25x <- profiles_pcoa_plot(df.long = df.long_50_25x)

# 100 MAGs 2x
df.long_100_2x <- abundance_profile(coverage_df = median_coverage_genomes, 
                                    rand_mags_n = 100, 
                                    seed = 123, 
                                    adjustment_factor = 2.0)

pcoa_100_2x <- profiles_pcoa_plot(df.long = df.long_100_2x)

# 100 MAGs 5x
df.long_100_5x <- abundance_profile(coverage_df = median_coverage_genomes, 
                                    rand_mags_n = 100, 
                                    seed = 123, 
                                    adjustment_factor = 5.0)

pcoa_100_5x <- profiles_pcoa_plot(df.long = df.long_100_5x)

# 100 MAGs 10x
df.long_100_10x <- abundance_profile(coverage_df = median_coverage_genomes, 
                                     rand_mags_n = 100, 
                                     seed = 123, 
                                     adjustment_factor = 10.0)

pcoa_100_10x <- profiles_pcoa_plot(df.long = df.long_100_10x)

# 100 MAGs 20x
df.long_100_20x <- abundance_profile(coverage_df = median_coverage_genomes, 
                                     rand_mags_n = 100, 
                                     seed = 123, 
                                     adjustment_factor = 20.0)

pcoa_100_20x <- profiles_pcoa_plot(df.long = df.long_100_20x)

# 100 MAGs 25x
df.long_100_25x <- abundance_profile(coverage_df = median_coverage_genomes, 
                                     rand_mags_n = 100, 
                                     seed = 123, 
                                     adjustment_factor = 25.0)

pcoa_100_25x <- profiles_pcoa_plot(df.long = df.long_100_25x)

# -----------------------------------------------------------------------------
# PCoA plots
# -----------------------------------------------------------------------------

# 10 MAGs 2x
pcoa_plot_10_2x <- pcoa_10_2x$pcoa_plot + labs(title = "2x") 

# 10 MAGs 5x
pcoa_plot_10_5x <- pcoa_10_5x$pcoa_plot + labs(title = "5x") 

# 10 MAGs 10x
pcoa_plot_10_10x <- pcoa_10_10x$pcoa_plot + labs(title = "10x") 

# 10 MAGs 20x
pcoa_plot_10_20x <- pcoa_10_20x$pcoa_plot + labs(title = "20x") 

# 10 MAGs 25x
pcoa_plot_10_25x <- pcoa_10_25x$pcoa_plot + labs(title = "25x") 

# 20 MAGs 2x
pcoa_plot_20_2x <- pcoa_20_2x$pcoa_plot + labs(title = "") 

# 20 MAGs 5x
pcoa_plot_20_5x <- pcoa_20_5x$pcoa_plot + labs(title = "") 

# 20 MAGs 10x
pcoa_plot_20_10x <- pcoa_20_10x$pcoa_plot + labs(title = "") 

# 20 MAGs 20x
pcoa_plot_20_20x <- pcoa_20_20x$pcoa_plot + labs(title = "") 

# 20 MAGs 25x
pcoa_plot_20_25x <- pcoa_20_25x$pcoa_plot + labs(title = "") 

# 50 MAGs 2x
pcoa_plot_50_2x <- pcoa_50_2x$pcoa_plot + labs(title = "") 

# 50 MAGs 5x
pcoa_plot_50_5x <- pcoa_50_5x$pcoa_plot + labs(title = "") 

# 50 MAGs 10x
pcoa_plot_50_10x <- pcoa_50_10x$pcoa_plot + labs(title = "") 

# 50 MAGs 20x
pcoa_plot_50_20x <- pcoa_50_20x$pcoa_plot + labs(title = "") 

# 50 MAGs 25x
pcoa_plot_50_25x <- pcoa_50_25x$pcoa_plot + labs(title = "") 

# 100 MAGs 2x
pcoa_plot_100_2x <- pcoa_100_2x$pcoa_plot + labs(title = "") 

# 100 MAGs 5x
pcoa_plot_100_5x <- pcoa_100_5x$pcoa_plot + labs(title = "") 

# 100 MAGs 10x
pcoa_plot_100_10x <- pcoa_100_10x$pcoa_plot + labs(title = "") 

# 100 MAGs 20x
pcoa_plot_100_20x <- pcoa_100_20x$pcoa_plot + labs(title = "") 

# 100 MAGs 25x
pcoa_plot_100_25x <- pcoa_100_25x$pcoa_plot + labs(title = "") 

# -----------------------------------------------------------------------------
# Gather all PCoA plot in one grid
# -----------------------------------------------------------------------------

# PCoA score plots
plots <- list(pcoa_plot_10_2x, pcoa_plot_10_5x, pcoa_plot_10_10x, pcoa_plot_10_20x, pcoa_plot_10_25x,
              pcoa_plot_20_2x, pcoa_plot_20_5x, pcoa_plot_20_20x, pcoa_plot_20_20x, pcoa_plot_20_25x,
              pcoa_plot_50_2x, pcoa_plot_50_5x, pcoa_plot_50_10x, pcoa_plot_50_20x, pcoa_plot_50_25x,
              pcoa_plot_100_2x, pcoa_plot_100_5x, pcoa_plot_100_10x, pcoa_plot_100_20x, pcoa_plot_100_25x)

# Define a custom theme function that increases legend size
increase_legend_size <- function(size) {
  theme(legend.text = element_text(size = size),
        legend.title = element_text(size = size))
}

# Apply the custom theme function to each plot
plots <- lapply(plots, function(p) p + increase_legend_size(14))

# Arrange all plots in a grid
grid <- ggarrange(
  plotlist = plots,
  ncol = 5,
  nrow = 4,
  common.legend = TRUE,
  legend = "bottom"
)

# Add the annotations to the grid
annotated_grid <- annotate_figure(grid,
                                  left = grobTree(
                                    textGrob("10 MAGs", y=0.88, rot = 90, gp=gpar(fontface="bold", fontsize = 12)),
                                    textGrob("20 MAGs", y=0.65, rot = 90, gp=gpar(fontface="bold", fontsize = 12)),
                                    textGrob("50 MAGs", y=0.40, rot = 90, gp=gpar(fontface="bold", fontsize = 12)),
                                    textGrob("100 MAGs", y=0.16, rot = 90, gp=gpar(fontface="bold", fontsize = 12))
                                  )
)
