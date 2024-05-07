# Imports
library(tidyverse)
library(vegan)
library(glue)
library(mclust)
library(ggpubr)
library(grid)
library(patchwork)

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
# Preparing profiles with different adjustments
# -----------------------------------------------------------------------------

# 10 MAGs 2x
df.long_10_2x <- abundance_profile(coverage_df = median_coverage_genomes, 
                                   rand_mags_n = 10, 
                                   seed = 123, 
                                   adjustment_factor = 2.0)

# 20 MAGs 5x
df.long_20_5x <- abundance_profile(coverage_df = median_coverage_genomes, 
                                   rand_mags_n = 20, 
                                   seed = 123, 
                                   adjustment_factor = 5.0)

# 50 MAGs 10x
df.long_50_10x <- abundance_profile(coverage_df = median_coverage_genomes, 
                                    rand_mags_n = 50, 
                                    seed = 123, 
                                    adjustment_factor = 10.0)

# 100 MAGs 20x
df.long_100_20x <- abundance_profile(coverage_df = median_coverage_genomes, 
                                     rand_mags_n = 100, 
                                     seed = 123, 
                                     adjustment_factor = 20.0)

# -----------------------------------------------------------------------------
# Preparing data with real variables
# -----------------------------------------------------------------------------

# Sex data

data_by_sample.tbl <- readRDS("/PATH/datasets/data_by_FIT_sample/data_by_sample.Rds") %>% 
  select(sample_id, deltaker_id, Prøvetype, Total_Bases_QC_ATLAS) %>% 
  filter(Prøvetype == "Baseline") %>% 
  filter(Total_Bases_QC_ATLAS > 1e9)

screening_proc <- readRDS("/ess/p1068/data/durable/007-f_smei/001-trro/CRCbiome/datasets/screening/241122_screening_proc.rds") %>% 
  select(deltaker_id, kjonn)

df <- left_join(data_by_sample.tbl, screening_proc, by = "deltaker_id") %>% 
  select(sample_id, kjonn)

df.long_crcbiome_kjonn <- read_delim(mag_path, delim = "\t") %>%
  rename(sample_id = `...1`) %>% 
  mutate(sample_id = gsub("-", "_", sample_id, fixed = TRUE)) %>% 
  filter(sample_id %in% data_by_sample.tbl$sample_id) %>%  # Only include baseline samples
  left_join(df, by = "sample_id") %>%
  rename(group = kjonn) %>% 
  select(sample_id, group, everything()) %>% 
  # Convert data to long format
  pivot_longer(
    cols = -c(sample_id, group), # exclude the id and group from the gathering
    names_to = "MAG",
    values_to = "median_coverage") %>% 
  group_by(sample_id) %>%
  mutate(sum_measurements = sum(median_coverage, na.rm = TRUE)) %>% # Calculate the sum of measurements for each sample
  mutate(rel_abundance = (median_coverage/ sum_measurements)) %>% # Calculate the relative abundance by total reads and sum of measurements
  ungroup() %>% 
  select(-c(sum_measurements))

# Screening center data

data_by_sample.tbl <- readRDS("/PATH/datasets/data_by_FIT_sample/data_by_sample.Rds") %>% 
  select(sample_id, deltaker_id, Prøvetype, Total_Bases_QC_ATLAS) %>% 
  filter(Prøvetype == "Baseline") %>% 
  filter(Total_Bases_QC_ATLAS > 1e9)

screening_proc <- readRDS("/ess/p1068/data/durable/007-f_smei/001-trro/CRCbiome/datasets/screening/241122_screening_proc.rds") %>% 
  select(deltaker_id, senter)

df <- left_join(data_by_sample.tbl, screening_proc, by = "deltaker_id") %>% 
  select(sample_id, senter)

df.long_crcbiome_senter <- read_delim(mag_path, delim = "\t") %>%
  rename(sample_id = `...1`) %>% 
  mutate(sample_id = gsub("-", "_", sample_id, fixed = TRUE)) %>% 
  filter(sample_id %in% data_by_sample.tbl$sample_id) %>%  # Only include baseline samples
  left_join(df, by = "sample_id") %>% 
  rename(group = senter) %>% 
  select(sample_id, group, everything()) %>% 
  # Convert data to long format
  pivot_longer(
    cols = -c(sample_id, group), # exclude the id and group from the gathering
    names_to = "MAG",
    values_to = "median_coverage") %>% 
  group_by(sample_id) %>%
  mutate(sum_measurements = sum(median_coverage, na.rm = TRUE)) %>% # Calculate the sum of measurements for each sample
  mutate(rel_abundance = (median_coverage/ sum_measurements)) %>% # Calculate the relative abundance by total reads and sum of measurements
  ungroup() %>% 
  select(-c(sum_measurements))

# -----------------------------------------------------------------------------
# Performing PCoA and plotting PCoA plots
# -----------------------------------------------------------------------------

# 10 MAGs 2x
pcoa_10_2x <- profiles_pcoa_plot(df.long = df.long_10_2x)
pcoa_plot_10_2x <- pcoa_10_2x$pcoa_plot + labs(title = "10 MAGs 2x")

# 20 MAGs 5x
pcoa_20_5x <- profiles_pcoa_plot(df.long = df.long_20_5x)
pcoa_plot_20_5x <- pcoa_20_5x$pcoa_plot + labs(title = "20 MAGs 5x")

# 50 MAGs 10x
pcoa_50_10x <- profiles_pcoa_plot(df.long = df.long_50_10x)
pcoa_plot_50_10x <- pcoa_50_10x$pcoa_plot + labs(title = "50 MAGs 10x")

# 100 MAGs 20x
pcoa_100_20x <- profiles_pcoa_plot(df.long = df.long_100_20x)
pcoa_plot_100_20x <- pcoa_100_20x$pcoa_plot + labs(title = "100 MAGs 20x")

# CRCbiome sex
pcoa_crcbiome_kjonn <- profiles_pcoa_plot(df.long = df.long_crcbiome_kjonn)
pcoa_plot_crcbiome_kjonn <- pcoa_crcbiome_kjonn$pcoa_plot + labs(color = "Sex")

# CRCbiome screening center
pcoa_crcbiome_senter <- profiles_pcoa_plot(df.long = df.long_crcbiome_senter)
pcoa_plot_crcbiome_senter <- pcoa_crcbiome_senter$pcoa_plot + labs(color = "Screening center")

# -----------------------------------------------------------------------------
# Gathering all PCoA score plots in one grid
# -----------------------------------------------------------------------------

# PCoA plots simulated data
plots <- list(pcoa_plot_10_2x + scale_y_continuous(limits = c(-0.5, 0.5)) + scale_x_continuous(limits = c(-0.4, 0.4)) + coord_fixed(), 
              pcoa_plot_20_5x + scale_y_continuous(limits = c(-0.5, 0.5)) + scale_x_continuous(limits = c(-0.4, 0.4)) + coord_fixed(), 
              pcoa_plot_50_10x + scale_y_continuous(limits = c(-0.5, 0.5)) + scale_x_continuous(limits = c(-0.4, 0.4)) + coord_fixed(), 
              pcoa_plot_100_20x + scale_y_continuous(limits = c(-0.5, 0.5)) + scale_x_continuous(limits = c(-0.4, 0.4)) + coord_fixed())


# Define a custom theme function that increases legend size
increase_legend_size <- function(size) {
  theme(legend.text = element_text(size = size),
        legend.title = element_text(size = size))
}

# Apply the custom theme function to each plot
plots <- lapply(plots, function(p) p + increase_legend_size(14))

# Arrange all plots
grid <- ggarrange(
  plotlist = plots,
  ncol = 4, 
  nrow = 1,
  common.legend = TRUE,
  legend = "bottom"
)

# PcoA plots real variables

sex_plot <- pcoa_plot_crcbiome_kjonn + 
  scale_y_continuous(limits = c(-0.5, 0.5)) + 
  scale_x_continuous(limits = c(-0.4, 0.4)) + 
  coord_fixed() +
  theme(
    text = element_text(size = 15, face = "bold"),  
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 14),  
    legend.title = element_text(size = 14),  
    legend.text = element_text(size = 14),  
    plot.title = element_text(size = 16, face = "bold"),  
    strip.text.x = element_text(size = 14)
  )

center_plot <- pcoa_plot_crcbiome_senter + 
  scale_y_continuous(limits = c(-0.5, 0.5)) + 
  scale_x_continuous(limits = c(-0.4, 0.4)) + 
  coord_fixed() +
  theme(
    text = element_text(size = 15, face = "bold"),  
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 14),  
    legend.title = element_text(size = 14),  
    legend.text = element_text(size = 14),  
    plot.title = element_text(size = 16, face = "bold"),  
    strip.text.x = element_text(size = 14)
  )

p <- sex_plot + center_plot + plot_annotation(tag_levels = 'A')

# -----------------------------------------------------------------------------
# Finding the number of PCoA components to work with in the clustering for 
# each dataset 
# HC with average linkage with bootstrapping
# -----------------------------------------------------------------------------

# Datasets to investigate
abundance_data_list <- list(df.long_10_2x,
                            df.long_20_5x,
                            df.long_50_10x,
                            df.long_100_20x,
                            df.long_crcbiome_kjonn,
                            df.long_crcbiome_senter)

n_bootstraps <- 50
dimension_range <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 50)

# Initialize an empty data frame to store results
results_df_HC_avg_PCoA <- data.frame(Dataset = integer(),
                                     Dimensions = integer(),
                                     Mean_Variance_Explained = integer(),
                                     Bootstrapped_ARIs = I(list()), 
                                     p_values = I(list()),
                                     Mean_Bootstrapped_ARI = numeric(),
                                     SD_Bootstrapped_ARI = numeric())

set.seed(123)
for (dataset_index in 1:length(abundance_data_list)) {
  
  wide_data <- abundance_data_list[[dataset_index]] %>% 
    pivot_wider(names_from = MAG, 
                values_from = rel_abundance, 
                id_cols = c(sample_id, group))
  
  # Initialize a matrix to store the bootstrapped ARI values for each dimension
  bootstrapped_ari_matrix <- matrix(NA, nrow = n_bootstraps, ncol = length(dimension_range))
  
  # Initialize a matrix to store the p_values for each dimension
  p_value_ari_matrix <- matrix(NA, nrow = n_bootstraps, ncol = length(dimension_range))
  
  # Initialize a list where each element holds the variances explained for each dimension.
  variance_explained_sums <- vector("list", length(dimension_range))
  for (i in seq_along(dimension_range)) {
    variance_explained_sums[[i]] <- numeric(n_bootstraps)
  }
  
  for (i_boot in 1:n_bootstraps) {
    
    df <- wide_data %>% 
      slice_sample(n = nrow(.), replace = TRUE)
    
    pcoa_res <- df %>% 
      select(-c(sample_id, group)) %>% 
      vegdist(method = "bray") %>% 
      cmdscale(k = max(dimension_range), eig = TRUE)
    
    # Extract coordinates and set column names
    coordinates <- as.data.frame(pcoa_res$points)
    colnames(coordinates) <- paste("PCo", seq_len(ncol(coordinates)), sep = "")
    
    # Extract explained variation for each axis
    percent_explained <- round(100 * pcoa_res$eig / sum(pcoa_res$eig), 2)
    
    coordinates_with_groups <- bind_cols(group = df$group, coordinates)
    
    true_clusters <- as.integer(coordinates_with_groups$group)
    
    for (num_dimensions in dimension_range) {

      hcl.obj <- coordinates_with_groups[, 1:(num_dimensions + 1)] %>% 
        select(-c(group)) %>% 
        dist(method = "euclidean") %>% 
        hclust(method = "average")
      
      # Choose an appropriate number of clusters for cutting the dendrogram
      predicted_clusters <- cutree(hcl.obj, k = 2)
      
      # Compute ARI value
      ari_val <- as.numeric(adjRand_test(true_clusters, predicted_clusters, perm = 1000)[1])
      
      p_val <- as.numeric(adjRand_test(true_clusters, predicted_clusters, perm = 1000)[2])
      
      # Find the position of the current_dimension in the dimension_range
      dimension_pos <- match(num_dimensions, dimension_range)
      
      # Store the ARI value in the matrix
      bootstrapped_ari_matrix[i_boot, dimension_pos] <- ari_val
      
      # Store the p value in the matrix
      p_value_ari_matrix[i_boot, dimension_pos] <- p_val
      
      variance_explained_sums[[which(dimension_range == num_dimensions)]][i_boot] <- sum(percent_explained[1:num_dimensions])
    }
  }
  
  # Now, compute the mean and SD for each dimension across all bootstraps
  mean_ari_per_dimension <- apply(bootstrapped_ari_matrix, 2, mean, na.rm = TRUE)
  sd_ari_per_dimension <- apply(bootstrapped_ari_matrix, 2, sd, na.rm = TRUE)
  mean_variance_explained <- sapply(variance_explained_sums, mean, na.rm = TRUE)
  
  # Append results to the data frame
  results_df_HC_avg_PCoA <- rbind(results_df_HC_avg_PCoA, data.frame(Dataset = dataset_index,
                                                                     Dimensions = dimension_range,
                                                                     Mean_Bootstrapped_ARI = mean_ari_per_dimension,
                                                                     SD_Bootstrapped_ARI = sd_ari_per_dimension,
                                                                     Bootstrapped_ARIs = I(split(bootstrapped_ari_matrix, col(bootstrapped_ari_matrix))),
                                                                     Mean_Variance_Explained = mean_variance_explained,
                                                                     p_values = I(split(p_value_ari_matrix, col(p_value_ari_matrix)))))
}

# Modify dataset labels in result_df
dataset_labels <- c("10 MAGs 2x", "20 MAGs 5x", "50 MAGs 10x", "100 MAGs 20x", "Sex", "Screening center", "Antibiotics")
results_df_HC_avg_PCoA$Dataset <- factor(results_df_HC_avg_PCoA$Dataset, 
                                         levels = 1:length(dataset_labels), 
                                         labels = dataset_labels)

# Count the number of significant p-values
results_df_HC_avg_PCoA$Count_P_Values_Below_0_01 <- sapply(results_df_HC_avg_PCoA$p_values, function(p_vec) {
  sum(p_vec < 0.01, na.rm = TRUE)
})

# Save results for later use
# saveRDS(results_df_HC_avg_PCoA,
#         "/PATH/scripts/profile_based_design_processing_results/data/results_df_pcoa_loop_HC_avg_bootstrap_50.rds")

# Read results
# results_df_HC_avg_PCoA <- readRDS("/PATH/scripts/profile_based_design_processing_results/data/results_df_pcoa_loop_HC_avg_bootstrap_50.rds")

# Finding the maximum ARI for each dataset
max_ari_per_dataset <- results_df_HC_avg_PCoA %>%
  group_by(Dataset) %>%
  summarize(Max_ARI = max(Mean_Bootstrapped_ARI))

# Join this back with the original data to get a flag for the max ARI
results_df_HC_avg_PCoA <- results_df_HC_avg_PCoA %>%
  left_join(max_ari_per_dataset, by = "Dataset") %>%
  mutate(Is_Max_ARI = Mean_Bootstrapped_ARI == Max_ARI)

# Finding which dimension has the maximum ARI for each dataset
max_ari_df <- results_df_HC_avg_PCoA %>% 
  filter(Is_Max_ARI == TRUE) %>% 
  select(Dataset, Dimensions, Mean_Bootstrapped_ARI)

# Plotting ARI score across the dimensions. Maximum ARI is marked in the plot
pcoa_looping_plot <- ggplot(results_df_HC_avg_PCoA, aes(x = Dimensions, y = Mean_Bootstrapped_ARI, group = Dataset, color = Dataset)) +
  geom_line(linewidth = 1) +
  geom_point(data = filter(results_df_HC_avg_PCoA, Is_Max_ARI), aes(x = Dimensions, y = Mean_Bootstrapped_ARI, color = Dataset), size = 3, shape = 19) +
  theme_bw() +
  labs(
    x = "Number of PCoA Dimensions",
    y = "Mean Adjusted Rand Index",
    color = "Dataset"
  ) +
  guides(color = guide_legend(override.aes = list(size=3))) +
  theme(
    text = element_text(size = 12, face = "bold"), 
    axis.title = element_text(size = 12), 
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12), 
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

# -----------------------------------------------------------------------------
# Finding the number of PCoA components to work with in the clustering 
# for each dataset 
# HC with complete linkage with bootstrapping
# -----------------------------------------------------------------------------

# Datasets to investigate
abundance_data_list <- list(df.long_10_2x,
                            df.long_20_5x,
                            df.long_50_10x,
                            df.long_100_20x,
                            df.long_crcbiome_kjonn,
                            df.long_crcbiome_senter)

n_bootstraps <- 50
dimension_range <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 50)

# Initialize an empty data frame to store results
results_df_HC_complete_PCoA <- data.frame(Dataset = integer(),
                                          Dimensions = integer(),
                                          Mean_Variance_Explained = integer(),
                                          Bootstrapped_ARIs = I(list()), 
                                          p_values = I(list()),
                                          Mean_Bootstrapped_ARI = numeric(),
                                          SD_Bootstrapped_ARI = numeric())

set.seed(123)
for (dataset_index in 1:length(abundance_data_list)) {
  
  wide_data <- abundance_data_list[[dataset_index]] %>% 
    pivot_wider(names_from = MAG, 
                values_from = rel_abundance, 
                id_cols = c(sample_id, group))
  
  # Initialize a matrix to store the bootstrapped ARI values for each dimension
  bootstrapped_ari_matrix <- matrix(NA, nrow = n_bootstraps, ncol = length(dimension_range))
  
  # Initialize a matrix to store the p_values for each dimension
  p_value_ari_matrix <- matrix(NA, nrow = n_bootstraps, ncol = length(dimension_range))
  
  # Initialize empty list where each element holds the variances explained for each dimension.
  variance_explained_sums <- vector("list", length(dimension_range))
  for (i in seq_along(dimension_range)) {
    variance_explained_sums[[i]] <- numeric(n_bootstraps)
  }
  
  for (i_boot in 1:n_bootstraps) {
    
    df <- wide_data %>% 
      slice_sample(n = nrow(.), replace = TRUE)
    
    pcoa_res <- df %>% 
      select(-c(sample_id, group)) %>% 
      vegdist(method = "bray") %>% 
      cmdscale(k = max(dimension_range), eig = TRUE)
    
    # Extract coordinates and set column names
    coordinates <- as.data.frame(pcoa_res$points)
    colnames(coordinates) <- paste("PCo", seq_len(ncol(coordinates)), sep = "")
    
    # Extract explained variation for each axis
    percent_explained <- round(100 * pcoa_res$eig / sum(pcoa_res$eig), 2)
    
    coordinates_with_groups <- bind_cols(group = df$group, coordinates)
    
    true_clusters <- as.integer(coordinates_with_groups$group)
    
    for (num_dimensions in dimension_range) {
      
      hcl.obj <- coordinates_with_groups[, 1:(num_dimensions + 1)] %>% 
        select(-c(group)) %>% 
        dist(method = "euclidean") %>% 
        hclust(method = "complete")
      
      # Choose an appropriate number of clusters for cutting the dendrogram
      predicted_clusters <- cutree(hcl.obj, k = 2)
      
      # Compute ARI value
      ari_val <- as.numeric(adjRand_test(true_clusters, predicted_clusters, perm = 1000)[1])
      
      p_val <- as.numeric(adjRand_test(true_clusters, predicted_clusters, perm = 1000)[2])
      
      # Find the position of the current_dimension in the dimension_range
      dimension_pos <- match(num_dimensions, dimension_range)
      
      # Store the ARI value in the matrix
      bootstrapped_ari_matrix[i_boot, dimension_pos] <- ari_val
      
      # Store the p value in the matrix
      p_value_ari_matrix[i_boot, dimension_pos] <- p_val
      
      variance_explained_sums[[which(dimension_range == num_dimensions)]][i_boot] <- sum(percent_explained[1:num_dimensions])
    }
  }
  
  # Compute the mean and SD for each dimension across all bootstraps
  mean_ari_per_dimension <- apply(bootstrapped_ari_matrix, 2, mean, na.rm = TRUE)
  sd_ari_per_dimension <- apply(bootstrapped_ari_matrix, 2, sd, na.rm = TRUE)
  mean_variance_explained <- sapply(variance_explained_sums, mean, na.rm = TRUE)
  
  # Append results to the data frame
  results_df_HC_complete_PCoA <- rbind(results_df_HC_complete_PCoA, data.frame(Dataset = dataset_index,
                                                                               Dimensions = dimension_range,
                                                                               Mean_Bootstrapped_ARI = mean_ari_per_dimension,
                                                                               SD_Bootstrapped_ARI = sd_ari_per_dimension,
                                                                               Bootstrapped_ARIs = I(split(bootstrapped_ari_matrix, col(bootstrapped_ari_matrix))),
                                                                               Mean_Variance_Explained = mean_variance_explained,
                                                                               p_values = I(split(p_value_ari_matrix, col(p_value_ari_matrix)))))
}

# Modify dataset labels in result_df
dataset_labels <- c("10 MAGs 2x", "20 MAGs 5x", "50 MAGs 10x", "100 MAGs 20x", "Sex", "Screening center")
results_df_HC_complete_PCoA$Dataset <- factor(results_df_HC_complete_PCoA$Dataset, 
                                              levels = 1:length(dataset_labels), 
                                              labels = dataset_labels)

# Count the number of significant p-values
results_df_HC_complete_PCoA$Count_P_Values_Below_0_01 <- sapply(results_df_HC_complete_PCoA$p_values, function(p_vec) {
  sum(p_vec < 0.01, na.rm = TRUE)
})

# Save results for later use
# saveRDS(results_df_HC_complete_PCoA,
#         "/PATH/scripts/profile_based_design_processing_results/data/results_df_pcoa_loop_HC_complete_bootstrap_50.rds")

# Read results
# results_df_HC_complete_PCoA <- readRDS("/PATH/scripts/profile_based_design_processing_results/data/results_df_pcoa_loop_HC_complete_bootstrap_50.rds")

# Finding the maximum ARI for each dataset
max_ari_per_dataset <- results_df_HC_complete_PCoA %>%
  group_by(Dataset) %>%
  summarize(Max_ARI = max(Mean_Bootstrapped_ARI))

# Join this back with the original data to get a flag for the max ARI
results_df_HC_complete_PCoA <- results_df_HC_complete_PCoA %>%
  left_join(max_ari_per_dataset, by = "Dataset") %>%
  mutate(Is_Max_ARI = Mean_Bootstrapped_ARI == Max_ARI)

# Finding which dimension has the maximum ARI for each dataset
max_ari_df <- results_df_HC_complete_PCoA %>% 
  filter(Is_Max_ARI == TRUE) %>% 
  select(Dataset, Dimensions, Mean_Bootstrapped_ARI)

# Plotting ARI score across the dimensions. Maximum ARI is marked in the plot
pcoa_looping_plot <- ggplot(results_df_HC_complete_PCoA, aes(x = Dimensions, y = Mean_Bootstrapped_ARI, group = Dataset, color = Dataset)) +
  geom_line(linewidth = 1) +
  geom_point(data = filter(results_df_HC_complete_PCoA, Is_Max_ARI), aes(x = Dimensions, y = Mean_Bootstrapped_ARI, color = Dataset), size = 3, shape = 19) +
  theme_bw() +
  labs(
    x = "Number of PCoA Dimensions",
    y = "Mean Adjusted Rand Index",
    color = "Dataset"
  ) +
  guides(color = guide_legend(override.aes = list(size=3))) +
  theme(
    text = element_text(size = 12, face = "bold"),  
    axis.title = element_text(size = 12), 
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12), 
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

# -----------------------------------------------------------------------------
# Finding the number of PCoA components to work with in the clustering for 
# each dataset 
# K means clustering with bootstrapping
# ------------------------------------------------------------------------------

# Datasets to investigate
abundance_data_list <- list(df.long_10_2x,
                            df.long_20_5x,
                            df.long_50_10x,
                            df.long_100_20x,
                            df.long_crcbiome_kjonn,
                            df.long_crcbiome_senter)

n_bootstraps <- 50
dimension_range <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 50)

# Initialize an empty data frame to store results
results_df_Kmeans_PCoA <- data.frame(Dataset = integer(),
                                     Dimensions = integer(),
                                     Mean_Variance_Explained = integer(),
                                     Bootstrapped_ARIs = I(list()), 
                                     p_values = I(list()),
                                     Mean_Bootstrapped_ARI = numeric(),
                                     SD_Bootstrapped_ARI = numeric())

set.seed(123)
for (dataset_index in 1:length(abundance_data_list)) {
  
  wide_data <- abundance_data_list[[dataset_index]] %>% 
    pivot_wider(names_from = MAG, 
                values_from = rel_abundance, 
                id_cols = c(sample_id, group))
  
  # Initialize a matrix to store the bootstrapped ARI values for each dimension
  bootstrapped_ari_matrix <- matrix(NA, nrow = n_bootstraps, ncol = length(dimension_range))
  
  # Initialize a matrix to store the p_values for each dimension
  p_value_ari_matrix <- matrix(NA, nrow = n_bootstraps, ncol = length(dimension_range))
  
  # Initialize empty list where each element holds the variances explained for each dimension.
  variance_explained_sums <- vector("list", length(dimension_range))
  for (i in seq_along(dimension_range)) {
    variance_explained_sums[[i]] <- numeric(n_bootstraps)
  }
  
  for (i_boot in 1:n_bootstraps) {
    
    df <- wide_data %>% 
      slice_sample(n = nrow(.), replace = TRUE)
    
    pcoa_res <- df %>% 
      select(-c(sample_id, group)) %>% 
      vegdist(method = "bray") %>% 
      cmdscale(k = max(dimension_range), eig = TRUE)
    
    # Extract coordinates and set column names
    coordinates <- as.data.frame(pcoa_res$points)
    colnames(coordinates) <- paste("PCo", seq_len(ncol(coordinates)), sep = "")
    
    # Extract explained variation for each axis
    percent_explained <- round(100 * pcoa_res$eig / sum(pcoa_res$eig), 2)
    
    coordinates_with_groups <- bind_cols(group = df$group, coordinates)
    
    true_clusters <- as.integer(coordinates_with_groups$group)
    
    for (num_dimensions in dimension_range) {
      
      # Perform k means clustering
      kmeans.obj <- coordinates_with_groups[, 1:(num_dimensions + 1)] %>% 
        select(-group) %>% 
        kmeans(centers = 2, 
               iter.max = 10, 
               nstart = 10)
      
      predicted_clusters <- kmeans.obj$cluster
      
      # Compute ARI value
      ari_val <- as.numeric(adjRand_test(true_clusters, predicted_clusters, perm = 1000)[1])
      
      p_val <- as.numeric(adjRand_test(true_clusters, predicted_clusters, perm = 1000)[2])
      
      # Find the position of the current_dimension in the dimension_range
      dimension_pos <- match(num_dimensions, dimension_range)
      
      # Store the ARI value in the matrix
      bootstrapped_ari_matrix[i_boot, dimension_pos] <- ari_val
      
      # Store the p value in the matrix
      p_value_ari_matrix[i_boot, dimension_pos] <- p_val
      
      variance_explained_sums[[which(dimension_range == num_dimensions)]][i_boot] <- sum(percent_explained[1:num_dimensions])
    }
  }
  
  # Compute the mean and SD for each dimension across all bootstraps
  mean_ari_per_dimension <- apply(bootstrapped_ari_matrix, 2, mean, na.rm = TRUE)
  sd_ari_per_dimension <- apply(bootstrapped_ari_matrix, 2, sd, na.rm = TRUE)
  mean_variance_explained <- sapply(variance_explained_sums, mean, na.rm = TRUE)
  
  # Append results to the data frame
  results_df_Kmeans_PCoA <- rbind(results_df_Kmeans_PCoA, data.frame(Dataset = dataset_index,
                                                                     Dimensions = dimension_range,
                                                                     Mean_Bootstrapped_ARI = mean_ari_per_dimension,
                                                                     SD_Bootstrapped_ARI = sd_ari_per_dimension,
                                                                     Bootstrapped_ARIs = I(split(bootstrapped_ari_matrix, col(bootstrapped_ari_matrix))),
                                                                     Mean_Variance_Explained = mean_variance_explained,
                                                                     p_values = I(split(p_value_ari_matrix, col(p_value_ari_matrix)))))
}

# Modify dataset labels in result_df
dataset_labels <- c("10 MAGs 2x", "20 MAGs 5x", "50 MAGs 10x", "100 MAGs 20x", "Sex", "Screening center")
results_df_Kmeans_PCoA$Dataset <- factor(results_df_Kmeans_PCoA$Dataset, 
                             levels = 1:length(dataset_labels), 
                             labels = dataset_labels)

# Count the number of significant p-values
results_df_Kmeans_PCoA$Count_P_Values_Below_0_01 <- sapply(results_df_Kmeans_PCoA$p_values, function(p_vec) {
  sum(p_vec < 0.01, na.rm = TRUE)
})

# Save results for later use
# saveRDS(results_df_Kmeans_PCoA,
#         "/PATH/scripts/profile_based_design_processing_results/data/results_df_pcoa_loop_kmeans_bootstrap_50.rds")

# Read results
# results_df_Kmeans_PCoA <- readRDS("/PATH/scripts/profile_based_design_processing_results/data/results_df_pcoa_loop_kmeans_bootstrap_50.rds")

# Finding the maximum ARI for each dataset
max_ari_per_dataset <- results_df_Kmeans_PCoA %>%
  group_by(Dataset) %>%
  summarize(Max_ARI = max(Mean_Bootstrapped_ARI))

# Join this back with the original data to get a flag for the max ARI
results_df_Kmeans_PCoA <- results_df_Kmeans_PCoA %>%
  left_join(max_ari_per_dataset, by = "Dataset") %>%
  mutate(Is_Max_ARI = Mean_Bootstrapped_ARI == Max_ARI)

# Finding which dimension has the maximum ARI for each dataset
max_ari_df <- results_df_Kmeans_PCoA %>% 
  filter(Is_Max_ARI == TRUE) %>% 
  select(Dataset, Dimensions, Mean_Bootstrapped_ARI)

# Plotting ARI score across the dimensions. Maximum ARI is marked in the plot
pcoa_looping_plot <- ggplot(results_df_Kmeans_PCoA, aes(x = Dimensions, y = Mean_Bootstrapped_ARI, group = Dataset, color = Dataset)) +
  geom_line(linewidth = 1) +
  geom_point(data = filter(results_df_Kmeans_PCoA, Is_Max_ARI), aes(x = Dimensions, y = Mean_Bootstrapped_ARI, color = Dataset), size = 3, shape = 19) +
  theme_bw() +
  labs(
    x = "Number of PCoA Dimensions",
    y = "Mean Adjusted Rand Index",
    color = "Dataset"
  ) +
  guides(color = guide_legend(override.aes = list(size=3))) +
  theme(
    text = element_text(size = 12, face = "bold"),  
    axis.title = element_text(size = 12), 
    axis.text = element_text(size = 12), 
    legend.title = element_text(size = 12),  
    legend.text = element_text(size = 12),  
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

# -----------------------------------------------------------------------------
# Gather all PCoA scree plots into one grid
# -----------------------------------------------------------------------------

scree_plot_10_2x <- pcoa_10_2x$scree_plot + labs(title = "10 MAGs 2x") + xlim(0, 20 + 1) + scale_y_continuous(limits = c(0, 15)) + coord_fixed()
scree_plot_20_5x <- pcoa_20_5x$scree_plot + labs(title = "20 MAGs 5x") + xlim(0, 20 + 1) + scale_y_continuous(limits = c(0, 15)) + coord_fixed()
scree_plot_50_10x <- pcoa_50_10x$scree_plot + labs(title = "50 MAGs 10x") + xlim(0, 20 + 1) + scale_y_continuous(limits = c(0, 15)) + coord_fixed()
scree_plot_100_20x <- pcoa_100_20x$scree_plot + labs(title = "100 MAGs 20x") + xlim(0, 20 + 1) + scale_y_continuous(limits = c(0, 15)) + coord_fixed()

# Grid scree plots
plots <- list(scree_plot_10_2x, scree_plot_20_5x, scree_plot_50_10x, scree_plot_100_20x)

# Arrange all plots
grid <- ggarrange(
  plotlist = plots,
  ncol = 4, 
  nrow = 1,
  common.legend = TRUE,
  legend = "bottom"
)

# Add the annotations to the combined plot
annotated_grid <- annotate_figure(grid,
                                  top = text_grob("Scree plot for selected data", face="bold", size = 20))

# -----------------------------------------------------------------------------
# Gather score plots and scree plots into one grid
# -----------------------------------------------------------------------------

scree_plot_10_2x <- pcoa_10_2x$scree_plot + labs(title = "") + xlim(0, 20 + 1) + scale_y_continuous(limits = c(0, 15)) + coord_fixed()
scree_plot_20_5x <- pcoa_20_5x$scree_plot + labs(title = "") + xlim(0, 20 + 1) + scale_y_continuous(limits = c(0, 15)) + coord_fixed()
scree_plot_50_10x <- pcoa_50_10x$scree_plot + labs(title = "") + xlim(0, 20 + 1) + scale_y_continuous(limits = c(0, 15)) + coord_fixed()
scree_plot_100_20x <- pcoa_100_20x$scree_plot + labs(title = "") + xlim(0, 20 + 1) + scale_y_continuous(limits = c(0, 15)) + coord_fixed()

plots <- list(pcoa_plot_10_2x + scale_y_continuous(limits = c(-0.5, 0.5)) + scale_x_continuous(limits = c(-0.4, 0.4)) + coord_fixed(), 
              pcoa_plot_20_5x + scale_y_continuous(limits = c(-0.5, 0.5)) + scale_x_continuous(limits = c(-0.4, 0.4)) + coord_fixed(), 
              pcoa_plot_50_10x + scale_y_continuous(limits = c(-0.5, 0.5)) + scale_x_continuous(limits = c(-0.4, 0.4)) + coord_fixed(), 
              pcoa_plot_100_20x + scale_y_continuous(limits = c(-0.5, 0.5)) + scale_x_continuous(limits = c(-0.4, 0.4)) + coord_fixed(),
              scree_plot_10_2x,
              scree_plot_20_5x,
              scree_plot_50_10x,
              scree_plot_100_20x)

# Define a custom theme function that increases legend size
increase_legend_size <- function(size) {
  theme(legend.text = element_text(size = size),
        legend.title = element_text(size = size))
}

# Apply the custom theme function to each plot
plots <- lapply(plots, function(p) p + increase_legend_size(14))

# Arrange all plots
grid <- ggarrange(
  plotlist = plots,
  ncol = 4, 
  nrow = 2,
  common.legend = TRUE,
  legend = "bottom",
  labels = "AUTO"
)

# -----------------------------------------------------------------------------
# Hierarchical clustering directly on the abundance profiles 
# (no dimension reduction)
# HC average linkage
# -----------------------------------------------------------------------------

# Datasets to investigate
abundance_data_list <- list(df.long_10_2x,
                            df.long_20_5x,
                            df.long_50_10x,
                            df.long_100_20x,
                            df.long_crcbiome_kjonn,
                            df.long_crcbiome_senter)

n_bootstraps <- 50

# Initialize an empty data frame to store results
results_df_HC_avg <- data.frame(Dataset = integer(),
                                Bootstrapped_ARIs = I(list()), 
                                p_values = I(list()),
                                Mean_Bootstrapped_ARI = numeric(),
                                SD_Bootstrapped_ARI = numeric())

# Initialize a matrix to store the bootstrapped ARI values for each dataset
bootstrapped_ari_matrix <- matrix(NA, nrow = n_bootstraps, ncol = length(abundance_data_list))

# Initialize a matrix to store the bootstrapped p values for each dataset
p_value_ari_matrix <- matrix(NA, nrow = n_bootstraps, ncol = length(abundance_data_list))


set.seed(123)
for (dataset_index in 1:length(abundance_data_list)) {
  
  wide_data <- abundance_data_list[[dataset_index]] %>% 
    pivot_wider(names_from = MAG, 
                values_from = rel_abundance, 
                id_cols = c(sample_id, group))

  for (i_boot in 1:n_bootstraps) {
    
    df <- wide_data %>% 
      slice_sample(n = nrow(.), replace = TRUE)
    
    true_clusters <- as.integer(df$group)
      
    hcl.obj <- df %>% 
      select(-c(sample_id, group)) %>% 
      dist(method = "euclidean") %>% 
      hclust(method = "average")
    
    # Choose an appropriate number of clusters for cutting the dendrogram
    predicted_clusters <- cutree(hcl.obj, k = 2)
    
    # Compute ARI value
    ari_val <- as.numeric(adjRand_test(true_clusters, predicted_clusters, perm = 1000)[1])

    p_val <- as.numeric(adjRand_test(true_clusters, predicted_clusters, perm = 1000)[2])
    
    # Store the ARI value in the matrix
    bootstrapped_ari_matrix[i_boot, dataset_index] <- ari_val

    # Store the p value in the matrix
    p_value_ari_matrix[i_boot, dataset_index] <- p_val
  }
  
  # Calculate the mean ARI for the current dataset across all bootstraps
  mean_ari <- apply(bootstrapped_ari_matrix[, dataset_index, drop = FALSE], 2, mean, na.rm = TRUE)
  
  # Calculate the standard deviation of ARI for the current dataset across all bootstraps
  sd_ari <- apply(bootstrapped_ari_matrix[, dataset_index, drop = FALSE], 2, sd, na.rm = TRUE)
  
  # Append results to the data frame
  results_df_HC_avg <- rbind(results_df_HC_avg, data.frame(Dataset = dataset_index,
                                                           Mean_Bootstrapped_ARI = mean_ari,
                                                           SD_Bootstrapped_ARI = sd_ari,
                                                           Bootstrapped_ARIs = I(list(bootstrapped_ari_matrix[, dataset_index])),
                                                           p_values = I(list(p_value_ari_matrix[, dataset_index]))
  ))
}

# Modify dataset labels in result_df
dataset_labels <- c("10 MAGs 2x", "20 MAGs 5x", "50 MAGs 10x", "100 MAGs 20x", "Sex", "Screening center")
results_df_HC_avg$Dataset <- factor(results_df_HC_avg$Dataset, 
                                    levels = 1:length(dataset_labels), 
                                    labels = dataset_labels)

# Count the number of significant p-values
results_df_HC_avg$Count_P_Values_Below_0_01 <- sapply(results_df_HC_avg$p_values, function(p_vec) {
  sum(p_vec < 0.01, na.rm = TRUE)
})

# Save results for later use
# saveRDS(results_df_HC_avg,
#         "/PATH/scripts/profile_based_design_processing_results/data/results_df_HC_avg_bootstrap_50_wo_pcoa.rds")

# Read results
# results_df_HC_avg <- readRDS("/PATH/scripts/profile_based_design_processing_results/data/results_df_HC_avg_bootstrap_50_wo_pcoa.rds")

# -----------------------------------------------------------------------------
# Hierarchical clustering directly on the abundance profiles 
# (no dimension reduction)
# HC complete linkage
# -----------------------------------------------------------------------------

# Datasets to investigate
abundance_data_list <- list(df.long_10_2x,
                            df.long_20_5x,
                            df.long_50_10x,
                            df.long_100_20x,
                            df.long_crcbiome_kjonn,
                            df.long_crcbiome_senter)

n_bootstraps <- 50

# Initialize an empty data frame to store results
results_df_HC_complete <- data.frame(Dataset = integer(),
                                     Bootstrapped_ARIs = I(list()),  
                                     p_values = I(list()),
                                     Mean_Bootstrapped_ARI = numeric(),
                                     SD_Bootstrapped_ARI = numeric())

# Initialize a matrix to store the bootstrapped ARI values for each dataset
bootstrapped_ari_matrix <- matrix(NA, nrow = n_bootstraps, ncol = length(abundance_data_list))

# Initialize a matrix to store the bootstrapped p values for each dataset
p_value_ari_matrix <- matrix(NA, nrow = n_bootstraps, ncol = length(abundance_data_list))


set.seed(123)
for (dataset_index in 1:length(abundance_data_list)) {
  
  wide_data <- abundance_data_list[[dataset_index]] %>% 
    pivot_wider(names_from = MAG, 
                values_from = rel_abundance, 
                id_cols = c(sample_id, group))
  
  for (i_boot in 1:n_bootstraps) {
    
    df <- wide_data %>% 
      slice_sample(n = nrow(.), replace = TRUE)
    
    true_clusters <- as.integer(df$group)
    
    hcl.obj <- df %>% 
      select(-c(sample_id, group)) %>% 
      dist(method = "euclidean") %>% 
      hclust(method = "complete")
    
    # Choose an appropriate number of clusters for cutting the dendrogram
    predicted_clusters <- cutree(hcl.obj, k = 2)
    
    # Compute ARI value
    ari_val <- as.numeric(adjRand_test(true_clusters, predicted_clusters, perm = 1000)[1])

    p_val <- as.numeric(adjRand_test(true_clusters, predicted_clusters, perm = 1000)[2])
    
    # Store the ARI value in the matrix
    bootstrapped_ari_matrix[i_boot, dataset_index] <- ari_val

    # Store the p value in the matrix
    p_value_ari_matrix[i_boot, dataset_index] <- p_val
  }
  
  # Calculate the mean ARI for the current dataset across all bootstraps
  mean_ari <- apply(bootstrapped_ari_matrix[, dataset_index, drop = FALSE], 2, mean, na.rm = TRUE)
  
  # Calculate the standard deviation of ARI for the current dataset across all bootstraps
  sd_ari <- apply(bootstrapped_ari_matrix[, dataset_index, drop = FALSE], 2, sd, na.rm = TRUE)
  
  # Append results to the data frame
  results_df_HC_complete <- rbind(results_df_HC_complete, data.frame(Dataset = dataset_index,
                                                                     Mean_Bootstrapped_ARI = mean_ari,
                                                                     SD_Bootstrapped_ARI = sd_ari,
                                                                     Bootstrapped_ARIs = I(list(bootstrapped_ari_matrix[, dataset_index])),
                                                                     p_values = I(list(p_value_ari_matrix[, dataset_index]))
  ))
}

# Modify dataset labels in result_df
dataset_labels <- c("10 MAGs 2x", "20 MAGs 5x", "50 MAGs 10x", "100 MAGs 20x", "Sex", "Screening center", "Antibiotics")
results_df_HC_complete$Dataset <- factor(results_df_HC_complete$Dataset, 
                                         levels = 1:length(dataset_labels), 
                                         labels = dataset_labels)

# Count the number of significant p-values
results_df_HC_complete$Count_P_Values_Below_0_01 <- sapply(results_df_HC_complete$p_values, function(p_vec) {
  sum(p_vec < 0.01, na.rm = TRUE)
})

# Save results for later use
# saveRDS(results_df_HC_complete,
#         "/PATH/scripts/profile_based_design_processing_results/data/results_df_HC_complete_bootstrap_50_wo_pcoa.rds")

# Read results
# results_df_HC_complete <- readRDS("/PATH/scripts/profile_based_design_processing_results/data/results_df_HC_complete_bootstrap_50_wo_pcoa.rds")

# -----------------------------------------------------------------------------
# K means clustering directly on the abundance profiles (no dimension reduction)
# -----------------------------------------------------------------------------

# Datasets to investigate
abundance_data_list <- list(df.long_10_2x,
                            df.long_20_5x,
                            df.long_50_10x,
                            df.long_100_20x,
                            df.long_crcbiome_kjonn,
                            df.long_crcbiome_senter,
                            df.long_crcbiome_antibiotics)

n_bootstraps <- 50

# Initialize an empty data frame to store results
results_df_Kmeans <- data.frame(Dataset = integer(),
                                Bootstrapped_ARIs = I(list()), 
                                p_values = I(list()),
                                Mean_Bootstrapped_ARI = numeric(),
                                SD_Bootstrapped_ARI = numeric())

# Initialize a matrix to store the bootstrapped ARI values for each dataset
bootstrapped_ari_matrix <- matrix(NA, nrow = n_bootstraps, ncol = length(abundance_data_list))

# Initialize a matrix to store the bootstrapped p values for each dataset
p_value_ari_matrix <- matrix(NA, nrow = n_bootstraps, ncol = length(abundance_data_list))


set.seed(123)
for (dataset_index in 1:length(abundance_data_list)) {
  
  wide_data <- abundance_data_list[[dataset_index]] %>% 
    pivot_wider(names_from = MAG, 
                values_from = rel_abundance, 
                id_cols = c(sample_id, group))
  
  for (i_boot in 1:n_bootstraps) {
    
    df <- wide_data %>% 
      slice_sample(n = nrow(.), replace = TRUE)
    
    true_clusters <- as.integer(df$group)
    
    # Perform k means clustering
    kmeans.obj <- df %>% 
      select(-c(sample_id, group)) %>% 
      kmeans(centers = 2, 
             iter.max = 10, 
             nstart = 10)
    
    predicted_clusters <- kmeans.obj$cluster
    
    # Compute ARI value
    ari_val <- as.numeric(adjRand_test(true_clusters, predicted_clusters, perm = 1000)[1])
    
    p_val <- as.numeric(adjRand_test(true_clusters, predicted_clusters, perm = 1000)[2])
    
    # Store the ARI value in the matrix
    bootstrapped_ari_matrix[i_boot, dataset_index] <- ari_val
    
    # Store the p value in the matrix
    p_value_ari_matrix[i_boot, dataset_index] <- p_val
  }
  
  # Calculate the mean ARI for the current dataset across all bootstraps
  mean_ari <- apply(bootstrapped_ari_matrix[, dataset_index, drop = FALSE], 2, mean, na.rm = TRUE)
  
  # Calculate the standard deviation of ARI for the current dataset across all bootstraps
  sd_ari <- apply(bootstrapped_ari_matrix[, dataset_index, drop = FALSE], 2, sd, na.rm = TRUE)
  
  # Append results to the data frame
  results_df_Kmeans <- rbind(results_df_Kmeans, data.frame(Dataset = dataset_index,
                                                           Mean_Bootstrapped_ARI = mean_ari,
                                                           SD_Bootstrapped_ARI = sd_ari,
                                                           Bootstrapped_ARIs = I(list(bootstrapped_ari_matrix[, dataset_index])),
                                                           p_values = I(list(p_value_ari_matrix[, dataset_index]))
  ))
}

# Modify dataset labels in result_df
dataset_labels <- c("10 MAGs 2x", "20 MAGs 5x", "50 MAGs 10x", "100 MAGs 20x", "Sex", "Screening center", "Antibiotics")
results_df_Kmeans$Dataset <- factor(results_df_Kmeans$Dataset, 
                                    levels = 1:length(dataset_labels), 
                                    labels = dataset_labels)

# Count the number of significant p-values
results_df_Kmeans$Count_P_Values_Below_0_01 <- sapply(results_df_Kmeans$p_values, function(p_vec) {
  sum(p_vec < 0.01, na.rm = TRUE)
})

# Save results for later use
# saveRDS(results_df_Kmeans,
#         "/PATH/scripts/profile_based_design_processing_results/data/results_df_kmeans_bootstrap_50_wo_pcoa.rds")

# Read results
# results_df_Kmeans <- readRDS("/PATH/scripts/profile_based_design_processing_results/data/results_df_kmeans_bootstrap_50_wo_pcoa.rds")

# -----------------------------------------------------------------------------
# UMAP plot grid with distance matrix calculated in Python and calculations
# -----------------------------------------------------------------------------

# Exporting data for UMAP to be used in Python (this only needs to be done once)

df_list <- list(df.long_10_2x, df.long_20_5x, df.long_50_10x, df.long_100_20x)

suffixes <- list("10_2x", "20_5x", "50_10x", "100_20x")

directory_path <- "/PATH/umap_input_data/"

export_data_for_umap(df_list, suffixes, directory_path)

# Defining original clusters for simulation data for significance testing (this is the same for all simulated data)
wide_data <- df.long_10_2x %>% 
  pivot_wider(names_from = MAG, 
              values_from = rel_abundance, 
              id_cols = c(sample_id, group))

original_clusters <- as.integer(wide_data$group)

# Defining original clusters for sex data for significance testing
wide_data_sex <- df.long_crcbiome_kjonn %>% 
  pivot_wider(names_from = MAG, 
              values_from = rel_abundance, 
              id_cols = c(sample_id, group))

original_clusters_sex <- as.integer(wide_data_sex$group)

# Defining original clusters for screening center data for significance testing
wide_data_screening_center <- df.long_crcbiome_senter %>% 
  pivot_wider(names_from = MAG, 
              values_from = rel_abundance, 
              id_cols = c(sample_id, group))

original_clusters_screening_center <- as.integer(wide_data_screening_center$group)

# (Run python script)

# Defining file paths to UMAP dimensions
files <- list(
  "/PATH/umap_embedding_data/default/data_10_2x.csv",
  "/PATH/umap_embedding_data/default/data_20_5x.csv",
  "/PATH/umap_embedding_data/default/data_50_10x.csv",
  "/PATH/umap_embedding_data/default/data_100_20x.csv",
  "/PATH/umap_embedding_data/default/data_sex.csv",
  "/PATH/umap_embedding_data/default/data_screening_center.csv")

# Creating grid with UMAP plots
umap_plots <- plotting_umap(files[1:4], df.long_10_2x, n_col = 2, n_row = 2, grid_plot_title = "Default parameters") # Remember to change title
umap_grid <- umap_plots$grid

# Calculating ARI values and p-values from clustering with UMAP embeddings 

# 10 MAGs 2x
umap_hc_10_2x <- HC_umap(embedding_file = files[[1]], groups = original_clusters, hc_method = "average") # Remember to change hc_method
umap_hc_10_2x$ari
umap_hc_10_2x$p_val

umap_kmeans_10_2x <- kmeans_umap(embedding_file = files[[1]], groups = original_clusters)
umap_kmeans_10_2x$ari
umap_kmeans_10_2x$p_val

# 20 MAGs 5x
umap_hc_20_5x <- HC_umap(embedding_file = files[[2]], groups = original_clusters, hc_method = "average") # Remember to change hc_method
umap_hc_20_5x$ari
umap_hc_20_5x$p_val

umap_kmeans_20_5x <- kmeans_umap(embedding_file = files[[2]], groups = original_clusters)
umap_kmeans_20_5x$ari
umap_kmeans_20_5x$p_val

# 50 MAGs 10x
umap_hc_50_10x <- HC_umap(embedding_file = files[[3]], groups = original_clusters, hc_method = "average") # Remember to change hc_method
umap_hc_50_10x$ari
umap_hc_50_10x$p_val

umap_kmeans_50_10x <- kmeans_umap(embedding_file = files[[3]], groups = original_clusters)
umap_kmeans_50_10x$ari
umap_kmeans_50_10x$p_val

# 100 MAGs 20x
umap_hc_100_20x <- HC_umap(embedding_file = files[[4]], groups = original_clusters, hc_method = "average") # Remember to change hc_method
umap_hc_100_20x$ari
umap_hc_100_20x$p_val

umap_kmeans_100_20x <- kmeans_umap(embedding_file = files[[4]], groups = original_clusters)
umap_kmeans_100_20x$ari
umap_kmeans_100_20x$p_val

# Variable = sex
umap_hc_sex <- HC_umap(embedding_file = files[[5]], groups = original_clusters_sex, hc_method = "complete") # Remember to change hc_method
umap_hc_sex$ari
umap_hc_sex$p_val

umap_kmeans_sex <- kmeans_umap(embedding_file = files[[5]], groups = original_clusters_sex)
umap_kmeans_sex$ari
umap_kmeans_sex$p_val

# Variable = screening center
umap_hc_screening_center <- HC_umap(embedding_file = files[[6]], groups = original_clusters_screening_center, hc_method = "complete") # Remember to change hc_method
umap_hc_screening_center$ari
umap_hc_screening_center$p_val

umap_kmeans_screening_center <- kmeans_umap(embedding_file = files[[6]], groups = original_clusters_screening_center)
umap_kmeans_screening_center$ari
umap_kmeans_screening_center$p_val

# -----------------------------------------------------------------------------
# UMAP ARI plots
# -----------------------------------------------------------------------------

# Read UMAP results from other sources
n_components_HC_average <- readRDS("/PATH/umap_plotting_data/n_components_HC_average.rds") %>% 
  mutate(Clustering_method = "HC with average linkage")

n_components_HC_complete <- readRDS("/PATH/umap_plotting_data/n_components_HC_complete.rds") %>% 
  mutate(Clustering_method = "HC with complete linkage")

n_components_kmeans <- readRDS("/PATH/umap_plotting_data/n_components_kmeans.rds") %>% 
  mutate(Clustering_method = "K-means clustering")

# Combine all data
n_components_df <- bind_rows(n_components_HC_average, n_components_HC_complete) %>% 
  bind_rows(n_components_kmeans) %>% 
  mutate(Clustering_method = factor(Clustering_method)) %>% 
  mutate(Dataset = gsub("MetaPhlAn 20 MAGs 5x", "20 MAGs 5x MetaPhlAn", Dataset)) %>% 
  mutate(Dataset = gsub("MetaPhlAn 50 MAGs 10x", "50 MAGs 10x MetaPhlAn", Dataset)) %>% 
  mutate(Dataset = factor(Dataset, levels = c("10 MAGs 2x", "20 MAGs 5x", "20 MAGs 5x MetaPhlAn", "50 MAGs 10x", "50 MAGs 10x MetaPhlAn", "100 MAGs 20x", "Sex", "Screening center")))


# Create bar plot with faceting
n_components_plot <- ggplot(n_components_df, aes(x = n_components, y = ARI, fill = Dataset)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black") +
  facet_wrap(~ Clustering_method, scales = "free", ncol = 1) +
  labs(x = "n_components", y = "ARI") +
  theme_bw() + 
  theme(
    text = element_text(size = 12, face = "bold"),   
    axis.title = element_text(size = 12, face = "bold"),   
    axis.text = element_text(size = 12, face = "bold"),  
    legend.title = element_text(size = 12, face = "bold"),   
    legend.text = element_text(size = 12, face = "bold"),   
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    strip.text.x = element_text(size = 14, face = "bold")
  )

# Read UMAP results from other sources
n_neighbors_HC_average <- readRDS("/PATH/umap_plotting_data/n_neighbors_HC_average.rds") %>% 
  mutate(Clustering_method = "HC with average linkage")

n_neighbors_HC_complete <- readRDS("/PATH/umap_plotting_data/n_neighbors_HC_complete.rds") %>% 
  mutate(Clustering_method = "HC with complete linkage")

n_neighbors_kmeans <- readRDS("/PATH/umap_plotting_data/n_neighbors_kmeans.rds") %>% 
  mutate(Clustering_method = "K-means clustering")

n_neighbors_df <- bind_rows(n_neighbors_HC_average, n_neighbors_HC_complete) %>% 
  bind_rows(n_neighbors_kmeans) %>% 
  mutate(Clustering_method = factor(Clustering_method)) %>% 
  mutate(Dataset = gsub("MetaPhlAn 20 MAGs 5x", "20 MAGs 5x MetaPhlAn", Dataset)) %>% 
  mutate(Dataset = gsub("MetaPhlAn 50 MAGs 10x", "50 MAGs 10x MetaPhlAn", Dataset)) %>% 
  mutate(Dataset = factor(Dataset, levels = c("10 MAGs 2x", "20 MAGs 5x", "20 MAGs 5x MetaPhlAn", "50 MAGs 10x", "50 MAGs 10x MetaPhlAn", "100 MAGs 20x", "Sex", "Screening center")))

# Create bar plot with faceting
n_neighbors_plot <- ggplot(n_neighbors_df, aes(x = n_neighbors, y = ARI, fill = Dataset)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black") +
  facet_wrap(~ Clustering_method, scales = "free", ncol = 1) +
  labs(x = "n_neighbors", y = "ARI") +
  theme_bw() + 
  theme(
    text = element_text(size = 12, face = "bold"),   
    axis.title = element_text(size = 12, face = "bold"),   
    axis.text = element_text(size = 12, face = "bold"),  
    legend.title = element_text(size = 12, face = "bold"),   
    legend.text = element_text(size = 12, face = "bold"),   
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    strip.text.x = element_text(size = 14, face = "bold")
  )

# Read UMAP results from other sources
min_dist_HC_average <- readRDS("/PATH/umap_plotting_data/min_dist_HC_average.rds") %>% 
  mutate(Clustering_method = "HC with average linkage")

min_dist_HC_complete <- readRDS("/PATH/umap_plotting_data/min_dist_HC_complete.rds") %>% 
  mutate(Clustering_method = "HC with complete linkage")

min_dist_kmeans <- readRDS("/PATH/umap_plotting_data/min_dist_kmeans.rds") %>% 
  mutate(Clustering_method = "K-means clustering")

min_dist_df <- bind_rows(min_dist_HC_average, min_dist_HC_complete) %>% 
  bind_rows(min_dist_kmeans) %>% 
  mutate(Clustering_method = factor(Clustering_method)) %>% 
  mutate(Dataset = gsub("MetaPhlAn 20 MAGs 5x", "20 MAGs 5x MetaPhlAn", Dataset)) %>% 
  mutate(Dataset = gsub("MetaPhlAn 50 MAGs 10x", "50 MAGs 10x MetaPhlAn", Dataset)) %>% 
  mutate(Dataset = factor(Dataset, levels = c("10 MAGs 2x", "20 MAGs 5x", "20 MAGs 5x MetaPhlAn", "50 MAGs 10x", "50 MAGs 10x MetaPhlAn", "100 MAGs 20x", "Sex", "Screening center")))


# Create bar plot with faceting
min_dist_plot <- ggplot(min_dist_df, aes(x = min_dist, y = ARI, fill = Dataset)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black") +
  facet_wrap(~ Clustering_method, scales = "free", ncol = 1) +
  labs(x = "min_dist", y = "ARI") +
  theme_bw() + 
  theme(
    text = element_text(size = 12, face = "bold"),   
    axis.title = element_text(size = 12, face = "bold"),   
    axis.text = element_text(size = 12, face = "bold"),  
    legend.title = element_text(size = 12, face = "bold"),   
    legend.text = element_text(size = 12, face = "bold"),   
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    strip.text.x = element_text(size = 14, face = "bold")
  )
