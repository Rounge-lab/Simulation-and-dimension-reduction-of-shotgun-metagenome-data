# Imports
library(tidyverse)
library(vegan)
library(glue)
library(mclust)
library(ggpubr)
library(grid)

# Sourcing functions for processing of data
source("/PATH/scripts/functions/all_functions.R")

# -----------------------------------------------------------------------------
# Collecting and preparing MetaPhlAn data from simulations
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

# PERMANOVA analysis
permanova_20_5x_metaphlan <- adonis2(formula = metaphlan_20_5x.wide[, -c(1:2)] ~ group,
                                     data = metaphlan_20_5x.wide,
                                     method = "bray",
                                     permutations = 1000)

permanova_50_10x_metaphlan <- adonis2(formula = metaphlan_50_10x.wide[, -c(1:2)] ~ group,
                                      data = metaphlan_50_10x.wide,
                                      method = "bray",
                                      permutations = 1000)

# -----------------------------------------------------------------------------
# PCoA plots
# -----------------------------------------------------------------------------

pcoa_20_5x_metaphlan <- metaphlan_pcoa_plot(df.wide = metaphlan_20_5x.wide,
                                            method = "bray",
                                            k = 200)

pcoa_50_10x_metaphlan <- metaphlan_pcoa_plot(df.wide = metaphlan_50_10x.wide,
                                             method = "bray",
                                             k = 200)

pcoa_plot_20_5x <- pcoa_20_5x_metaphlan$pcoa_plot + labs(title = "20 MAGs 5x MetaPhlAn")
pcoa_plot_50_10x <- pcoa_50_10x_metaphlan$pcoa_plot + labs(title = "50 MAGs 10x MetaPhlAn")

scree_plot_20_5x <- pcoa_20_5x_metaphlan$scree_plot + labs(title = "") + xlim(0, 20 + 1) + scale_y_continuous(limits = c(0, 15)) + coord_fixed()
scree_plot_50_10x <- pcoa_50_10x_metaphlan$scree_plot + labs(title = "") + xlim(0, 20 + 1) + scale_y_continuous(limits = c(0, 15)) + coord_fixed()

plots <- list(pcoa_plot_20_5x + scale_y_continuous(limits = c(-0.5, 0.5)) + scale_x_continuous(limits = c(-0.4, 0.4)) + coord_fixed(), 
              pcoa_plot_50_10x + scale_y_continuous(limits = c(-0.5, 0.5)) + scale_x_continuous(limits = c(-0.4, 0.4)) + coord_fixed(), 
              scree_plot_20_5x,
              scree_plot_50_10x)

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
  ncol = 2, 
  nrow = 2,
  common.legend = TRUE,
  legend = "bottom",
  labels = "AUTO"
)

# -----------------------------------------------------------------------------
# Finding the number of PCoA components to work with in the clustering for each dataset 
# (HC with average linkage) with bootstrapping
# -----------------------------------------------------------------------------

# Datasets to investigate
metaphlan_data_list <- list(metaphlan_20_5x.wide,
                            metaphlan_50_10x.wide)

n_bootstraps <- 50
dimension_range <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 50)

# Initialize an empty data frame to store results
results_df_HC_avg_PCoA_metaphlan <- data.frame(Dataset = integer(),
                                               Dimensions = integer(),
                                               Mean_Variance_Explained = integer(),
                                               Bootstrapped_ARIs = I(list()),
                                               p_values = I(list()),
                                               Mean_Bootstrapped_ARI = numeric(),
                                               SD_Bootstrapped_ARI = numeric())

set.seed(123)
for (dataset_index in 1:length(metaphlan_data_list)) {
  
  wide_data <- metaphlan_data_list[[dataset_index]]
  
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
  results_df_HC_avg_PCoA_metaphlan <- rbind(results_df_HC_avg_PCoA_metaphlan, data.frame(Dataset = dataset_index,
                                                                                         Dimensions = dimension_range,
                                                                                         Mean_Bootstrapped_ARI = mean_ari_per_dimension,
                                                                                         SD_Bootstrapped_ARI = sd_ari_per_dimension,
                                                                                         Bootstrapped_ARIs = I(split(bootstrapped_ari_matrix, col(bootstrapped_ari_matrix))),
                                                                                         Mean_Variance_Explained = mean_variance_explained,
                                                                                         p_values = I(split(p_value_ari_matrix, col(p_value_ari_matrix)))))
}

# Modify dataset labels in result_df
dataset_labels <- c("20 MAGs 5x MetaPhlAn", "50 MAGs 10x MetaPhlAn")
results_df_HC_avg_PCoA_metaphlan$Dataset <- factor(results_df_HC_avg_PCoA_metaphlan$Dataset, 
                                                   levels = 1:length(dataset_labels), 
                                                   labels = dataset_labels)

# Count number of significant p-values
results_df_HC_avg_PCoA_metaphlan$Count_P_Values_Below_0_01 <- sapply(results_df_HC_avg_PCoA_metaphlan$p_values, function(p_vec) {
  sum(p_vec < 0.01, na.rm = TRUE)
})

# Save results for later use
# saveRDS(results_df_HC_avg_PCoA_metaphlan,
#         "/PATH/scripts/profile_based_design_processing_results/data/results_df_pcoa_loop_HC_avg_bootstrap_50_metaphlan.rds")

# Reading results
# results_df_HC_avg_PCoA_metaphlan <- readRDS("/PATH/scripts/profile_based_design_processing_results/data/results_df_pcoa_loop_HC_avg_bootstrap_50_metaphlan.rds")

# Finding the maximum ARI for each dataset
max_ari_per_dataset <- results_df_HC_avg_PCoA_metaphlan %>%
  group_by(Dataset) %>%
  summarize(Max_ARI = max(Mean_Bootstrapped_ARI))

# Join this back with the original data to get a flag for the max ARI
results_df_HC_avg_PCoA_metaphlan <- results_df_HC_avg_PCoA_metaphlan %>%
  left_join(max_ari_per_dataset, by = "Dataset") %>%
  mutate(Is_Max_ARI = Mean_Bootstrapped_ARI == Max_ARI)

# Finding which dimension has the maximum ARI for each dataset
max_ari_df <- results_df_HC_avg_PCoA_metaphlan %>% 
  filter(Is_Max_ARI == TRUE) %>% 
  select(Dataset, Dimensions, Mean_Bootstrapped_ARI)

# Plotting ARI score across the dimensions. Maximum ARI is marked in the plot
pcoa_looping_plot <- ggplot(results_df_HC_avg_PCoA_metaphlan, aes(x = Dimensions, y = Mean_Bootstrapped_ARI, group = Dataset, color = Dataset)) +
  geom_line(linewidth = 1) +
  geom_point(data = filter(results_df_HC_avg_PCoA_metaphlan, Is_Max_ARI), aes(x = Dimensions, y = Mean_Bootstrapped_ARI, color = Dataset), size = 3, shape = 19) +
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
# Finding the number of PCoA components to work with in the clustering for each dataset 
# (HC with complete linkage) with bootstrapping
# -----------------------------------------------------------------------------

# Datasets to investigate
metaphlan_data_list <- list(metaphlan_20_5x.wide,
                            metaphlan_50_10x.wide)

n_bootstraps <- 50
dimension_range <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 50)

# Initialize an empty data frame to store results
results_df_HC_complete_PCoA_metaphlan <- data.frame(Dataset = integer(),
                                                    Dimensions = integer(),
                                                    Mean_Variance_Explained = integer(),
                                                    Bootstrapped_ARIs = I(list()),
                                                    p_values = I(list()),
                                                    Mean_Bootstrapped_ARI = numeric(),
                                                    SD_Bootstrapped_ARI = numeric())

set.seed(123)
for (dataset_index in 1:length(metaphlan_data_list)) {
  
  wide_data <- metaphlan_data_list[[dataset_index]]
  
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
  
  # Now, compute the mean and SD for each dimension across all bootstraps
  mean_ari_per_dimension <- apply(bootstrapped_ari_matrix, 2, mean, na.rm = TRUE)
  sd_ari_per_dimension <- apply(bootstrapped_ari_matrix, 2, sd, na.rm = TRUE)
  mean_variance_explained <- sapply(variance_explained_sums, mean, na.rm = TRUE)
  
  # Append results to the data frame
  results_df_HC_complete_PCoA_metaphlan <- rbind(results_df_HC_complete_PCoA_metaphlan, data.frame(Dataset = dataset_index,
                                                                                                   Dimensions = dimension_range,
                                                                                                   Mean_Bootstrapped_ARI = mean_ari_per_dimension,
                                                                                                   SD_Bootstrapped_ARI = sd_ari_per_dimension,
                                                                                                   Bootstrapped_ARIs = I(split(bootstrapped_ari_matrix, col(bootstrapped_ari_matrix))),
                                                                                                   Mean_Variance_Explained = mean_variance_explained,
                                                                                                   p_values = I(split(p_value_ari_matrix, col(p_value_ari_matrix)))))
}

# Modify dataset labels in result_df
dataset_labels <- c("20 MAGs 5x MetaPhlAn", "50 MAGs 10x MetaPhlAn")
results_df_HC_complete_PCoA_metaphlan$Dataset <- factor(results_df_HC_complete_PCoA_metaphlan$Dataset, 
                                                        levels = 1:length(dataset_labels), 
                                                        labels = dataset_labels)

# Count number of significant p-values
results_df_HC_complete_PCoA_metaphlan$Count_P_Values_Below_0_01 <- sapply(results_df_HC_complete_PCoA_metaphlan$p_values, function(p_vec) {
  sum(p_vec < 0.01, na.rm = TRUE)
})

# Save results for later use
# saveRDS(results_df_HC_complete_PCoA_metaphlan,
#         "/PATH/scripts/profile_based_design_processing_results/data/results_df_pcoa_loop_HC_complete_bootstrap_50_metaphlan.rds")

# Read results
# results_df_HC_complete_PCoA_metaphlan <- readRDS("/PATH/scripts/profile_based_design_processing_results/data/results_df_pcoa_loop_HC_complete_bootstrap_50_metaphlan.rds")

# Finding the maximum ARI for each dataset
max_ari_per_dataset <- results_df_HC_complete_PCoA_metaphlan %>%
  group_by(Dataset) %>%
  summarize(Max_ARI = max(Mean_Bootstrapped_ARI))

# Join this back with the original data to get a flag for the max ARI
results_df_HC_complete_PCoA_metaphlan <- results_df_HC_complete_PCoA_metaphlan %>%
  left_join(max_ari_per_dataset, by = "Dataset") %>%
  mutate(Is_Max_ARI = Mean_Bootstrapped_ARI == Max_ARI)

# Finding which dimension has the maximum ARI for each dataset
max_ari_df <- results_df_HC_complete_PCoA_metaphlan %>% 
  filter(Is_Max_ARI == TRUE) %>% 
  select(Dataset, Dimensions, Mean_Bootstrapped_ARI)

# Plotting ARI score across the dimensions. Maximum ARI is marked in the plot
pcoa_looping_plot <- ggplot(results_df_HC_complete_PCoA_metaphlan, aes(x = Dimensions, y = Mean_Bootstrapped_ARI, group = Dataset, color = Dataset)) +
  geom_line(linewidth = 1) +
  geom_point(data = filter(results_df_HC_complete_PCoA_metaphlan, Is_Max_ARI), aes(x = Dimensions, y = Mean_Bootstrapped_ARI, color = Dataset), size = 3, shape = 19) +
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
# Finding the number of PCoA components to work with in the clustering for each dataset 
# K means clustering with bootstrapping
# ------------------------------------------------------------------------------
# Datasets to investigate
metaphlan_data_list <- list(metaphlan_20_5x.wide,
                            metaphlan_50_10x.wide)

n_bootstraps <- 50
dimension_range <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 50)

# Initialize an empty data frame to store results
results_df_Kmeans_PCoA_metaphlan <- data.frame(Dataset = integer(),
                                               Dimensions = integer(),
                                               Mean_Variance_Explained = integer(),
                                               Bootstrapped_ARIs = I(list()),
                                               p_values = I(list()),
                                               Mean_Bootstrapped_ARI = numeric(),
                                               SD_Bootstrapped_ARI = numeric())

set.seed(123)
for (dataset_index in 1:length(metaphlan_data_list)) {
  
  wide_data <- metaphlan_data_list[[dataset_index]]
  
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
  results_df_Kmeans_PCoA_metaphlan <- rbind(results_df_Kmeans_PCoA_metaphlan, data.frame(Dataset = dataset_index,
                                                                                         Dimensions = dimension_range,
                                                                                         Mean_Bootstrapped_ARI = mean_ari_per_dimension,
                                                                                         SD_Bootstrapped_ARI = sd_ari_per_dimension,
                                                                                         Bootstrapped_ARIs = I(split(bootstrapped_ari_matrix, col(bootstrapped_ari_matrix))),
                                                                                         Mean_Variance_Explained = mean_variance_explained,
                                                                                         p_values = I(split(p_value_ari_matrix, col(p_value_ari_matrix)))))
}

# Modify dataset labels in result_df
dataset_labels <- c("20 MAGs 5x MetaPhlAn", "50 MAGs 10x MetaPhlAn")
results_df_Kmeans_PCoA_metaphlan$Dataset <- factor(results_df_Kmeans_PCoA_metaphlan$Dataset, 
                                                   levels = 1:length(dataset_labels), 
                                                   labels = dataset_labels)

# Count the number of significant p-values
results_df_Kmeans_PCoA_metaphlan$Count_P_Values_Below_0_01 <- sapply(results_df_Kmeans_PCoA_metaphlan$p_values, function(p_vec) {
  sum(p_vec < 0.01, na.rm = TRUE)
})

# Save results for later use
# saveRDS(results_df_Kmeans_PCoA_metaphlan,
#         "/PATH/scripts/profile_based_design_processing_results/data/results_df_pcoa_loop_kmeans_bootstrap_50_metaphlan.rds")

# Read results
# results_df_Kmeans_PCoA_metaphlan <- readRDS("/PATH/scripts/profile_based_design_processing_results/data/results_df_pcoa_loop_kmeans_bootstrap_50_metaphlan.rds")

# Finding the maximum ARI for each dataset
max_ari_per_dataset <- results_df_Kmeans_PCoA_metaphlan %>%
  group_by(Dataset) %>%
  summarize(Max_ARI = max(Mean_Bootstrapped_ARI))

# Join this back with the original data to get a flag for the max ARI
results_df_Kmeans_PCoA_metaphlan <- results_df_Kmeans_PCoA_metaphlan %>%
  left_join(max_ari_per_dataset, by = "Dataset") %>%
  mutate(Is_Max_ARI = Mean_Bootstrapped_ARI == Max_ARI)

# Finding which dimension has the maximum ARI for each dataset
max_ari_df <- results_df_Kmeans_PCoA_metaphlan %>% 
  filter(Is_Max_ARI == TRUE) %>% 
  select(Dataset, Dimensions, Mean_Bootstrapped_ARI)

# Plotting ARI score across the dimensions. Maximum ARI is marked in the plot
pcoa_looping_plot <- ggplot(results_df_Kmeans_PCoA_metaphlan, aes(x = Dimensions, y = Mean_Bootstrapped_ARI, group = Dataset, color = Dataset)) +
  geom_line(linewidth = 1) +
  geom_point(data = filter(results_df_Kmeans_PCoA_metaphlan, Is_Max_ARI), aes(x = Dimensions, y = Mean_Bootstrapped_ARI, color = Dataset), size = 3, shape = 19) +
  theme_bw() +
  labs(
    x = "Number of PCoA Dimensions",
    y = "Mean Adjusted Rand Index",
    color = "Dataset"
  ) +
  guides(color = guide_legend(override.aes = list(size=3))) +
  theme(
    text = element_text(size = 12, face = "bold"),  # Base size for all text in the plot; you can adjust this as needed
    axis.title = element_text(size = 12),  # Adjust axis title size
    axis.text = element_text(size = 12),  # Adjust axis text size
    legend.title = element_text(size = 12),  # Adjust legend title size
    legend.text = element_text(size = 12),  # Adjust legend text size
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  # Adjust plot title size and make it bold
  )

# -----------------------------------------------------------------------------
# Hierarchical clustering directly on the MetaPhlAn data (no dimension reduction)
# HC average linkage
# -----------------------------------------------------------------------------

# Datasets to investigate
metaphlan_data_list <- list(metaphlan_20_5x.wide,
                            metaphlan_50_10x.wide)

n_bootstraps <- 50

# Initialize an empty data frame to store results
results_df_HC_avg_metaphlan <- data.frame(Dataset = integer(),
                                          Bootstrapped_ARIs = I(list()),
                                          p_values = I(list()),
                                          Mean_Bootstrapped_ARI = numeric(),
                                          SD_Bootstrapped_ARI = numeric())

# Initialize a matrix to store the bootstrapped ARI values for each dataset
bootstrapped_ari_matrix <- matrix(NA, nrow = n_bootstraps, ncol = length(metaphlan_data_list))

# Initialize a matrix to store the bootstrapped p values for each dataset
p_value_ari_matrix <- matrix(NA, nrow = n_bootstraps, ncol = length(metaphlan_data_list))

set.seed(123)
for (dataset_index in 1:length(metaphlan_data_list)) {
  
  wide_data <- metaphlan_data_list[[dataset_index]]
  
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
  results_df_HC_avg_metaphlan <- rbind(results_df_HC_avg_metaphlan, data.frame(Dataset = dataset_index,
                                                                               Mean_Bootstrapped_ARI = mean_ari,
                                                                               SD_Bootstrapped_ARI = sd_ari,
                                                                               Bootstrapped_ARIs = I(list(bootstrapped_ari_matrix[, dataset_index])),
                                                                               p_values = I(list(p_value_ari_matrix[, dataset_index]))
  ))
}

# Modify dataset labels in result_df
dataset_labels <- c("20 MAGs 5x", "50 MAGs 10x")
results_df_HC_avg_metaphlan$Dataset <- factor(results_df_HC_avg_metaphlan$Dataset, 
                                              levels = 1:length(dataset_labels), 
                                              labels = dataset_labels)

# Count the number of significant p-values
results_df_HC_avg_metaphlan$Count_P_Values_Below_0_01 <- sapply(results_df_HC_avg_metaphlan$p_values, function(p_vec) {
  sum(p_vec < 0.01, na.rm = TRUE)
})

# Save results for later use
# saveRDS(results_df_HC_avg_metaphlan,
#         "/PATH/scripts/profile_based_design_processing_results/data/results_df_HC_avg_bootstrap_50_wo_pcoa_metaphlan.rds")

# Read results
# results_df_HC_avg_metaphlan <- readRDS("/PATH/scripts/profile_based_design_processing_results/data/results_df_HC_avg_bootstrap_50_wo_pcoa_metaphlan.rds")

# -----------------------------------------------------------------------------
# Hierarchical clustering directly on the MetaPhlAn data (no dimension reduction)
# HC complete linkage
# -----------------------------------------------------------------------------

# Datasets to investigate
metaphlan_data_list <- list(metaphlan_20_5x.wide,
                            metaphlan_50_10x.wide)

n_bootstraps <- 50

# Initialize an empty data frame to store results
results_df_HC_complete_metaphlan <- data.frame(Dataset = integer(),
                                               Bootstrapped_ARIs = I(list()),
                                               p_values = I(list()),
                                               Mean_Bootstrapped_ARI = numeric(),
                                               SD_Bootstrapped_ARI = numeric())

# Initialize a matrix to store the bootstrapped ARI values for each dataset
bootstrapped_ari_matrix <- matrix(NA, nrow = n_bootstraps, ncol = length(metaphlan_data_list))

# Initialize a matrix to store the bootstrapped p values for each dataset
p_value_ari_matrix <- matrix(NA, nrow = n_bootstraps, ncol = length(metaphlan_data_list))

set.seed(123)
for (dataset_index in 1:length(metaphlan_data_list)) {
  
  wide_data <- metaphlan_data_list[[dataset_index]]
  
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
  results_df_HC_complete_metaphlan <- rbind(results_df_HC_complete_metaphlan, data.frame(Dataset = dataset_index,
                                                                                         Mean_Bootstrapped_ARI = mean_ari,
                                                                                         SD_Bootstrapped_ARI = sd_ari,
                                                                                         Bootstrapped_ARIs = I(list(bootstrapped_ari_matrix[, dataset_index])),
                                                                                         p_values = I(list(p_value_ari_matrix[, dataset_index]))
  ))
}

# Modify dataset labels in result_df
dataset_labels <- c("20 MAGs 5x", "50 MAGs 10x")
results_df_HC_complete_metaphlan$Dataset <- factor(results_df_HC_complete_metaphlan$Dataset, 
                                                   levels = 1:length(dataset_labels), 
                                                   labels = dataset_labels)

# Count the number of significant p-values
results_df_HC_complete_metaphlan$Count_P_Values_Below_0_01 <- sapply(results_df_HC_complete_metaphlan$p_values, function(p_vec) {
  sum(p_vec < 0.01, na.rm = TRUE)
})

# Save results for later use
# saveRDS(results_df_HC_complete_metaphlan,
#         "/PATH/scripts/profile_based_design_processing_results/data/results_df_HC_complete_bootstrap_50_wo_pcoa_metaphlan.rds")

# Read results
# results_df_HC_complete_metaphlan <- readRDS("/PATH/scripts/profile_based_design_processing_results/data/results_df_HC_complete_bootstrap_50_wo_pcoa_metaphlan.rds")

# -----------------------------------------------------------------------------
# K means clustering directly on the abundance profiles (no dimension reduction)
# -----------------------------------------------------------------------------

# Datasets to investigate
metaphlan_data_list <- list(metaphlan_20_5x.wide,
                            metaphlan_50_10x.wide)

n_bootstraps <- 50

# Initialize an empty data frame to store results
results_df_Kmeans_metaphlan <- data.frame(Dataset = integer(),
                                          Bootstrapped_ARIs = I(list()),
                                          p_values = I(list()),
                                          Mean_Bootstrapped_ARI = numeric(),
                                          SD_Bootstrapped_ARI = numeric())

# Initialize a matrix to store the bootstrapped ARI values for each dataset
bootstrapped_ari_matrix <- matrix(NA, nrow = n_bootstraps, ncol = length(metaphlan_data_list))

# Initialize a matrix to store the bootstrapped p values for each dataset
p_value_ari_matrix <- matrix(NA, nrow = n_bootstraps, ncol = length(metaphlan_data_list))


set.seed(123)
for (dataset_index in 1:length(metaphlan_data_list)) {
  
  wide_data <- metaphlan_data_list[[dataset_index]] 
  
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
  results_df_Kmeans_metaphlan <- rbind(results_df_Kmeans_metaphlan, data.frame(Dataset = dataset_index,
                                                                               Mean_Bootstrapped_ARI = mean_ari,
                                                                               SD_Bootstrapped_ARI = sd_ari,
                                                                               Bootstrapped_ARIs = I(list(bootstrapped_ari_matrix[, dataset_index])),
                                                                               p_values = I(list(p_value_ari_matrix[, dataset_index]))
  ))
}

# Modify dataset labels in result_df
dataset_labels <- c("20 MAGs 5x", "50 MAGs 10x")
results_df_Kmeans_metaphlan$Dataset <- factor(results_df_Kmeans_metaphlan$Dataset, 
                                              levels = 1:length(dataset_labels), 
                                              labels = dataset_labels)

# Count the number of significant p-values
results_df_Kmeans_metaphlan$Count_P_Values_Below_0_01 <- sapply(results_df_Kmeans_metaphlan$p_values, function(p_vec) {
  sum(p_vec < 0.01, na.rm = TRUE)
})

# Save results for later use
# saveRDS(results_df_Kmeans_metaphlan,
#         "/PATH/scripts/profile_based_design_processing_results/data/results_df_kmeans_bootstrap_50_wo_pcoa_metaphlan.rds")

# -----------------------------------------------------------------------------
# UMAP
# -----------------------------------------------------------------------------

# Extracting the group variable (which is the same for both MetaPhlAn datasets)
groups <- metaphlan_20_5x.wide$group

# Define the files where the UMAP dimensions can be found
files <- list("/PATH/umap_embedding_data_metaphlan/some_directory/data.csv",
              "/PATH/umap_embedding_data_metaphlan/some_directory/data.csv")

# Plotting UMAP plots
umap_plots_metaphlan <- plotting_umap_metaphlan(files, groups, n_col = 2, n_row = 1, grid_plot_title = "Plot title")
parameters_plot_metaphlan <- umap_plots_metaphlan$grid

# Calculating ARI value and its p-value for each dataset for each clustering method, using the UMAP dimensions
HC_umap_avg_20_5x <- HC_umap(groups, files[[1]], hc_method = "average")
HC_umap_avg_20_5x$ari
HC_umap_avg_20_5x$p_val

HC_umap_complete_20_5x <- HC_umap(groups, files[[1]], hc_method = "complete")
HC_umap_complete_20_5x$ari
HC_umap_complete_20_5x$p_val

kmeans_20_5x <- kmeans_umap(groups, files[[1]])
kmeans_20_5x$ari
kmeans_20_5x$p_val

HC_umap_avg_50_10x <- HC_umap(groups, files[[2]], hc_method = "average")
HC_umap_avg_50_10x$ari
HC_umap_avg_50_10x$p_val

HC_umap_complete_50_10x <- HC_umap(groups, files[[2]], hc_method = "complete")
HC_umap_complete_50_10x$ari
HC_umap_complete_50_10x$p_val

kmeans_50_10x <- kmeans_umap(groups, files[[2]])
kmeans_50_10x$ari
kmeans_50_10x$p_val

