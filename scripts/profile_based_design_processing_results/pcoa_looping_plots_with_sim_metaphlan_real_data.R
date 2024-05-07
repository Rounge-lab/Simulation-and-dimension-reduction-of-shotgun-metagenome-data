# Imports
library(tidyverse)

# ------------------------------------------------------------------------------
# PCoA looping results with hierarchical clustering with average linkage
# (All results in one plot)
# ------------------------------------------------------------------------------

# Reading different results from hierarchical clustering with average linkage on PCoA dimensions
results_df_HC_avg_PCoA <- readRDS("/PATH/scripts/profile_based_design_processing_results/data/results_df_pcoa_loop_HC_avg_bootstrap_50.rds")

results_df_HC_avg_PCoA_metaphlan <- readRDS("/PATH/scripts/profile_based_design_processing_results/data/results_df_pcoa_loop_HC_avg_bootstrap_50_metaphlan.rds")

# Combining the data
results_df_HC_avg_PCoA <- bind_rows(results_df_HC_avg_PCoA, results_df_HC_avg_PCoA_metaphlan) %>% 
  mutate(Dataset = factor(Dataset, levels = c("10 MAGs 2x", "20 MAGs 5x", "20 MAGs 5x MetaPhlAn", "50 MAGs 10x", "50 MAGs 10x MetaPhlAn", "100 MAGs 20x", "Sex", "Screening center")))


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
    y = "Mean ARI",
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

# ------------------------------------------------------------------------------
# PCoA looping results with hierarchical clustering with complete linkage
# (All results in one plot)
# ------------------------------------------------------------------------------

# Reading different results from hierarchical clustering with complete linkage on PCoA dimensions
results_df_HC_complete_PCoA <- readRDS("/PATH/scripts/profile_based_design_processing_results/data/results_df_pcoa_loop_HC_complete_bootstrap_50.rds")

results_df_HC_complete_PCoA_metaphlan <- readRDS("/PATH/scripts/profile_based_design_processing_results/data/results_df_pcoa_loop_HC_complete_bootstrap_50_metaphlan.rds")

# Combining the data
results_df_HC_complete_PCoA <- bind_rows(results_df_HC_complete_PCoA, results_df_HC_complete_PCoA_metaphlan) %>% 
  mutate(Dataset = factor(Dataset, levels = c("10 MAGs 2x", "20 MAGs 5x", "20 MAGs 5x MetaPhlAn", "50 MAGs 10x", "50 MAGs 10x MetaPhlAn", "100 MAGs 20x", "Sex", "Screening center")))

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
    y = "Mean ARI",
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

# ------------------------------------------------------------------------------
# PCoA looping results with K-means clustering (All results in one plot)
# ------------------------------------------------------------------------------



# Reading different results from K-means clustering on PCoA dimensions
results_df_Kmeans_PCoA <- readRDS("/PATH/scripts/profile_based_design_processing_results/data/results_df_pcoa_loop_kmeans_bootstrap_50.rds")

results_df_Kmeans_PCoA_metaphlan <- readRDS("/PATH/scripts/profile_based_design_processing_results/data/results_df_pcoa_loop_kmeans_bootstrap_50_metaphlan.rds")

# Combining the data
results_df_Kmeans_PCoA <- bind_rows(results_df_Kmeans_PCoA, results_df_Kmeans_PCoA_metaphlan) %>% 
  mutate(Dataset = factor(Dataset, levels = c("10 MAGs 2x", "20 MAGs 5x", "20 MAGs 5x MetaPhlAn", "50 MAGs 10x", "50 MAGs 10x MetaPhlAn", "100 MAGs 20x", "Sex", "Screening center")))

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
    y = "Mean ARI",
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
