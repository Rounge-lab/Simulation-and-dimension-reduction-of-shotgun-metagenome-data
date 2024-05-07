# process_metaphlan_table_rel_abundance_filtering
# Function to process raw metaphlan tables to contain relative abundance values and only species level
#
# Args: 
#   filepath: String, the absolute or relative filepath to the metaphlan table (e.g., "tables/profiles.tsv").
#             The file should be a tab-delimited file.
#
# Returns:
#   A data frame containing the processed metaphlan table with only species level ('s') and
#   relative abundances, with species as rows and samples as columns.
#
# Example:
#   processed_data <- process_metaphlan_table_rel_abundance_filtering("path/to/metaphlan_output.tsv")

process_metaphlan_table_rel_abundance_filtering <- function(filepath) {
  
  # Filtering for species-level data by selecting rows where clade name contains "s__",
  # and excluding rows containing "t__"
  metaphlan_table_filt <- read.table(filepath, sep = "\t", header = TRUE)%>% 
    filter(str_detect(clade_name, "s__")) %>% 
    filter(!str_detect(clade_name, "t__"))
  
  # Calculating the sum of abuncance values for each sample
  sample_sum <- colSums(metaphlan_table_filt[, -1])
  
  # Convert absolute abundance counts to relative abundances by dividing by the sample sum
  rel_abundance_data <- sweep(metaphlan_table_filt[, -1], 2, sample_sum, FUN = "/") %>% 
    mutate(clade_name = metaphlan_table_filt$clade_name) %>% 
    select(clade_name, everything())
  
  return(rel_abundance_data) # Species as rows and samples as columns
}

# process_metaphlan_table_pcoa
# Prepares a MetaPhlAn table for Principal Coordinates Analysis (PCoA).
# It transposes the table, converts abundance values to numeric, and adds a data source identifier.
#
# Args: 
#   metaphlan_df: Data frame, in the same format as process_metaphlan_table_rel_abundance_filtering produces 
#                 The data frame should have species as rows and samples as columns with relative abundance values (0-1).
#   data_source : String, a discriptor indicating the data's origin or grouping variable for the analysis.
#
# Returns:
#   A data frame ready for PCoA, with samples as rows and species as columns, including a 'sample_id' column
#   and a 'data_source' column for grouping or identifying the data.
#
# Example:
#   metaphlan_pcoa.tbl <- process_metaphlan_table_pcoa(metaphlan_df, "Dataset1")


process_metaphlan_table_pcoa <- function(metaphlan_df, data_source) {
  # Transpose the table to fit cmdscale() requirements with samples as rows and species as columns
  metaphlan_table <- t(metaphlan_df) %>% 
    as.data.frame()
  
  # Assign column names and create a 'sample_id' column while removing row names
  colnames(metaphlan_table) <- metaphlan_table[1, ]
  metaphlan_table <- rownames_to_column(metaphlan_table, "sample_id")
  
  # Remove the first row which contains clade names, as it's not needed for analysis
  metaphlan_table <- metaphlan_table[-1, ]
  
  # Store 'sample_id' for later use and convert all abundance values to numeric type
  sample_ids <- metaphlan_table$sample_id
  metaphlan_table <- data.frame(lapply(metaphlan_table[, -1], as.numeric), sample_id = sample_ids)
  
  # Re-introduce 'sample_id' at the beginning, add 'data_source', and select all columns for output
  metaphlan_table_proc <- metaphlan_table %>%
    mutate(data_source = data_source) %>%
    select(sample_id, data_source, everything())
  
  # Return the prepared table
  return(metaphlan_table_proc)
}

# get_random_samples
# Extracts a specified number of random rows (samples) from a MetaPhlAn data frame prepared for PCoA analysis.
#
# Args:
#   metaphlan_data: A data frame containing MetaPhlAn output data specifically structured
#                   by the process_metaphlan_table_pcoa function, with samples as rows and species
#                   as columns.
#   num_rows: Integer, the number of random rows to extract from 'metaphlan_data'.
#   seed: Integer, the seed value for the random number generator to ensure reproducibility. Default is 123.
#
# Returns:
#   A data frame containing a random subset of rows from the 'metaphlan_data'. The number of rows
#   in the returned data frame equals 'num_rows', provided that 'metaphlan_data' contains a sufficient
#   number of rows.
#
# Example:
#   subset_metaphlan_data <- get_random_samples(metaphlan_data, 10, seed = 42)
#
# Note:
#   The function will stop and throw an error if the 'metaphlan_data' has fewer rows than 'num_rows' requested.
#
get_random_samples <- function(metaphlan_data, num_rows, seed = 123) {
  # Set the seed for reproducibility
  set.seed(seed)
  
  # Ensure the data frame has a sufficient number of rows
  if (nrow(metaphlan_data) < num_rows) {
    stop("The data frame has fewer rows than the number requested.")
  }
  
  # Sample 'num_rows' random row indices
  random_indices <- sample(nrow(metaphlan_data), num_rows)
  
  # Subset the data frame to get the random rows (and all columns)
  random_rows_df <- metaphlan_data[random_indices, ]
  
  # Return the subset data frame
  return(random_rows_df)
}

# joining_data
# Joins two MetaPhlAn data frames, which have been pre-processed for PCoA analysis. It combines
# the data frames row-wise and replaces any resulting NA values with zeros. 
#
# Args:
#   df1: A data frame, structured by the process_metaphlan_table_pcoa function, containing MetaPhlAn
#        output data with samples as rows and species as columns.
#   df2: Another data frame with the same structure as 'df1', to be combined with 'df1'.
#
# Returns:
#   A single data frame where the rows from 'df1' and 'df2' are combined. NA values resulting from
#   the join are replaced with zeros.
#
# Example:
#   combined_metaphlan_data <- joining_data(processed_metaphlan_data1, processed_metaphlan_data2)
#
# Note:
#   The two data frames must have the same species as columns. If they do not, the combined data frame will
#   have NA values where species are missing in either data frame, which will then be replaced with zeros.
#

joining_data <- function(df1, df2) {
  # Combine the two data frames row-wise
  joined_data <- bind_rows(df1, df2)
  
  # Replacing NAs with zero to handle absent data points
  joined_data[is.na(joined_data)] <- 0
  
  # Return the combined data frame with NAs replaced by zeros
  return(joined_data)
}

# perform_pcoa_and_plot
# Conducts a principal coordinates analysis (PCoA) on a given MetaPhlAn data frame preprocessed
# for PCoA and generates two plots: a PCoA plot and a scree plot. 
#
# Args:
#   df: A data frame containing processed MetaPhlAn output data, structured by the
#       process_metaphlan_table_pcoa function, with samples as rows and species as columns.
#   method: String, the method used to calculate distances between samples, default is "bray".
#   k: Integer, the number of principal coordinates to compute, default is 2.
#   remove_columns: Numeric vector, indices of columns to remove before PCoA, default is c(1:2) (sample_id and data source).
#   log_mu: Numeric, a parameter used in the simulation of data.
#   log_sigma: Numeric, a parameter used in the simulation of data.
#   mode: String, the mode of the simulation, default is "differential abundance".
#
# Returns:
#   A list containing two ggplot objects: 
#     - 'pcoa_plot': The PCoA plot showing the samples in the space defined by the first two principal coordinates.
#     - 'scree_plot': The scree plot showing the percentage of variance explained by each principal coordinate.
#
# Example:
#   list_of_plots <- perform_pcoa_and_plot(processed_metaphlan_df, method = "bray", k = 2, 
#                                          remove_columns = c(1,2), log_mu = 1, log_sigma = 2,
#                                          mode = "differential abundance")
#
# Note:
#   The input data frame should have species as columns and samples as rows, with the first two columns
#   typically containing metadata that can be removed before PCoA. The function assumes that the necessary
#   libraries (vegan, dplyr, ggplot2, and glue) are installed and loaded in the R session.
#
perform_pcoa_and_plot <- function(df, 
                                  method = "bray", 
                                  k = 2, 
                                  remove_columns = c(1:2), 
                                  log_mu, 
                                  log_sigma, 
                                  mode = "differential abundance") {
  # Calculate distance matrix
  dist.mat <- vegdist(df[, -remove_columns], method = method)
  
  # Perform PCoA
  pcoa_mat <- cmdscale(dist.mat, k = k, eig = TRUE)
  
  # Extract coordinates and set column names
  positions <- pcoa_mat$points
  colnames(positions) <- c("PCo1", "PCo2")
  
  # Extract explained variation for each axis
  percent_explained <- 100 * pcoa_mat$eig / sum(pcoa_mat$eig)
  pretty_pe <- round(percent_explained, 2)
  
  # Create labels for the plots
  labs <- c(glue("PCo1 ({pretty_pe[1]}%)"),
            glue("PCo2 ({pretty_pe[2]}%)"))
  
  pcoa_plot_title <- glue("PCoA plot with {mode} mode (mu = {log_mu}, sigma = {log_sigma})")
  
  if (any(grepl("Replicates", df$data_source, ignore.case = TRUE))) {
    source <- "Replicates mode (n = 10)"
    colour <- "#00BA38"
    
  } else {
    source <- "Differential abundance mode (n = 10)"
    colour <- "#619CFF"
  }
  
  # Define colors
  source_colors <- setNames(c("#F8766D", colour), 
                            c("CRCbiome (n = 50)", source))
  
  # Create PCoA plot
  pcoa_plot <- ggplot(bind_cols(df, as.data.frame(positions)),
                      aes(x = PCo1, y = PCo2, color = data_source)) +
    scale_color_manual(values = source_colors) +
    geom_point(size = 6) +
    labs(x = labs[1], y = labs[2], title = pcoa_plot_title) +
    theme_bw()+
    theme(
      text = element_text(size = 20, face = "bold"),  
      axis.title = element_text(size = 20), 
      axis.text = element_text(size = 20),  
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 19),
      plot.title = element_text(size = 16, face = "bold")  
    )
  
  scree_plot_title <- glue("Scree plot with {mode} mode (mu = {log_mu}, sigma = {log_sigma})")
  
  # Create scree plot
  scree_plot <- ggplot(tibble(pe = percent_explained),
                       aes(x = 1:length(pe), y = pe)) +
    geom_col() +
    labs(x = "PCoA axis", y = "Percent explained by axis",
         title = scree_plot_title) +
    theme_bw()
  
  # Return both plots in list
  return(list(pcoa_plot = pcoa_plot, scree_plot = scree_plot))
}


# perform_pcoa_and_plot_3_df
# Performs Principal Coordinates Analysis (PCoA) on 3 MetaPhlAn data frames and creates visualizations.
# The function takes three data frames preprocessed for relative abundance and one preprocessed for PCoA.
# It then joins these data frames, performs PCoA, and generates a PCoA plot and a scree plot.
#
# Args:
#   metaphlan_df_1: A data frame preprocessed by process_metaphlan_table_rel_abundance_filtering, representing the first dataset.
#   metaphlan_df_2: A data frame preprocessed by process_metaphlan_table_rel_abundance_filtering, representing the second dataset.
#   metaphlan_df_3: A data frame preprocessed by process_metaphlan_table_rel_abundance_filtering, representing the third dataset.
#   metaphlan_df_crc: A data frame preprocessed by process_metaphlan_table_pcoa, representing the CRCbiome dataset.
#   source_1: String, a descriptor for the first data source.
#   source_2: String, a descriptor for the second data source.
#   source_3: String, a descriptor for the third data source.
#   method: String, the method used to calculate distances for PCoA, default is "bray".
#   k: Integer, the number of principal coordinates to compute, default is 2.
#   remove_columns: Numeric vector, indices of columns to remove before joining and PCoA, default is c(1:2).
#   alpha: Numeric, transparency level for the points in the PCoA plot, default is 0.8.
#
# Returns:
#   A list containing two elements:
#     - 'pcoa_plot': A ggplot object representing the PCoA plot.
#     - 'scree_plot': A ggplot object representing the scree plot with the percentage of variance explained.
#
# Example:
#   plots <- perform_pcoa_and_plot_3_df(metaphlan_df_1, metaphlan_df_2, metaphlan_df_3,
#                                              metaphlan_df_crc, "Source 1", "Source 2", "Source 3",
#                                              method = "bray", k = 2, remove_columns = c(1:2), alpha = 0.8)
#
# Note:
#   The first three data frames should have species as columns and samples as rows with relative abundance values.
#   The last data frame (metaphlan_df_crc) should be processed specifically for PCoA with species as rows and 
#   samples as columns. The function assumes that all necessary libraries are installed and loaded in the R session.
#
perform_pcoa_and_plot_3_df <- function(metaphlan_df_1, 
                                       metaphlan_df_2, 
                                       metaphlan_df_3, 
                                       metaphlan_df_crc, 
                                       source_1, 
                                       source_2, 
                                       source_3, 
                                       method = "bray", 
                                       k = 2, 
                                       remove_columns = c(1:2),
                                       alpha = 0.8) {
  
  # Processing metaphlan tables for PCoA
  metaphlan_df_1 <- process_metaphlan_table_pcoa(metaphlan_df_1, 
                                                 data_source = source_1)
  metaphlan_df_2 <- process_metaphlan_table_pcoa(metaphlan_df_2, 
                                                 data_source = source_2)
  metaphlan_df_3 <- process_metaphlan_table_pcoa(metaphlan_df_3, 
                                                 data_source = source_3)
  
  # Joining data sets
  joined_data_1 <- joining_data(metaphlan_df_crc, metaphlan_df_1) # Joining df1 and CRCbiome df
  joined_data_1_2 <- joining_data(joined_data_1, metaphlan_df_2) # Joining previous and df2
  joined_data <- joining_data(joined_data_1_2, metaphlan_df_3) # Joining previous and df3
  
  # Calculate distance matrix
  dist.mat <- vegdist(joined_data[, -remove_columns], method = method)
  
  # Perform PCoA
  pcoa_mat <- cmdscale(dist.mat, k = k, eig = TRUE)
  
  # Extract coordinates and set column names
  positions <- pcoa_mat$points
  colnames(positions) <- c("PCo1", "PCo2")
  
  # Extract explained variation for each axis
  percent_explained <- 100 * pcoa_mat$eig / sum(pcoa_mat$eig)
  pretty_pe <- round(percent_explained, 2)
  
  # Create labels for the plots
  labs <- c(glue("PCo1 ({pretty_pe[1]}%)"),
            glue("PCo2 ({pretty_pe[2]}%)"))
  
  # Define colors for plotting
  source_colors <- setNames(c("#F8766D", "#619CFF", "#00BA38", "#C77CFF"), 
                            c("CRCbiome", source_1, source_2, source_3))
  
  # Create PCoA plot
  pcoa_plot <- ggplot(bind_cols(joined_data, as.data.frame(positions)), 
                      aes(x = PCo1, y = PCo2, color = data_source, shape = data_source)) +
    geom_point(size = 4, alpha = alpha) +
    scale_color_manual(values = source_colors) +
    labs(x = labs[1], y = labs[2], color = "", shape = "") + 
    theme_bw() + 
    theme(
      text = element_text(size = 15, face = "bold"), 
      axis.title = element_text(size = 16), 
      axis.text = element_text(size = 14),  
      legend.title = element_text(size = 14), 
      legend.text = element_text(size = 16),  
      plot.title = element_text(size = 16, face = "bold")
    )
  
  # Create scree plot
  scree_plot <- ggplot(tibble(pe = percent_explained), 
                       aes(x = 1:length(pe), y = pe)) +
    geom_col() +
    labs(x = "PCoA axis", y = "Percent explained by axis") + 
    theme_bw()
  
  # Return both plots in list
  return(list(pcoa_plot = pcoa_plot, scree_plot = scree_plot))
}

# samples_per_species_histograms
# This function calculates the number of samples per species within multiple MetaPhlAn data frames
# and creates a histogram grid with data from each of the data sources. The function also executes a Wilcoxon test. 
#
# Args:
#   metaphlan_df_1: A data frame preprocessed by process_metaphlan_table_rel_abundance_filtering, representing the first dataset.
#   metaphlan_df_2: A data frame preprocessed by process_metaphlan_table_rel_abundance_filtering, representing the second dataset.
#   metaphlan_df_3: A data frame preprocessed by process_metaphlan_table_rel_abundance_filtering, representing the third dataset.
#   metaphlan_df_crc: A data frame preprocessed by process_metaphlan_table_pcoa, representing the CRCbiome dataset.
#   source_1: String, a descriptor for the first data source.
#   source_2: String, a descriptor for the second data source.
#   source_3: String, a descriptor for the third data source.
#
# Returns:
#   A list containing one ggplot object:
#     - 'histogram_grid': Combined histogram plot for all data sources.
#
# Example:
#   example <- samples_per_species(metaphlan_df_1, 
#                                  metaphlan_df_2, 
#                                  metaphlan_df_3, 
#                                  metaphlan_df_crc,
#                                  "Source 1", 
#                                  "Source 2", 
#                                  "Source 3")
#   histogram <- example$histogram_grid
#   wilcoxon <- example$wilcoxon_res
# Note:
#   Each input data frames should contain species abundance data with samples as columns and species as rows.
#   The data frames must be compatible in terms of species columns for proper comparison and visualization.
#   The function assumes that the 'ggplot2' library is installed and loaded in the R session.
#
samples_per_species_histograms <- function(metaphlan_df_1, 
                                           metaphlan_df_2, 
                                           metaphlan_df_3, 
                                           metaphlan_df_crc,
                                           source_1, 
                                           source_2, 
                                           source_3){
  
  # Change sample ids to correspond with the simulated samples
  colnames(metaphlan_df_crc) <- colnames(metaphlan_df_1)
  
  # Finding number of samples per species for each data frame
  metaphlan_df_1$samples_per_species <- rowSums(metaphlan_df_1[,2:ncol(metaphlan_df_1)] > 0)
  metaphlan_df_2$samples_per_species <- rowSums(metaphlan_df_2[,2:ncol(metaphlan_df_2)] > 0)
  metaphlan_df_3$samples_per_species <- rowSums(metaphlan_df_3[,2:ncol(metaphlan_df_3)] > 0)
  metaphlan_df_crc$samples_per_species <- rowSums(metaphlan_df_crc[,2:ncol(metaphlan_df_crc)] > 0)
  
  # Filter data frames so that it does not include 0 number of samples for any species
  metaphlan_df_1 <- metaphlan_df_1 %>% 
    filter(samples_per_species != 0)
  
  metaphlan_df_2 <- metaphlan_df_2 %>% 
    filter(samples_per_species != 0)
  
  metaphlan_df_3 <- metaphlan_df_3 %>% 
    filter(samples_per_species != 0)
  
  metaphlan_df_crc <- metaphlan_df_crc %>% 
    filter(samples_per_species != 0)
  
  # Add a new column to each data frame indicating the source
  metaphlan_df_1$source <- source_1
  metaphlan_df_2$source <- source_2
  metaphlan_df_3$source <- source_3
  metaphlan_df_crc$source <-"CRCbiome"
  
  # Combine the data frames
  combined_df <- bind_rows(metaphlan_df_1, metaphlan_df_2, metaphlan_df_3, metaphlan_df_crc)
  
  # Create colors for combined plots
  source_colors <- setNames(c("#F8766D", "#619CFF", "#00BA38", "#C77CFF"), 
                            c("CRCbiome", source_1, source_2, source_3))
  
  histogram_grid <- combined_df %>% 
    select(clade_name, samples_per_species, source) %>% 
    mutate(source = as.factor(source)) %>% 
    ggplot(aes(x = samples_per_species, fill = source)) + 
    geom_histogram(binwidth = 1) + 
    scale_fill_manual(values = source_colors) +
    facet_grid(. ~ source, scales = "free_y") +
    labs(x = "Number of samples",
         y = "Number of species",
         title = "Samples per species",
         fill = "") +  # Set the legend title
    theme_bw() +
    theme(
      text = element_text(size = 15, face = "bold"),  
      axis.title = element_text(size = 16), 
      axis.text = element_text(size = 14),  
      legend.title = element_text(size = 14),  
      legend.text = element_text(size = 14),  
      plot.title = element_text(size = 16, face = "bold"),  
      strip.text.x = element_text(size = 14)
    )
  
  # Wilcoxon test
  wilcoxon_data <- combined_df %>% 
    select(samples_per_species, source) %>% 
    mutate(source = factor(source))
  
  wilcoxon_res <- pairwise.wilcox.test(wilcoxon_data$samples_per_species, wilcoxon_data$source, exact = FALSE, p.adjust.method = "none")
  
  return(list(histogram_grid = histogram_grid, wilcoxon_res = wilcoxon_res))
}

# species_per_sample_boxplot
# This function aggregates the count of species per sample from multiple MetaPhlAn data frames and 
# creates a box plot visualizing the distribution of species counts across different data sources. 
# It processes each MetaPhlAn data frame to calculate the number of species present in each sample 
# and combines this information into a single data frame for plotting. The function also executes a Wilcoxon test.
#
# Args:
#   metaphlan_df_1: A data frame preprocessed by process_metaphlan_table_rel_abundance_filtering, representing the first dataset.
#   metaphlan_df_2: A data frame preprocessed by process_metaphlan_table_rel_abundance_filtering, representing the second dataset.
#   metaphlan_df_3: A data frame preprocessed by process_metaphlan_table_rel_abundance_filtering, representing the third dataset.
#   metaphlan_df_crc: A data frame preprocessed by process_metaphlan_table_pcoa, representing the CRCbiome dataset.
#   source_1: String, a label for the first data source to be used in the plot.
#   source_2: String, a label for the second data source to be used in the plot.
#   source_3: String, a label for the third data source to be used in the plot.
#
# Returns:
#   A list containing a single ggplot object 'box_plot' which is a box plot of the number of species 
#   per sample across the provided data sources.
#
# Example:
#   example <- species_per_sample(metaphlan_df_1, 
#                                 metaphlan_df_2, 
#                                 metaphlan_df_3, 
#                                 metaphlan_df_crc,
#                                 "Source 1", 
#                                 "Source 2", 
#                                 "Source 3")
#   box_plot <- example$boxplot
#   wilcoxon <- example$wilcoxon_res
#
# Note:
#   The input data frames should have samples as columns and species as rows. 
#
species_per_sample_boxplot <- function(metaphlan_df_1, 
                                       metaphlan_df_2, 
                                       metaphlan_df_3, 
                                       metaphlan_df_crc,
                                       source_1, 
                                       source_2, 
                                       source_3){
  
  species_per_sample_1 <- colSums(metaphlan_df_1[, -c(1, ncol(metaphlan_df_1))] > 0) %>% 
    as.data.frame() %>% 
    rownames_to_column("sample_id") %>% 
    mutate(source = source_1)
  colnames(species_per_sample_1)[2] <- c("num_species")
  
  species_per_sample_2 <- colSums(metaphlan_df_2[, -c(1, ncol(metaphlan_df_2))] > 0) %>% 
    as.data.frame() %>% 
    rownames_to_column("sample_id") %>% 
    mutate(source = source_2)
  colnames(species_per_sample_2)[2] <- c("num_species")
  
  species_per_sample_3 <- colSums(metaphlan_df_3[, -c(1, ncol(metaphlan_df_3))] > 0) %>% 
    as.data.frame() %>% 
    rownames_to_column("sample_id") %>% 
    mutate(source = source_3)
  colnames(species_per_sample_3)[2] <- c("num_species")
  
  species_per_sample_crc <- colSums(metaphlan_df_crc[, -c(1, ncol(metaphlan_df_crc))] > 0) %>% 
    as.data.frame() %>% 
    rownames_to_column("sample_id") %>% 
    mutate(source = "CRCbiome")
  colnames(species_per_sample_crc)[2] <- c("num_species")
  
  
  df <- bind_rows(species_per_sample_1, species_per_sample_2, species_per_sample_3, species_per_sample_crc)
  
  # Calculate mean and sd for each group
  stats <- df %>%
    group_by(source) %>%
    summarize(
      mean = mean(num_species),
      sd = sd(num_species),
      y_position = min(num_species) - 80 # Position the text just above the lowest value
    )
  
  # Define colours for groups
  source_colors <- setNames(c("#F8766D", "#619CFF", "#00BA38", "#C77CFF"), 
                            c("CRCbiome", source_1, source_2, source_3))
  
  # Create a box plot
  box_plot <- ggplot(df, aes(x = source, y = num_species, fill = source)) + 
    geom_boxplot(outlier.colour = "red", outlier.shape = 8) + 
    geom_text(data = stats, aes(x = source, y = y_position, label = paste("Mean:", round(mean), "\nSD:", round(sd))),
              position = position_dodge(width = 0.75),
              vjust = -0.5,
              size = 4,
              fontface = "bold") +
    labs(x = "",
         y = "Number of species per sample", 
         fill = "") + 
    scale_fill_manual(values = source_colors) +
    theme_bw() +
    theme(
      text = element_text(size = 15, face = "bold"),  
      axis.title = element_text(size = 16), 
      axis.text = element_text(size = 14),  
      legend.title = element_text(size = 14),  
      legend.text = element_text(size = 14), 
      plot.title = element_text(size = 16, face = "bold")
    )
  
  # Wilcoxon test
  wilcoxon_data <- df %>% 
    select(-sample_id) %>% 
    mutate(source = factor(source))
  
  wilcoxon_res <- pairwise.wilcox.test(wilcoxon_data$num_species, wilcoxon_data$source, exact = FALSE, p.adjust.method = "none")
  
  return(list(box_plot = box_plot, wilcoxon_res = wilcoxon_res))
}

# calculate_species_abundances
# This function processes a MetaPhlAn data frame, extracting the species names and calculating the mean 
# relative abundance of each species across all samples. The data frame is expected to have samples as 
# columns and species as rows. The function then returns a data frame with species sorted by their mean 
# relative abundance in descending order.
#
# Args:
#   metaphlan_df: A data frame that has been preprocessed by process_metaphlan_table_rel_abundance
#                 with species as rows and samples as columns.
#
# Returns:
#   A data frame where each row corresponds to a species. The data frame contains the species names 
#   and their corresponding mean relative abundances, sorted in descending order of the mean abundance.
#
# Example:
#   sorted_species_abundances <- calculate_species_abundances(processed_metaphlan_df)
#
# Note:
#   The input data frame is expected to have a 'clade_name' column with taxonomic information from 
#   which the species name is extracted. The function assumes that the first column of the input data 
#   frame is 'clade_name', and this column will not be included in the mean abundance calculation.
#
calculate_species_abundances <- function(metaphlan_df) {
  # Extract the last part of the clade_name starting with 's__'
  metaphlan_df$clade_name <- sub(".*s__", "", metaphlan_df$clade_name)
  
  # Calculate the mean relative abundance for each species
  mean_abundances <- rowMeans(metaphlan_df[, -1])  # Exclude the first column (clade_name)
  
  # Combine the species names with their mean and sum abundances
  species_abundances <- data.frame(clade_name = metaphlan_df$clade_name, 
                                   mean_abundance = mean_abundances)
  
  # Sort by mean abundance in descending order and select the top 10
  sorted_species <- head(species_abundances[order(-species_abundances$mean_abundance), ], 10)
  
  # Return the sorted data frame
  return(list(all_species_abundances = species_abundances, sorted_species = sorted_species))
}

# plot_abundance
# This function plots relative abundance of taxa across multiple sources. It takes abundance data from 
# multiple sources along with a list of taxa sorted by abundance. 
# It requires abundance and sorted species data for each source, which are outputs from the 
# `calculate_species_abundances` function. The function merges the data and calculates the mean relative abundance for 
# each taxon, then plots the top taxa based on their mean abundance.
#
# Args:
#   all_abundances_a, all_abundances_b, all_abundances_c: DataFrames, containing abundance data for each taxon from three different sources.
#   sorted_species_a, sorted_species_b, sorted_species_c: DataFrames, with sorted taxa names sorted by abundance.
#   source_a, source_b, source_c: Strings, labels for the sources to be used in the plot legend (needs to start with "seed" as of now).
#
# Returns:
#   A ggplot object representing the bar plot of mean relative abundance of taxa.
#
# Example:
#   all_abundances_a <- calculate_species_abundances(metaphlan.sim5_a)$all_species_abundances
#   sorted_species_a <- calculate_species_abundances(metaphlan.sim5_a)$sorted_species
#    
#   all_abundances_b <- calculate_species_abundances(metaphlan.sim5_b)$all_species_abundances
#   sorted_species_b <- calculate_species_abundances(metaphlan.sim5_b)$sorted_species
#    
#   all_abundances_c <- calculate_species_abundances(metaphlan.sim5_c)$all_species_abundances
#   sorted_species_c <- calculate_species_abundances(metaphlan.sim5_c)$sorted_species
#
#   abundance_plot <- plot_abundance(all_abundances_a = all_abundances_a,
#                                    all_abundances_b = all_abundances_b,
#                                    all_abundances_c = all_abundances_c,
#                                    sorted_species_a = sorted_species_a,
#                                    sorted_species_b = sorted_species_b,
#                                    sorted_species_c = sorted_species_c,
#                                    source_a = "seed = 632741178", 
#                                    source_b = "seed = 584925030", 
#                                    source_c = "seed = 470267527")
plot_abundance <- function(all_abundances_a,
                           all_abundances_b,
                           all_abundances_c,
                           sorted_species_a,
                           sorted_species_b,
                           sorted_species_c,
                           source_a, 
                           source_b, 
                           source_c) {
  
  all_species <- sorted_species_a %>% 
    full_join(sorted_species_b, by = "clade_name") %>% 
    full_join(sorted_species_c, by = "clade_name") %>% 
    select(clade_name)
  
  df_a <- all_abundances_a %>% 
    filter(clade_name %in% all_species$clade_name) %>% 
    rename_with(~source_a, .cols = mean_abundance) 
  
  df_b <- all_abundances_b %>% 
    filter(clade_name %in% all_species$clade_name) %>% 
    rename_with(~source_b, .cols = mean_abundance) 
  
  df_c <- all_abundances_c %>% 
    filter(clade_name %in% all_species$clade_name) %>% 
    rename_with(~source_c, .cols = mean_abundance)
  
  data_long <- df_a %>% 
    full_join(df_b, by = "clade_name") %>% 
    full_join(df_c, by = "clade_name") %>% 
    pivot_longer(
      cols = starts_with("seed"), # Select columns that start with 'seed'
      names_to = "simulation",   # New column for simulation identifiers
      values_to = "mean_abundance"
    )
  
  # Define colours for groups
  source_colors <- setNames(c("#619CFF", "#00BA38", "#C77CFF"), 
                            c(source_a, source_b, source_c))
  
  # Plotting bar plot
  bar_plot <- ggplot(data_long, aes(x = reorder(clade_name, mean_abundance), 
                                    y = mean_abundance, 
                                    fill = as.factor(simulation))) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = source_colors) +
    theme_bw() +
    theme(
      text = element_text(size = 15, face = "bold"), 
      axis.title = element_text(size = 16),  
      axis.text = element_text(size = 12),  
      legend.title = element_text(size = 14), 
      legend.text = element_text(size = 16),  
      plot.title = element_text(size = 16, face = "bold")) + 
    labs(x = "Species", 
         y = "Mean Relative Abundance",
         title = "Top 10 Most Abundant Species",
         fill = "") +
    coord_flip()
  
  return(bar_plot)
}

# calculate_prevalence
# This function calculates the prevalence of taxa in a metaphlan dataset based on a specified presence threshold.
# The output includes the overall prevalence of each taxon and the prevalence of the top N taxa.
#
# Args:
#   df: DataFrame, the processed metagenomic abundance table from `process_metaphlan_table_rel_abundance_filtering` 
#       with relative abundances, where the first column is taxa names and the remaining columns are samples.
#   presence_threshold: Numeric, the minimum abundance value to consider a taxon as present in a sample (default is 0).
#   top_n: Integer, the number of top taxa to return based on prevalence (default is 10).
#
# Returns:
#   A list containing two elements:
#   - taxa_prevalence: A DataFrame with taxa names and their corresponding prevalence across samples.
#   - top_taxa_prevalence: A DataFrame of the top N taxa by prevalence.
#
# Example:
#   metaphlan_path <- "path/to/metaphlan_table.tsv"
#   processed_data <- process_metaphlan_table_rel_abundance_filtering(metaphlan_path)
#   prevalence_results <- calculate_prevalence(processed_data, presence_threshold = 0.01, top_n = 5)

calculate_prevalence <- function(df, presence_threshold = 0, top_n = 10) {
  # Ensure all columns except the first are numeric
  df[, -1] <- lapply(df[, -1], function(x) as.numeric(as.character(x)))
  
  # Extract the last part of the clade_name starting with 's__'
  df$clade_name <- sub(".*s__", "", df$clade_name)
  
  # Convert abundance data to presence/absence (1/0)
  presence_absence_df <- as.data.frame(lapply(df[, -1], function(x) as.integer(x > presence_threshold)))
  
  # Calculate prevalence for each taxon (row-wise sum divided by the number of columns)
  prevalence <- rowSums(presence_absence_df) / ncol(presence_absence_df)
  
  # Combine prevalence with the taxa names
  taxa_prevalence <- data.frame(clade_name = df$clade_name, prevalence = prevalence)
  
  # Sort by prevalence in descending order
  taxa_prevalence <- taxa_prevalence[order(-taxa_prevalence$prevalence), ]
  
  # Sort by prevalence in descending order and select the top 10
  top_taxa_prevalence <- head(taxa_prevalence[order(-taxa_prevalence$prevalence), ], top_n)
  
  return(list(taxa_prevalence = taxa_prevalence, top_taxa_prevalence = top_taxa_prevalence))
}

# plot_prevalence
# This function plots the prevalence of taxa across multiple sources. It requires prevalence 
# data from each source along with a list of taxa sorted by prevalence. The function merges this data to calculate 
# and visualize the prevalence of each taxon across different sources. It produces a bar plot showing the taxa with 
# the highest prevalence across the provided sources.
#
# Args:
#   all_prevalences_a, all_prevalences_b, all_prevalences_c: DataFrames, the prevalence data for taxa from three 
#       different sources. These should be outputs from `calculate_prevalence`.
#   sorted_species_a, sorted_species_b, sorted_species_c: DataFrames, sorted taxa names from three sources, which 
#       are outputs from `calculate_prevalence`.
#   source_a, source_b, source_c: Strings, the labels for the data sources to be used in the plot legend.
#
# Returns:
#   A ggplot object representing the bar plot of the prevalence of taxa.
#
# Example:
#   prevalence_data_a <- calculate_prevalence(processed_metaphlan_data_a)
#   all_prevalences_a <- calculate_prevalence(metaphlan.sim5_a)$taxa_prevalence
#   sorted_species_a <- calculate_prevalence(metaphlan.sim5_a)$top_taxa_prevalence
# 
#   all_prevalences_b <- calculate_prevalence(metaphlan.sim5_b)$taxa_prevalence
#   sorted_species_b <- calculate_prevalence(metaphlan.sim5_b)$top_taxa_prevalence
# 
#   all_prevalences_c <- calculate_prevalence(metaphlan.sim5_c)$taxa_prevalence
#   sorted_species_c <- calculate_prevalence(metaphlan.sim5_c)$top_taxa_prevalence
# 
#   prevalence_plot <- plot_prevalence(all_prevalences_a = all_prevalences_a,
#                                      all_prevalences_b = all_prevalences_b,
#                                      all_prevalences_c = all_prevalences_c,
#                                      sorted_species_a = sorted_species_a,
#                                      sorted_species_b = sorted_species_b,
#                                      sorted_species_c = sorted_species_c,
#                                      source_a = "seed = 632741178",
#                                      source_b = "seed = 584925030",
#                                      source_c = "seed = 470267527")

plot_prevalence <- function(all_prevalences_a,
                            all_prevalences_b,
                            all_prevalences_c,
                            sorted_species_a,
                            sorted_species_b,
                            sorted_species_c,
                            source_a,
                            source_b,
                            source_c) {
  
  all_species <- sorted_species_a %>% 
    full_join(sorted_species_b, by = "clade_name") %>% 
    full_join(sorted_species_c, by = "clade_name") %>% 
    select(clade_name)
  
  df_a <- all_prevalences_a %>% 
    filter(clade_name %in% all_species$clade_name) %>% 
    rename_with(~source_a, .cols = prevalence) 
  
  df_b <- all_prevalences_b %>% 
    filter(clade_name %in% all_species$clade_name) %>% 
    rename_with(~source_b, .cols = prevalence) 
  
  df_c <- all_prevalences_c %>% 
    filter(clade_name %in% all_species$clade_name) %>% 
    rename_with(~source_c, .cols = prevalence)
  
  data_long <- df_a %>% 
    full_join(df_b, by = "clade_name") %>% 
    full_join(df_c, by = "clade_name") %>% 
    pivot_longer(
      cols = starts_with("seed"), # Select columns that start with 'seed'
      names_to = "simulation",   # New column for simulation identifiers
      values_to = "prevalence"
    )
  
  # Define colours for groups
  source_colors <- setNames(c("#619CFF", "#00BA38", "#C77CFF"), 
                            c(source_a, source_b, source_c))
  
  # Plotting
  p <- ggplot(data_long, aes(x = reorder(clade_name, prevalence), y = prevalence, fill = simulation)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = source_colors) +
    theme_bw() +
    theme(
      text = element_text(size = 15, face = "bold"), 
      axis.title = element_text(size = 16),  
      axis.text = element_text(size = 12),  
      legend.title = element_text(size = 14), 
      legend.text = element_text(size = 16), 
      plot.title = element_text(size = 16, face = "bold")) + 
    labs(x = "Species", 
         y = "Prevalence", 
         title = paste("Top 10 Species by Prevalence"),
         fill = "") +
    coord_flip()
  
  return(p)
}

# This function is a part of the GeometricMorphometricsMix R package (https://github.com/fruciano/GeometricMorphometricsMix)
#' Test the significance of the adjusted Rand index
#'
#' Permutation test of the adjusted Rand index, which quantifies the level of agreement
#' between two partitions (e.g., two schemes of classification of the same individuals obtained with two methods)
#'
#' The adjusted Rand index (Hubert and Arabie 1985), is an adjusted for chance version of the Rand index (Rand 1971).
#' The adjusted Rand index has an expected value of zero in the case of random partitions,
#' and values approaching one as the two partitions become more similar to each other
#' (with one being perfect match of the classification).
#'  This function implements the permutation test proposed by Qannari et al. (2014)
#' to obtain a p value against the null hypothesis of independence of the two partitions.
#'
#'
#' This function is useful in various contexts, such as in integrative taxonomy
#' when comparing the classification of individual specimens obtained using different data
#' (e.g., sequence data and morphometric data).
#' For an example of the application of this technique with the classification obtained with genetic data and
#' morphometric data for multiple traits, see Fruciano et al. 2016.
#'
#'
#' @section Notice:
#' The function requires internally the package mclust.
#'
#' @param A,B numerical or character vectors reflecting the assignment of individual observations to groups
#' @param perm number of permutations
#' @return The function outputs a vector with the adjusted Rand index and the p value obtained from the permutation test
#'
#' @section Citation:
#' If you use this function in the context of integrative taxonomy or similar
#' (comparison of classification/unsupervised clustering with biological data), please cite all the papers in the references
#' (otherwise, please use the relevant citations for the context).
#'
#' @references Fruciano C, Franchini P, Raffini F, Fan S, Meyer A. 2016. Are sympatrically speciating Midas cichlid fish special? Patterns of morphological and genetic variation in the closely related species Archocentrus centrarchus. Ecology and Evolution 6:4102-4114.
#' @references Hubert L, Arabie P. 1985. Comparing partitions. Journal of Classification 2:193-218.
#' @references Qannari EM, Courcoux P, Faye P. 2014. Significance test of the adjusted Rand index. Application to the free sorting task. Food Quality and Preference 32, Part A:93-97.
#' @references Rand WM. 1971. Objective Criteria for the Evaluation of Clustering Methods. Journal of the American Statistical Association 66:846-850.
#'
#' @examples
#' library(mclust)
#' set.seed(123)
#'
#' irisBIC = mclustBIC(iris[,-5])
#' mclustBIC_classification = summary(irisBIC,iris[,-5])$classification
#' original_classification = iris[,5]
#' # This is one of the examples in the package mclust
#' # Here a classification algorithm is used on the iris dataset
#'
#' adjustedRandIndex(mclustBIC_classification, original_classification)
#' # The mclust package allows computing the adjusted Rand index
#' # which quantifies the agreement between the original (correct) classification
#' # and the one obtained with the algorithm.
#' # However, it is not clear whether the adjusted Rand index is "large enough"
#' # compared to the null hypothesis of independence between the two classification schemes
#'
#' adjRand_test(mclustBIC_classification, original_classification, perm = 999)
#' # For that, we use the function adjRand_test, which performs the permutation test
#' # of Qannari et al. 2014 (in this case p<0.001, as 1000 permutations have been used).
#'
#' adjRand_test(original_classification, original_classification, perm = 999)
#' # As it can be seen, in the ideal case of the exact same grouping,
#' # the adjusted Rand index takes a value of 1 (which is obviously significant)
#'
#'
#' @export
adjRand_test=function(A, B, perm=999) {
  if (length(A)!=length(B)) { stop("A and B should have the same length") }
  # Make sure that the two groups of partitions have the same length
  
  ARIo=mclust::adjustedRandIndex(A, B)
  # Observed adjusted Rand index
  
  Aperm=lapply(seq_len(perm), function(X) sample(A, length(A), replace=FALSE))
  Bperm=lapply(seq_len(perm), function(X) sample(B, length(B), replace=FALSE))
  # Generate permuted samples
  
  ARIperm=unlist(lapply(seq_len(perm),
                        function(i) mclust::adjustedRandIndex(Aperm[[i]], Bperm[[i]])))
  # Compute adjusted Rand index for the permuted samples
  
  m=mean(ARIperm); v=var(ARIperm)
  # compute mean and variance of the permuted samples
  
  NARI=(ARIperm-m)/sqrt(v)
  # compute NARI (normalized ARI)
  
  NARIo=(ARIo-m)/sqrt(v)
  # compute observed NARI
  
  p_value=(length(which(NARI>NARIo))+1)/(perm+1)
  # Compute p value as proportion of permuted NARI larger than the observed
  
  Results=c(ARIo, p_value)
  names(Results)=c("Adjusted_Rand_index", "p_value")
  return(Results)
}

# abundance_profile
# This function simulates a differential abundance profile based on median coverage data of metagenomic samples. 
# It divides the samples into two groups, calculates the prevalence of MAGs in each group, 
# selects metagenomic assembled genomes (MAGs) randomly to be differentially abundant, and adjusts coverage values to 
# create a simulated abundance profile, considering the set prevalence threshold and an adjustment factor.
#
# Args:
#   coverage_df: DataFrame, median coverage data for metagenomic samples, with samples as rows and features as columns.
#   rand_mags_n: Integer, number of random MAGs to select for differential abundance.
#   seed: Integer, seed for random number generation to ensure reproducibility.
#   adjustment_factor: Numeric, factor by which to adjust the coverage of selected MAGs in group 2.
#
# Returns:
#   profiles: DataFrame, relative abundance profiles with adjusted coverage for the selected differentially abundant MAGs.
#
# Example:
#   df <- abundance_profile(coverage_df = median_coverage_genomes, 
#                           rand_mags_n = 10, 
#                           seed = 123, 
#                           adjustment_factor = 2.0)

abundance_profile <- function(coverage_df = median_coverage_genomes, rand_mags_n, seed = 123, adjustment_factor){
  
  # Dividing samples into 2 groups
  set.seed(seed)
  coverage_df.grops <- coverage_df %>% 
    mutate(group = rep(2, n())) %>% # Start with all 2s for group
    mutate(group = replace(group, sample(row_number(), n() / 2), 1)) %>% # Replace half with 1s randomly
    mutate(group = as.factor(group)) %>%  # Convert group to a factor for analysis purposes
    select(sample_id, group, everything())
  
  # Set the prevalence threshold
  prevalence_threshold <- 0.1
  
  # Function to calculate prevalence for a given group
  calculate_proportions <- function(data_subset, threshold) {
    sapply(data_subset, function(feature) {
      mean(feature > threshold, na.rm = TRUE)
    })
  }
  
  # Split the data by group
  group1_data <- coverage_df.grops[coverage_df.grops$group == 1, -c(1:2)]
  group2_data <- coverage_df.grops[coverage_df.grops$group == 2, -c(1:2)]
  
  # Calculate prevalence for each group
  proportions_group1 <- calculate_proportions(group1_data, prevalence_threshold)
  proportions_group2 <- calculate_proportions(group2_data, prevalence_threshold)
  
  # Find MAGs with prevalence above 10% in both groups
  prevalent_mags <- names(proportions_group1)[proportions_group1 > 0.1 & proportions_group2 > 0.1]
  
  # Randomly choose n random mags from the prevalent mags to be differential abundant
  set.seed(seed)
  random_mags <- prevalent_mags %>% 
    sample(size = rand_mags_n)
  
  # Convert data to long format
  coverage_df.long <- coverage_df.grops %>%
    pivot_longer(
      cols = -c(sample_id, group), # exclude the id and group from the gathering
      names_to = "MAG",
      values_to = "median_coverage")
  
  # Finding the randomly selected MAGs in the data frame. If the MAG is found in the vector DA = 1, else DA = 0
  mag_in_vector <- coverage_df.long$MAG %in% random_mags
  coverage_df.long$DA <- as.integer(mag_in_vector)
  
  # Adjust coverage values
  coverage_df.long$median_coverage_adjusted <- 
    ifelse(coverage_df.long$DA == 1 & coverage_df.long$group == 2, 
           coverage_df.long$median_coverage * adjustment_factor, 
           coverage_df.long$median_coverage) # median_coverage_adjusted stays the same as median_coverage if conditions are not fulfilled
  
  # Create relative abundance values
  profiles <- coverage_df.long %>% 
    group_by(sample_id) %>%
    mutate(sum_measurements = sum(median_coverage_adjusted, na.rm = TRUE)) %>% # Calculate the sum of measurements for each sample
    mutate(rel_abundance = (median_coverage_adjusted / sum_measurements)) %>% # Calculate the relative abundance by total reads and sum of measurements
    ungroup() %>% 
    select(-c(sum_measurements))
  
  return(profiles)
}

# write_profiles_to_file
# This function writes abundance profiles to files suitable for CAMISIM.
# It takes a data frame of relative abundance profiles and writes them to individual TSV files
# for each sample. These files are formatted to be compatible with CAMISIM. It also
# creates two strings listing the TSV files for each group, which can be used in CAMISIM configuration files.
#
# Args:
#   df.long: DataFrame, the long format of abundance profiles, which is output from the abundance_profile function.
#   out_dir_tsv: String, the directory path where the TSV files will be saved.
#
# Returns:
#   A list containing three elements:
#   - gr1_string: A string with paths to the group 1 TSV files concatenated together, separated by commas.
#   - group1_tsv_files: A vector containing the paths to the TSV files for group 1.
#   - gr2_string: A string with paths to the group 2 TSV files concatenated together, separated by commas.
#   - group2_tsv_files: A vector containing the paths to the TSV files for group 2.
#
# Example:
#   relative_abundance_df <- abundance_profile(...)
#   df.long <- abundance_profile(coverage_df = median_coverage_genomes, 
#                                rand_mags_n = 50, 
#                                seed = 123, 
#                                adjustment_factor = 10.0)
#   out_dir <- "/path/to/somewhere/"
#   writing_files <- write_profiles_to_file(df.long, out_dir)

write_profiles_to_file <- function(df.long, out_dir_tsv){
  
  # Make correct format for CAMISIM
  profiles <- df.long %>% 
    mutate(MAG = paste0(MAG, ".0")) %>%  # Add .0 to match requirements for CAMISIM
    mutate(rel_abundance = as.numeric(format(rel_abundance, scientific = FALSE)))
  
  # Writing to file
  profiles_list <- profiles %>%
    select(sample_id, MAG, rel_abundance) %>%
    group_by(sample_id) %>%
    group_split()
  
  lapply(profiles_list, function(x){
    file_path <- paste0(out_dir_tsv, x$sample_id[1], ".tsv")
    write_tsv(x[, -1], file = file_path, col_names = FALSE)
  })
  
  # Extracting sample ids for each group
  
  cleaned_path <- sub("/$", "", out_dir_tsv)
  
  all_tsv_files <- list.files(path = cleaned_path, full.names = TRUE, recursive = TRUE)
  
  # Extract files for group 1
  
  samples_gr1 <- profiles %>%
    select(group, sample_id) %>%
    filter(group == 1) %>%
    distinct(sample_id)
  
  group1_sample_ids <- samples_gr1$sample_id
  
  # Initialize a vector to hold the TSV files for group 1
  group1_tsv_files <- c()
  
  # Loop through each sample ID and check if there is a matching TSV file
  for(sample_id in group1_sample_ids) {
    # Create a pattern to match the file name
    pattern <- paste0(sample_id, ".tsv")
    
    # Find matching TSV files
    matching_files <- grep(pattern, all_tsv_files, value = TRUE)
    
    # Append the matching files to the group1_tsv_files vector
    group1_tsv_files <- c(group1_tsv_files, matching_files)
  }
  
  # Create listing of files to use for config file
  gr1_string <- paste(group1_tsv_files, collapse = ",")
  
  # Extract files for group 2
  
  samples_gr2 <- profiles %>%
    select(group, sample_id) %>%
    filter(group == 2) %>%
    distinct(sample_id)
  
  group2_sample_ids <- samples_gr2$sample_id
  
  # Initialize a vector to hold the TSV files for group 2
  group2_tsv_files <- c()
  
  # Loop through each sample ID and check if there is a matching TSV file
  for(sample_id in group2_sample_ids) {
    # Create a pattern to match the file name
    pattern <- paste0(sample_id, ".tsv")
    
    # Find matching TSV files
    matching_files <- grep(pattern, all_tsv_files, value = TRUE)
    
    # Append the matching files to the group2_tsv_files vector
    group2_tsv_files <- c(group2_tsv_files, matching_files)
  }
  
  # Create listing of files to use for config file
  gr2_string <- paste(group2_tsv_files, collapse = ",")
  
  return(list(gr1_string = gr1_string, 
              gr2_string = gr2_string, 
              group1_tsv_files = group1_tsv_files, 
              group2_tsv_files = group2_tsv_files))
}

# profiles_pcoa_plot
# This function performs PCoA based on abundance profiles and plots the results. It prepares the data for PCoA,
# calculates distance matrices, runs the PCoA, and then generates a PCoA plot, a scree plot, and a box plot of PCoA scores
# for each dimension. The function returns these plots along with the PCoA coordinates and the percentage of 
# variation explained by each PCoA axis.
#
# Args:
#   df.long: DataFrame, a long-format DataFrame of relative abundance profiles.
#   method: String, the method to be used for distance calculation in PCoA (default is "bray").
#   k: Integer, the number of principal coordinates to calculate (default is 200).
#
# Returns:
#   A list containing:
#   - pcoa_plot: A ggplot object of the PCoA plot.
#   - scree_plot: A ggplot object of the scree plot.
#   - box_plot: A ggplot object of the box plot for PCoA scores.
#   - coordinates_with_groups: A DataFrame with PCoA coordinates and group assignments.
#   - explained_var: A vector of the percentage of variation explained by each PCoA axis.
#
# Example:
#   abundance_profiles_long <- read.csv("path/to/abundance_profiles_long.csv")
#   pcoa <- profiles_pcoa_plot(df.long)
#   pcoa_plot <- pcoa$pcoa_plot

profiles_pcoa_plot <- function(df.long, method = "bray", k = 200){
  
  # Process data to correct format for PCoA
  wide_data <- df.long %>% 
    pivot_wider(names_from = MAG, 
                values_from = rel_abundance, 
                id_cols = c(sample_id, group))
  
  # Calculate distance matrix
  dist.mat <- vegdist(wide_data[, -c(1:2)], method = method)
  
  # Perform PCoA
  pcoa_mat <- cmdscale(dist.mat, k = k, eig = TRUE)
  
  # Extract coordinates and set column names
  coordinates <- as.data.frame(pcoa_mat$points)
  colnames(coordinates) <- paste("PCo", seq_len(k), sep = "")
  
  # Extract explained variation for each axis
  percent_explained <- 100 * pcoa_mat$eig / sum(pcoa_mat$eig)
  pretty_pe <- round(percent_explained, 2)
  
  # Create labels for the plots
  labs <- c(glue("PCo1 ({pretty_pe[1]}%)"),
            glue("PCo2 ({pretty_pe[2]}%)"))
  
  # Process data
  coordinates_with_groups <- bind_cols(group = wide_data$group, coordinates)
  
  factor_levels <- levels(factor(coordinates_with_groups$group))
  colorVector <- c("#619CFF", "#00BA38")
  source_colors <- setNames(colorVector, factor_levels)
  
  # Create PCoA plot
  pcoa_plot <- ggplot(coordinates_with_groups[, c(1:3)], 
                      aes(x = PCo1, y = PCo2, color = group)) +
    geom_point(size = 2) +
    scale_color_manual(values = source_colors) +
    guides(shape = "none") +
    labs(x = labs[1], y = labs[2], color = "Group") + 
    theme_bw() + 
    theme(
      text = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  # theme(
  #   text = element_text(size = 6, face = "bold"),
  #   axis.title = element_text(size = 6),
  #   axis.text = element_text(size = 6),
  #   legend.title = element_text(size = 6),
  #   legend.text = element_text(size = 6),
  #   plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  # )
  
  # Create scree plot
  scree_plot <- ggplot(tibble(pe = percent_explained),
                       aes(x = 1:length(pe), y = pe)) +
    geom_bar(stat = "identity") +
    labs(x = "PCoA axis", y = "Percent explained by axis",
         title = "Scree plot") +
    theme_bw() +
    theme(
      text = element_text(size = 12, face = "bold"),  
      axis.title = element_text(size = 12),  
      axis.text = element_text(size = 12),  
      legend.title = element_text(size = 12),  
      legend.text = element_text(size = 12),  
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  
    )
  # theme(
  #   text = element_text(size = 6, face = "bold"), 
  #   axis.title = element_text(size = 6),  
  #   axis.text = element_text(size = 6),  
  #   legend.title = element_text(size = 6),  
  #   legend.text = element_text(size = 6), 
  #   plot.title = element_text(size = 16, face = "bold", hjust = 0.5) 
  # )
  
  # Create box plot of PCoA scores for each dimension
  coordinates.long <- gather(coordinates_with_groups, key = "PCo", value = "Score", -group) %>% 
    mutate(PCo = as.factor(PCo))
  
  # Reordering PCos for plotting
  # Extract numbers from the factor levels, convert to numeric
  group_levels <- as.numeric(gsub("PCo", "", levels(coordinates.long$PCo)))
  
  # Create an ordered factor with the levels in the correct order
  coordinates.long$PCo <- factor(coordinates.long$PCo, levels = levels(coordinates.long$PCo)[order(group_levels)])
  
  # Create the boxplot with ggplot2
  box_plot <- ggplot(coordinates.long, aes(x = PCo, y = Score, fill = group)) +
    geom_boxplot() +
    scale_fill_manual(values = source_colors) +
    labs(x = "",
         y = "Score",
         fill = "Group") +
    theme_bw() +
    theme(
      text = element_text(size = 12, face = "bold"),   
      axis.title = element_text(size = 12),   
      axis.text = element_text(size = 12),  
      legend.title = element_text(size = 12),   
      legend.text = element_text(size = 12),   
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)   
    )
  # theme(
  #   text = element_text(size = 6, face = "bold"),   
  #   axis.title = element_text(size = 6),   
  #   axis.text = element_text(size = 6),  
  #   axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
  #   legend.title = element_text(size = 6),   
  #   legend.text = element_text(size = 6),   
  #   plot.title = element_text(size = 16, face = "bold", hjust = 0.5)   
  # )
  
  return(list(pcoa_plot = pcoa_plot, 
              scree_plot = scree_plot, 
              box_plot = box_plot, 
              coordinates_with_groups = coordinates_with_groups,
              explained_var = percent_explained))
}

# export_data_for_umap
# This function processes a list of abundance profile data frames and exports them to CSV files in a specified output 
# directory to be used in UMAP. Each file is named with a provided suffix corresponding to each data frame. The data is pivoted to a 
# wide format required for UMAP analysis, with MAGs as columns and relative abundance as values.
#
# Args:
#   df_list: List of DataFrames, where each DataFrame contains long-format abundance profile data.
#   suffixes: Vector of strings, suffixes to be appended to each output CSV file corresponding to the data frames in df_list.
#   output_dir: String, the directory path where the CSV files will be saved.
#
# Returns:
#   CSV files are written to the specified directory, one for each DataFrame in df_list.
#
# Example:
#   abundance_profiles <- list(df1, df2, df3)  # replace with actual data frame variables
#   suffixes <- c("timepoint1", "timepoint2", "timepoint3")  # replace with relevant suffixes
#   output_directory <- "path/to/output_directory/"
#   export_data_for_umap(df_list = abundance_profiles, suffixes = suffixes, output_dir = output_directory)
#
# No value is returned by the function; it writes files directly to the filesystem.

export_data_for_umap <- function(df_list, suffixes, output_dir) {
  # Create the output directory if it does not exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Loop through the list of data frames and corresponding suffixes
  for (i in seq_along(df_list)) {
    df_long <- df_list[[i]]
    suffix <- suffixes[i]
    
    # Pivot the data frame
    wide_data <- df_long %>%
      pivot_wider(names_from = MAG, values_from = rel_abundance, id_cols = c(sample_id, group))
    
    # Construct the output file name
    csv_name <- paste0("data_", suffix, ".csv")
    file_path <- file.path(output_dir, csv_name)
    
    # Write to CSV
    write.csv(wide_data[, -c(1:2)], file_path, row.names = FALSE)
  }
}

# plotting_umap
# This function creates UMAP plots for visualizing high-dimensional abundance profiles grouped by sample categories. 
# It processes the data to extract the necessary information for UMAP plotting, reads the UMAP embeddings from CSV files, 
# generates individual UMAP plots for each file, and compiles them into a grid layout. It also adds annotations for the 
# total number of MAGs in the title.
#
# Args:
#   files: Vector of strings, paths to CSV files containing UMAP embeddings.
#   groups_df: DataFrame, containing group information for each sample.
#   n_col: Integer, number of columns in the plot grid.
#   n_row: Integer, number of rows in the plot grid.
#   grid_plot_title: String, title for the entire grid plot (optional, defaults to an empty string).
#
# Returns:
#   A list with the following elements:
#   - plot_grid: A ggplot object of the UMAP plots arranged in a grid.
#   - annotated_grid: The grid plot with annotations.
#
# Example:
#   files_vector <- c("umap_embeddings_timepoint1.csv", "umap_embeddings_timepoint2.csv")
#   umap_plots <- plotting_umap(files = files_vector, groups_df = df.long,
#                               n_col = 2, n_row = 2, grid_plot_title = "UMAP Visualization")
#
# Note: Ensure that UMAP embeddings have been generated and saved as CSV files before running this function.

plotting_umap <- function(files, groups_df, n_col = 5, n_row = 4, grid_plot_title = ""){
  
  # Titles for the first five plots
  # titles <- c("2x", "5x", "10x", "20x", "25x")
  titles <- c("10 MAGs 2x", "20 MAGs 5x", "50 MAGs 10x", "100 MAGs 20x")
  
  # Initialize an empty list to store plots
  plot_list <- list()
  
  # Initialize a counter for the plots
  plot_counter <- 1
  
  # Process data to extract group variable
  wide_data <- groups_df %>% 
    pivot_wider(names_from = MAG, 
                values_from = rel_abundance, 
                id_cols = c(sample_id, group))
  
  for (file in files){
    
    # Read the embedding
    embedding <- read.csv(file) %>%
      mutate(group = wide_data$group)
    
    # Define colors for plotting
    source_colors <- setNames(c("#619CFF", "#00BA38"), c("1", "2"))
    
    # Determine the title for the current plot
    plot_title <- if (plot_counter <= length(titles)) titles[plot_counter] else ""
    
    # Plot UMAP results
    umap_plot <- ggplot(embedding, aes(x=UMAP_1, y=UMAP_2, color = group)) +
      geom_point(size = 2) +
      scale_color_manual(values = source_colors) +
      labs(x = "UMAP 1",
           y = "UMAP 2",
           color = "Group",
           title = plot_title) + # Set the title here
      theme_bw() +
      theme(
        text = element_text(size = 6, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
      coord_fixed(ratio = 1)
    
    # Append the plot to the list
    plot_list[[length(plot_list) + 1]] <- umap_plot
    
    # Increment the plot counter
    plot_counter <- plot_counter + 1
  }
  
  plot_grid <- ggarrange(
    plotlist = plot_list,
    ncol = n_col, 
    nrow = n_row,
    common.legend = TRUE,
    legend = "bottom",
    labels = "AUTO"
  )
  
  plot_grid <- annotate_figure(plot_grid,
                               top = text_grob(grid_plot_title, face = "bold", size = 16))
  
  annotated_grid <- annotate_figure(plot_grid,
                                    left = grobTree(
                                      textGrob("10 MAGs", y=0.88, rot = 90, gp=gpar(fontface="bold", fontsize = 12)),
                                      textGrob("20 MAGs", y=0.65, rot = 90, gp=gpar(fontface="bold", fontsize = 12)),
                                      textGrob("50 MAGs", y=0.40, rot = 90, gp=gpar(fontface="bold", fontsize = 12)),
                                      textGrob("100 MAGs", y=0.16, rot = 90, gp=gpar(fontface="bold", fontsize = 12))
                                    )
  )
  
  return(list(grid = plot_grid, annotated_grid = annotated_grid))
  
}

# HC_umap
# This function reads a CSV file with UMAP embeddings, performs hierarchical clustering based on the specified method, 
# and calculates the Adjusted Rand Index (ARI) along with its p-value to compare the predicted clusters against a 
# given grouping variable.
#
# Args:
#   embedding_file: String, the file path to the CSV containing UMAP embeddings.
#   groups: Vector, contains the true group assignments of samples for comparison.
#   hc_method: String, hierarchical clustering method to use (default: "average").
#
# Returns:
#   A list with the following components:
#   - ari: Numeric, Adjusted Rand Index indicating the similarity between true groups and clusters.
#   - p_val: Numeric, p-value indicating the statistical significance of the ARI.
#
# Example:
#   embedding_csv <- "path/to/umap_embeddings.csv"
#   clustering_results <- HC_umap(embedding_file = embedding_csv, groups = true_groups, hc_method = "average")
#   print(paste("ARI:", clustering_results$ari))
#   print(paste("p-value:", clustering_results$p_val))

HC_umap <- function(embedding_file, groups, hc_method = "average"){
  
  embedding <- read.csv(embedding_file) %>% 
    mutate(group = groups) %>% 
    select(group, everything())
  
  # Making a distance matrix for clustering
  dist_mat <- dist(embedding[, -1], method = "euclidean") 
  
  # Executing hierarchical clustering with average linkage
  hcl.obj <- hclust(dist_mat, method = hc_method)
  
  # Predicted clusters
  predicted_clusters <- cutree(hcl.obj, k = 2)
  
  # Adjusted rand index
  set.seed(123)
  ari <- as.numeric(adjRand_test(as.integer(embedding$group), predicted_clusters, perm = 1000)[1])
  
  set.seed(123)
  p_val <- as.numeric(adjRand_test(as.integer(embedding$group), predicted_clusters, perm = 1000)[2])
  
  return(list(ari = ari, p_val = p_val))
  
}

# kmeans_umap
# This function takes UMAP embedding, applies k-means clustering to identify clusters, 
# and computes the Adjusted Rand Index (ARI) to evaluate the agreement between the clustering result 
# and a pre-defined grouping.
#
# Args:
#   embedding_file: String, the path to the CSV file containing UMAP embedding coordinates.
#   groups: Vector or Factor, containing the true group labels for each sample in the dataset.
#
# Returns:
#   A list containing:
#   - ari: The Adjusted Rand Index, assessing the similarity between the k-means clustering output and the true groups.
#   - p_val: The p-value for the ARI, evaluating the significance of the index.
#
# Example:
#   umap_embeddings_csv <- "path/to/umap_embeddings.csv"
#   true_group_labels <- c("Group1", "Group1", "Group2", "Group2")  # Replace with actual group labels
#   kmeans_results <- umap_kmeans(embedding_file = umap_embeddings_csv, groups = true_group_labels)
#   print(paste("ARI:", kmeans_results$ari))
#   print(paste("p-value:", kmeans_results$p_val))

kmeans_umap <- function(embedding_file, groups){
  
  embedding <- read.csv(embedding_file) %>% 
    mutate(group = groups) %>% 
    select(group, everything())
  
  # Perform k means clustering
  kmeans.obj <- kmeans(embedding[, -1], centers = 2, iter.max = 10, nstart = 10)
  
  predicted_clusters <- kmeans.obj$cluster
  
  # Adjusted rand index
  set.seed(123)
  ari <- as.numeric(adjRand_test(as.integer(embedding$group), predicted_clusters, perm = 1000)[1])
  
  set.seed(123)
  p_val <- as.numeric(adjRand_test(as.integer(embedding$group), predicted_clusters, perm = 1000)[2])
  
  return(list(ari = ari, p_val = p_val))
  
}

# process_metaphlan_table_pcoa_profile_based_sim
# This function further processes a MetaPhlAn table that has already been processed by `process_metaphlan_table_rel_abundance_filtering`
# to fit the format required for PCoA analysis. The function transposes the table to have samples as rows and species as columns, 
# and ensures all abundance values are numeric.
#
# Args:
#   metaphlan_df: DataFrame, the processed MetaPhlAn table from `process_metaphlan_table_rel_abundance_filtering` with relative abundances.
#
# Returns:
#   metaphlan_table_proc: DataFrame, the prepared MetaPhlAn table ready for PCoA analysis with 'sample_id' as the first column.
#
# Example:
#   metaphlan_filtered <- process_metaphlan_table_rel_abundance_filtering(metaphlan_path)
#   metaphlan_prepared_for_pcoa <- process_metaphlan_table_pcoa_profile_based_sim(metaphlan_filtered)

process_metaphlan_table_pcoa_profile_based_sim <- function(metaphlan_df) {
  
  # Transpose the table to fit cmdscale() requirements with samples as rows and species as columns
  metaphlan_table <- t(metaphlan_df) %>% 
    as.data.frame()
  
  # Assign column names and create a 'sample_id' column while removing row names
  colnames(metaphlan_table) <- metaphlan_table[1, ]
  metaphlan_table <- rownames_to_column(metaphlan_table, "sample_id")
  
  # Remove the first row which contains clade names, as it's not needed for analysis
  metaphlan_table <- metaphlan_table[-1, ]
  
  # Store 'sample_id' for later use and convert all abundance values to numeric type
  sample_ids <- metaphlan_table$sample_id
  metaphlan_table <- data.frame(lapply(metaphlan_table[, -1], as.numeric), sample_id = sample_ids)
  
  # Re-introduce 'sample_id' at the beginning and select all columns for output
  metaphlan_table_proc <- metaphlan_table %>%
    select(sample_id, everything())
  
  # Return the prepared table
  return(metaphlan_table_proc)
}

# metaphlan_pcoa_plot
#
# This function conducts a Principal Coordinate Analysis (PCoA) on a processed MetaPhlAn table. It computes the distance matrix,
# executes the PCoA, extracts the explained variance for each axis, and creates a PCoA plot, a scree plot to visualize the variance
# explained by each principal coordinate, and a box plot showing the distribution of PCoA scores for each dimension. It's tailored for
# abundance data with samples as rows and species as columns after being processed by `process_metaphlan_table_pcoa_profile_based_sim`.
#
# Args:
#   df.wide: DataFrame, a wide-format MetaPhlAn table ready for PCoA analysis, typically the output of the
#            `process_metaphlan_table_rel_abundance_filtering` and `process_metaphlan_table_pcoa_profile_based_sim` functions.
#   method: String, the method used to calculate the distance matrix for PCoA (default is "bray").
#   k: Integer, the number of dimensions to compute for the PCoA (default is 200).
#
# Returns:
#   A list containing:
#   - pcoa_plot: A ggplot object representing the PCoA plot.
#   - scree_plot: A ggplot object representing the scree plot.
#   - box_plot: A ggplot object representing the box plot of PCoA scores.
#   - coordinates_with_groups: A DataFrame combining group assignments with PCoA coordinates.
#   - explained_var: A numeric vector with the percentage of variance explained by each PCoA axis.
#
# Example:
#   metaphlan_processed <- process_metaphlan_table_rel_abundance_filtering(metaphlan_path)
#   metaphlan_processed <- process_metaphlan_table_pcoa_profile_based_sim(metaphlan_processed)
#   pcoa_results <- metaphlan_pcoa_plot(df.wide = metaphlan_processed, method = "bray", k = 20)
#   print(pcoa_results$pcoa_plot)
#   print(pcoa_results$scree_plot)
#   print(pcoa_results$box_plot)

metaphlan_pcoa_plot <- function(df.wide, method = "bray", k = 200){
  
  # Calculate distance matrix
  dist.mat <- vegdist(df.wide[, -c(1:2)], method = method)
  
  # Perform PCoA
  pcoa_mat <- cmdscale(dist.mat, k = k, eig = TRUE)
  
  # Extract coordinates and set column names
  coordinates <- as.data.frame(pcoa_mat$points)
  colnames(coordinates) <- paste("PCo", seq_len(k), sep = "")
  
  # Extract explained variation for each axis
  percent_explained <- 100 * pcoa_mat$eig / sum(pcoa_mat$eig)
  pretty_pe <- round(percent_explained, 2)
  
  # Create labels for the plots
  labs <- c(glue("PCo1 ({pretty_pe[1]}%)"),
            glue("PCo2 ({pretty_pe[2]}%)"))
  
  # Process data
  coordinates_with_groups <- bind_cols(group = df.wide$group, coordinates)
  
  factor_levels <- levels(factor(coordinates_with_groups$group))
  colorVector <- c("#619CFF", "#00BA38")
  source_colors <- setNames(colorVector, factor_levels)
  
  # Create PCoA plot
  pcoa_plot <- ggplot(coordinates_with_groups[, c(1:3)], 
                      aes(x = PCo1, y = PCo2, color = group)) +
    geom_point(size = 2) +
    scale_color_manual(values = source_colors) +
    guides(shape = "none") +
    labs(x = labs[1], y = labs[2], color = "Group") + 
    theme_bw() + 
    theme(
      text = element_text(size = 12, face = "bold"),   
      axis.title = element_text(size = 12),   
      axis.text = element_text(size = 12),  
      legend.title = element_text(size = 12),   
      legend.text = element_text(size = 12),   
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)   
    )
  # theme(
  #   text = element_text(size = 6, face = "bold"),   
  #   axis.title = element_text(size = 6),   
  #   axis.text = element_text(size = 6),  
  #   legend.title = element_text(size = 6),   
  #   legend.text = element_text(size = 6),   
  #   plot.title = element_text(size = 16, face = "bold", hjust = 0.5)   
  # )
  
  # Create scree plot
  scree_plot <- ggplot(tibble(pe = percent_explained),
                       aes(x = 1:length(pe), y = pe)) +
    geom_bar(stat = "identity") +
    #geom_col() +
    labs(x = "PCoA axis", y = "Percent explained by axis",
         title = "Scree plot") +
    theme_bw() +
    theme(
      text = element_text(size = 12, face = "bold"),   
      axis.title = element_text(size = 12),   
      axis.text = element_text(size = 12),  
      legend.title = element_text(size = 12),   
      legend.text = element_text(size = 12),   
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)   
    )
  # theme(
  #   text = element_text(size = 6, face = "bold"),   
  #   axis.title = element_text(size = 6),   
  #   axis.text = element_text(size = 6),  
  #   legend.title = element_text(size = 6),   
  #   legend.text = element_text(size = 6),   
  #   plot.title = element_text(size = 16, face = "bold", hjust = 0.5)   
  # )
  
  # Create box plot of PCoA scores for each dimension
  coordinates.long <- gather(coordinates_with_groups, key = "PCo", value = "Score", -group) %>% 
    mutate(PCo = as.factor(PCo))
  
  # Reordering PCos for plotting
  # Extract numbers from the factor levels, convert to numeric
  group_levels <- as.numeric(gsub("PCo", "", levels(coordinates.long$PCo)))
  
  # Create an ordered factor with the levels in the correct order
  coordinates.long$PCo <- factor(coordinates.long$PCo, levels = levels(coordinates.long$PCo)[order(group_levels)])
  
  # Create the boxplot with ggplot2
  box_plot <- ggplot(coordinates.long, aes(x = PCo, y = Score, fill = group)) +
    geom_boxplot() +
    scale_fill_manual(values = source_colors) +
    labs(x = "",
         y = "Score",
         fill = "Group") +
    theme_bw() +
    theme(
      text = element_text(size = 12, face = "bold"),   
      axis.title = element_text(size = 12),   
      axis.text = element_text(size = 12),  
      legend.title = element_text(size = 12),   
      legend.text = element_text(size = 12),   
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)   
    )
  # theme(
  #   text = element_text(size = 6, face = "bold"),   
  #   axis.title = element_text(size = 6),   
  #   axis.text = element_text(size = 6),  
  #   axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
  #   legend.title = element_text(size = 6),   
  #   legend.text = element_text(size = 6),   
  #   plot.title = element_text(size = 16, face = "bold", hjust = 0.5)   
  # )
  
  return(list(pcoa_plot = pcoa_plot, 
              scree_plot = scree_plot, 
              box_plot = box_plot, 
              coordinates_with_groups = coordinates_with_groups,
              explained_var = percent_explained))
}

# export_data_for_umap_metaphlan
# This function exports a list of wide-format MetaPhlAn data frames to CSV files for UMAP analysis. 
# It is designed to take the output data frames from the MetaPhlAn analysis pipeline and save them with 
# corresponding suffixes to distinguish between different conditions or time points.
#
# Args:
#   df_list: List of DataFrames, the wide-format MetaPhlAn tables to be exported.
#   suffixes: Vector of strings, the suffixes for the output file names, corresponding to each data frame.
#   output_dir: String, the directory path where the CSV files will be saved.
#
# Returns:
#   CSV files are written to the specified directory, one for each DataFrame in df_list.
#
# Example:
#   metaphlan_tables <- list(df1, df2)  # df1, df2,... are data frames from MetaPhlAn analysis
#   suffixes <- c("s1", "s2")
#   output_dir <- "/path/to/somewhere/"
#   export_data_for_umap_metaphlan(df_list = metaphlan_tables, suffixes = suffixes, output_dir = output_dir)

export_data_for_umap_metaphlan <- function(df_list, suffixes, output_dir) {
  
  # Create the output directory if it does not exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Loop through the list of data frames and corresponding suffixes
  for (i in seq_along(df_list)) {
    df_wide <- df_list[[i]]
    suffix <- suffixes[i]
    
    # Construct the output file name
    csv_name <- paste0("data_", suffix, ".csv")
    file_path <- file.path(output_dir, csv_name)
    
    # Write to CSV
    write.csv(df_wide[, -c(1:2)], file_path, row.names = FALSE)
  }
}

# plotting_umap_metaphlan
# This function creates UMAP plots. 
# It processes the data to extract the necessary information for UMAP plotting, reads the UMAP embeddings from CSV files, 
# generates individual UMAP plots for each file, and compiles them into a grid layout. It also adds annotations for the 
# total number of MAGs in the title.
#
# Args:
#   files: Vector of strings, paths to CSV files containing UMAP embeddings.
#   groups_df: DataFrame, containing group information for each sample.
#   n_col: Integer, number of columns in the plot grid.
#   n_row: Integer, number of rows in the plot grid.
#   grid_plot_title: String, title for the entire grid plot (optional, defaults to an empty string).
#
# Returns:
#   A list with the following elements:
#   - plot_grid: A ggplot object of the UMAP plots arranged in a grid.
#   - annotated_grid: The grid plot with annotations.
#
# Example:
#   files_vector <- c("umap_embeddings_timepoint1.csv", "umap_embeddings_timepoint2.csv")
#   umap_plots <- plotting_umap(files = files_vector, groups_df = df.long,
#                               n_col = 2, n_row = 2, grid_plot_title = "UMAP Visualization")
#
# Note: Ensure that UMAP embeddings have been generated and saved as CSV files before running this function.

plotting_umap_metaphlan <- function(files, groups, n_col = 2, n_row = 1, grid_plot_title = ""){
  
  titles <- c("20 MAGs 5x", "50 MAGs 10x")
  
  # Initialize an empty list to store plots
  plot_list <- list()
  
  # Initialize a counter for the plots
  plot_counter <- 1
  
  
  for (file in files){
    
    # Read the embedding
    embedding <- read.csv(file) %>%
      mutate(group = groups)
    
    factor_levels <- levels(factor(groups))
    colorVector <- c("#619CFF", "#00BA38")
    source_colors <- setNames(colorVector, factor_levels)
    
    # Determine the title for the current plot
    plot_title <- if (plot_counter <= length(titles)) titles[plot_counter] else ""
    
    # Plot UMAP results
    umap_plot <- ggplot(embedding, aes(x=UMAP_1, y=UMAP_2, color = group)) +
      geom_point(size = 2) +
      scale_color_manual(values = source_colors) +
      labs(x = "UMAP 1",
           y = "UMAP 2",
           color = "Group",
           title = plot_title) + # Set the title here
      theme_bw() +
      theme(
        text = element_text(size = 6, face = "bold"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
      )
    
    # Append the plot to the list
    plot_list[[length(plot_list) + 1]] <- umap_plot
    
    # Increment the plot counter
    plot_counter <- plot_counter + 1
  }
  
  plot_grid <- ggarrange(
    plotlist = plot_list,
    ncol = n_col, 
    nrow = n_row,
    common.legend = TRUE,
    legend = "bottom"
  )
  
  plot_grid <- annotate_figure(plot_grid,
                               top = text_grob(grid_plot_title, face = "bold", size = 16))
  
  annotated_grid <- annotate_figure(plot_grid,
                                    left = grobTree(
                                      textGrob("20 MAGs", y=0.65, rot = 90, gp=gpar(fontface="bold", fontsize = 12)),
                                      textGrob("50 MAGs", y=0.40, rot = 90, gp=gpar(fontface="bold", fontsize = 12)))
  )
  
  return(list(grid = plot_grid, annotated_grid = annotated_grid))
  
}