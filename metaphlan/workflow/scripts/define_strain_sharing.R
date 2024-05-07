#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(cutpointr))

set_strain_sharing_threshold <- function(dists_file, sample_file, min_samples = 50, vc_thresholds_file, plot_path) {
  
  samples <- read_tsv(sample_file, col_names = TRUE, show_col_types = FALSE)
  vc_thresholds <- read_tsv(vc_thresholds_file, col_names = TRUE, show_col_types = FALSE)
  
  training_data <-
    read_tsv(dists_file, col_names = FALSE, show_col_types = FALSE) %>% 
    rename(sampleID_1 = X1, sampleID_2 = X2) %>% 
    left_join(samples %>% select(sampleID_1 = sample_id,
                            subjectID_1 = participant_id,
                            sample_type_1 = sample_type), by = "sampleID_1") %>% 
    left_join(samples %>% select(sampleID_2 = sample_id,
                            subjectID_2 = participant_id,
                            sample_type_2 = sample_type), by = "sampleID_2") %>% 
    mutate(same_individual = case_when(subjectID_1 == subjectID_2 ~ "same_individual",
                                       TRUE ~ "different_individual") %>% 
             factor(levels = c("same_individual", "different_individual")),
           clade = basename(dists_file) %>% str_remove("_dists.tsv")) %>% 
    ## Only consider comparisons of FU samples from the same individual
    filter(sample_type_1 != "Baseline" | same_individual %in% "different_individual",
           sample_type_2 != "Baseline" | same_individual %in% "different_individual") %>% 
    group_by(subjectID_1, subjectID_2) %>% 
    ## Only one comparison per pair of individuals
    slice_head(n = 1) %>% 
    ungroup()

  if (str_remove(training_data$clade[1], "t__") %in% vc_thresholds$SGB) {
    vc_threshold <- vc_thresholds$used_nGD_score_percentile[ vc_thresholds$SGB == str_remove(training_data$clade[1], "t__")][1]
  } else {
    vc_threshold <- NA
  }
  
  if (sum(training_data$same_individual == "same_individual") < min_samples) {
    tmp <-
      tibble(optimal_cutpoint = NA, youden = NA, AUC = NA)  
  } else {
    tmp <-
      training_data %>% 
      cutpointr(x = X3, class = same_individual, method = maximize_metric, metric = youden, silent = TRUE) %>% 
      select(optimal_cutpoint, youden, AUC)  
  }
  
  tmp <- 
    tmp %>% 
    mutate(quantile_5pc = training_data %>% filter(same_individual != "same_individual") %>% pull(X3) %>% quantile(0.05),
           n_within = training_data %>% filter(same_individual == "same_individual") %>% nrow(),
           n_between = training_data %>% filter(same_individual != "same_individual") %>% nrow(),
           vc_threshold = vc_threshold,
           threshold = case_when(n_within >= min_samples ~ min(optimal_cutpoint, quantile_5pc),
                                 n_within < min_samples ~ min(vc_threshold, quantile_5pc, na.rm = TRUE)),
           threshold_method = case_when(n_within >= min_samples & optimal_cutpoint <= quantile_5pc ~ "Youden",
                                        n_within >= min_samples & optimal_cutpoint > quantile_5pc ~ "percentile - Y",
                                        n_within < min_samples & quantile_5pc < vc_threshold & !is.na(vc_threshold) ~ "percentile - VC",
                                        n_within < min_samples & quantile_5pc > vc_threshold & !is.na(vc_threshold) ~ "VC",
                                        n_within < min_samples & is.na(vc_threshold) ~ "percentile")) %>% 
    mutate(clade = training_data$clade[1]) %>% 
    select(clade, everything())
  
  tmp_plt <-
    training_data %>% 
    ggplot(aes(x = X3, fill = same_individual)) +
    geom_density() +
    geom_vline(aes(x = NULL, xintercept = value, color = name), 
               data = tmp %>% pivot_longer(c(optimal_cutpoint, quantile_5pc, vc_threshold)), linetype = 2)
    
  ggsave(filename = plot_path, plot = tmp_plt)
  
  tmp
}

get_samples_with_shared_strains <- function(dists, threshold) {
  
  read_tsv(dists, col_names = FALSE, show_col_types = FALSE) %>% 
    rename(sample_1 = X1,
           sample_2 = X2,
           strain_dist = X3) %>% 
    filter(strain_dist <= threshold)
}

threshold <- 
  set_strain_sharing_threshold(dists_file = snakemake@input[["dists"]],
                               sample_file = snakemake@input[["sample_info"]],
                               min_samples = snakemake@params[["min_obs"]],
                               vc_thresholds_file = snakemake@params[["vc_thresholds"]],
                               plot_path = snakemake@output[["dist_dist"]]) 
threshold %>%
  write_tsv(snakemake@output[["strain_sharing_threshold"]])

get_samples_with_shared_strains(dists = snakemake@input[["dists"]],
                                threshold = threshold %>% pull(threshold)) %>% 
  write_tsv(snakemake@output[["strain_sharing_samples"]])
