#!/usr/bin/env Rscript

print("getting clade prevalence")

# env ---------------------------------------------------------------------

suppressPackageStartupMessages(library("tidyverse"))

save.image("tmp.RData")
# load data ---------------------------------------------------------------

mp <- read_tsv(snakemake@input[["profile"]], skip = 1, show_col_types = FALSE)

db <- read_tsv(snakemake@params[["SGB_species"]], col_names = FALSE, show_col_types = FALSE) %>% 
  separate_rows(X2, sep = ",")


# Get prevalence ----------------------------------------------------------

clades_presence <-
  mp %>% 
  filter(grepl("s__", clade_name)) %>% 
  filter(!grepl("t__", clade_name)) %>% 
  pivot_longer(-clade_name) %>% 
  inner_join(db %>% 
               select(clade = X1,
                      clade_name = X2), by = "clade_name") %>% 
  group_by(clade, name) %>% 
  summarize(pres = any(value > 0)) %>% 
  ungroup() %>% 
  group_by(clade) %>% 
  summarize(n_obs = sum(pres)) %>% 
  ungroup()

clades_presence %>% 
  write_tsv(snakemake@output[["full_clade_list"]])

clades_presence %>% 
  filter(n_obs > snakemake@params[["min_obs"]]) %>% 
  slice_head(n = 20) %>% 
  write_tsv(snakemake@output[["filtered_clade_list"]])
