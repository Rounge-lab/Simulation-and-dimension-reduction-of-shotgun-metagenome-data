# Imports
library(tidyverse)

# Path to sample ids to use in MetaPhlAn pipeline
sample_dir <- "/PATH/metaphlan_input/fastq/some_directory"

# Get the name of the samples
sample_id <- list.files(path = sample_dir) %>% 
  str_remove_all(pattern = ".fastq.gz") %>% 
  str_remove_all(pattern = "_R1") %>% 
  str_remove_all(pattern = "_R2") %>% 
  unique()

# Create data frame for samples.tsv for MetaPhlAn
samples.tbl <- tibble(sample_id) %>% 
  mutate(sample_number = as.integer(sub("sample_", "", sample_id))) %>% # Extract numeric part and convert to integer
  arrange(sample_number) %>% # Order the table by the numeric values
  select(-sample_number)

# Write sample.tbl to .tsv
setwd("/PATH/metaphlan_input")

write.table(samples.tbl,
            "samples.tsv", # Remember to change the name
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)

