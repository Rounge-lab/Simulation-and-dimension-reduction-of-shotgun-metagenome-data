# Imports
library(tidyverse)

# Define path to MAG directory
mag_dir <- "/PATH/CRCbiome/datasets/metagenome/MAGs/raw/genomes"

# Get the list of files
fasta_files <- list.files(path = mag_dir, full.names = TRUE)

# Get the names of the MAGs
mag_names <- list.files(path = mag_dir) %>% 
  str_remove_all(pattern = ".fasta") %>% 
  lapply(function(x) paste0(as.character(x), ".0")) %>% 
  unlist()

# Create data frame for genome_to_id file for CAMISIM
genomes_to_id <- data_frame(mag_names, fasta_files)

# Create data frame for metadata.tsv for CAMISIM
genome_ID <- mag_names
OTU <- 1:length(mag_names)
NCBI_ID <- rep(2, times = length(mag_names))
novelty_category <- rep("known_strain", times = length(mag_names))

metadata <- data.frame(genome_ID, OTU, NCBI_ID, novelty_category)

# Write genome_to_id table to .tsv (skipping the column names)
write.table(genomes_to_id,
            "genomes_to_id_all_MAGs.tsv",
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

# Write metadata table to .tsv (skipping the column names)
write.table(metadata,
            "metadata_all_MAGs.tsv",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
