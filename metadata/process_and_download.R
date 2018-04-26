library(stringr)
library(tidyverse)
library(magrittr)


#
data_path <- "/home/michael/data/costea2017"
git_path <- "/home/michael/active_research/metagenomics_calibration"
metadata_path <- file.path(git_path, "costea2017", "metadata")

## Sample metadata
# The file sample_description.xlsx is a version of Supplementary Data 3
# reformatted for easier reading into R
phase_dfs <- map(1:3, 
    ~readxl::read_xlsx(file.path(metadata_path, "sample_description_mod"), 
        sheet = .))
sam_df <- bind_rows(phase_dfs, .id = "Phase")
## Sequence data info
seq_df <- readr::read_tsv(file.path(metadata_path, "PRJEB14847.tsv"))

##
# Next, connect the two so that I can download the files of the samples I wish

seq_df %>% names
seq_df %>% select(experiment_alias, sample_alias, library_name)

seq_df$sample_alias %>% str_subset("BYQ")
seq_df$library_name %>% str_subset("BYQ")

seq_df %>% select(experiment_alias, sample_alias, library_name) %>%
    filter(str_detect(library_name, "BYQ"))

# Seems like "library_name" is the field that contains the sample name as used
# in the sample description file for Phase III, but perhaps not for the other
# phases

seq_df$library_name %>% str_subset("AWF")
seq_df$sample_alias %>% str_subset("A1")

seq_df %>% select(experiment_alias, sample_alias, library_name) %>%
    filter(str_detect(sample_alias, "058"))

# sample_alias for Phase 1 has a dash, and for phase 2 has an underscore, between the 


## Start with Phase 3, since that's what I'm most interested in and it's
## easiest to understand.
ph3 <- seq_df %>% select(library_name, base_count, fastq_bytes, fastq_ftp) %>% 
    filter(str_detect(library_name, "BYQ"))
files <- ph3$fastq_ftp %>% str_split(";", simplify = TRUE)
ph3 <- ph3 %>%
    mutate(ftp1 = files[, 1], ftp2 = files[,2])
# Get the (combined for reads1 + reads2) file size of each sample in GB
file_sizes <- ph3$fastq_bytes %>%
    str_split(";", simplify = TRUE) %>%
    as_tibble %>%
    mutate(V1 = as.numeric(V1), V2 = as.numeric(V2),
        file_size = (V1 + V2) / 1e9) %$%
    file_size
# Get the sample names from the library names
sample_names <- ph3$library_name %>%
    str_match("(.*)_DA") %>%
    {.[,2]}
# Build a simplied df with Sample info and file info
ph3 <- ph3 %>%
    mutate(file_size = file_sizes, Sample = sample_names) %>%
    select(Sample, file_size, ftp1, ftp2)
ph3 <- sam_df %>% filter(str_detect(Sample, "BYQ")) %>%
    left_join(ph3, by = "Sample")
#
ph3 %>% filter(Individual == "Mock-only")

seq_df %>% filter(str_detect(library_name, "BYQ_AAAT"))

seq_df %>% filter(str_detect(library_name, "BYQ")) %>%
    select(library_name, submitted_ftp)

seq_df %>%
    select(library_name, submitted_ftp)

## Download phase 3 files

ph3 <- ph3 %>%
    mutate(rsync1 = paste0('rsync://', ftp1),
        rsync2 = paste0('rsync://', ftp2))
ph3$ftp1 %>% head
ph3$rsync1 %>% head
all_rsync_paths <- c(ph3$rsync1, ph3$rsync2)
sum(ph3$file_size)
# A pretty large download. Let's start with just the mock only samples

mock_rsync_paths <- ph3 %>%
    filter(Individual == "Mock-only") %$%
    c(rsync1, rsync2)

# Download the files with rsync
dir.create(file.path(data_path, "reads"), recursive = TRUE)
commands <- paste("rsync", "--progress", "--copy-links", "--times",
    "--verbose", mock_rsync_paths, file.path(data_path, "reads"))
walk(commands, system)
