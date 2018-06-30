# Process (munge) the metadata and download the sequence data for the
# Costea2017 Phase 3 experiment. This script was originally used to download
# the sequence data on the BRC cluster, but is no longer preferred.

library(tidyverse)
library(magrittr)

dotenv::load_dot_env("../../.env")
script_path <- getwd()
data_path <- file.path(Sys.getenv("DATA_DIR"),
    "costea2017")
metadata_path <- file.path(script_path)
aspera_path <- Sys.getenv("ASCP_PATH") %>%
    str_split(pattern = "\\|") %>%
    {.[[1]][1]} %>%
    dirname %>%
    dirname

# Study in the ENA at https://www.ebi.ac.uk/ena/data/view/PRJEB14847
# Accession info downloaded by choosing columns (with "Select columns" link)
# and downloading a tsv file ("Text" link) and saved as "PRJEB14847.tsv".

## Sample metadata
# The file sample_description.xlsx is a version of Supplementary Data 3
# reformatted for easier reading into R
phase_dfs <- map(1:3, 
    ~readxl::read_xlsx(file.path(metadata_path, "sample_description_mod.xlsx"), 
        sheet = .))
sam_df <- bind_rows(phase_dfs, .id = "Phase")
## Sequence data info
seq_df <- readr::read_tsv(file.path(metadata_path, "PRJEB14847.tsv"))
# Get the (combined for reads1 + reads2) file size of each sample in GB
file_sizes <- seq_df$fastq_bytes %>%
    str_split(";", simplify = TRUE) %>%
    as_tibble %>%
    mutate(V1 = as.numeric(V1), V2 = as.numeric(V2),
        file_size = (V1 + V2) / 1e9) %$%
    file_size
seq_df <- seq_df %>%
    mutate(fastq_GB = file_sizes)

## Next, connect the two so that I can download the files of the samples I wish
# Start with Phase 3, since that's what I'm most interested in and it's
# easiest to understand.
#
# Seems like "library_name" is the field that contains the sample name as used
# in the sample description file for Phase III, but perhaps not for the other
# phases
ph3 <- seq_df %>%
    filter(str_detect(library_name, "BYQ"))
sample_names <- ph3$library_name %>%
    str_match("(.*)_DA") %>%
    {.[,2]}
ph3 <- ph3 %>%
    mutate(Sample = sample_names) %>%
    left_join(sam_df, by="Sample")

## Download with ascp (aspera connect command line tool)
aspera_urls <- ph3$fastq_aspera %>% str_split(";", simplify=TRUE) %>% c
commands <- paste(
    file.path(aspera_path, "bin/ascp"),
    "-QT -l 300m -P33001 -i", 
    file.path(aspera_path, "etc/asperaweb_id_dsa.openssh"),
    paste0("era-fasp@", aspera_urls),
    file.path(data_path, "reads")
    )
walk(commands, system)
# For testing the commands:
system(commands[5])

## Check dowloaded files against md5sums
downloads <- aspera_urls %>% 
    str_extract("ERR[0-9]*_[1-2]\\.fastq\\.gz") %>%
    file.path(data_path, "reads", .)
md5sums_expected <- ph3$fastq_md5 %>% str_split(";", simplify=TRUE) %>% c
md5sums_actual <- tools::md5sum(downloads)
# Fraction of files that were successfully downloaded and gave an md5
mean(!is.na(md5sums_actual))
# Check that the downloaded files match the expected md5
all(md5sums_expected == md5sums_actual, na.rm = TRUE)

# ## Download with wget:
# ftp_urls <- ph3$fastq_ftp %>% str_split(";", simplify=TRUE) %>% c
# dir.create(file.path(data_path, "reads"), recursive = TRUE)
# commands <- paste("wget", "-P", file.path(data_path, "reads"), 
#     paste0("ftp://", ftp_urls)
# walk(commands, system)


#########################################################################
## Poking around the sample data

# seq_df %>% names
# seq_df %>% select(experiment_alias, sample_alias, library_name)
#
# seq_df %>% select(experiment_alias, sample_alias, library_name) %>%
#     filter(str_detect(sample_alias, "058"))
#
# seq_df$sample_alias %>% str_subset("BYQ")
# seq_df$library_name %>% str_subset("BYQ")
#
# seq_df %>% select(experiment_alias, sample_alias, library_name) %>%
#     filter(str_detect(library_name, "BYQ"))
#
# seq_df$library_name %>% str_subset("AWF")
# seq_df$sample_alias %>% str_subset("A1")
#
#
# # sample_alias for Phase 1 has a dash, and for phase 2 has an underscore, between the 
#
# ph3 %>% filter(Individual == "Mock-only")
#
# seq_df %>% filter(str_detect(library_name, "BYQ_AAAT"))
#
# seq_df %>% filter(str_detect(library_name, "BYQ")) %>%
#     select(library_name, submitted_ftp)
#
# seq_df %>%
#     select(library_name, submitted_ftp)
#
