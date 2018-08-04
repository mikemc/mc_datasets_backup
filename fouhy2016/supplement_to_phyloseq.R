library(tidyverse)
library(magrittr)
library(phyloseq)

dotenv::load_dot_env("../.env")
script_path <- getwd()
data_path <- file.path(Sys.getenv("DATA_DIR"),
    "fouhy2016")

## Load the relative abundance table

# First, download the supplementary file with the relative abundance info.
# "Additional file 1: Table S1. Percentage relative abundance of expected
# species detected in the mock cell DNA. (DOCX 21 kb)"
si_url <- "https://static-content.springer.com/esm/art%3A10.1186%2Fs12866-016-0738-z/MediaObjects/12866_2016_738_MOESM1_ESM.docx"
fn <- file.path(data_path, "supplementary_files", "additional_file_1.docx")
dir.create(file.path(data_path, "supplementary_files"), recursive = TRUE)
download.file(si_url, fn)

# I then created a csv version of the relative-abundance table in the file
# "table_s1.csv" by converting the table in addtional_file_1.docx to a a csv
# table, deleting empty rows, and replacing the table caption in the first cell
# of the header row with the column name "Sample"
st <- readr::read_csv(file.path(script_path, "table_s1.csv"))
names(st) <- names(st) %>%
    str_replace(" ", "_") %>%
    str_replace("Escherichia/Shigella_coli", "Escherichia_coli")
st <- st %>%
    select(-Sample) %>%
    data.frame(row.names = st$Sample) %>%
    otu_table(taxa_are_rows = FALSE)

## Parse the sample metadata from the sample names
# Primer Sequencer Extraction Solvent
template <- "(\\S{5}+(?: Deg)?+) (PGM|Miseq) (Qiagen|RBB) (PBS|Glycerol)"
sam <- sample_names(st) %>%
    str_match(template) %>%
    as_tibble
names(sam) <- c("Sample", "Primer", "Sequencer", "Extraction", "Solvent")
# Adjust the metadata variables to match format of sample_metadata.R
sam <- sam %>%
    mutate(Primer = case_when(
            Primer == "V1-V2" ~ "V12",
            Primer == "V1-V2 Deg" ~ "V12d",
            Primer == "V4-V5" ~ "V45",
            ),
        Type = "Cells", 
        Description = Sample)
# Convert to a `sample_data` object
sam <- sam %>%
    select(-Sample) %>%
    data.frame(row.names = sam$Sample) %>%
    sample_data
# Merge into a phyloseq object
ps <- phyloseq(st, sam)
# Create abbreviated sample names
sample_names(ps) <- sample_data(ps) %$%
    paste0(
        str_sub(Type, 1, 1), 
        Primer,
        str_sub(Sequencer, 1, 1), 
        str_sub(Extraction, 1, 1),
        str_sub(Solvent, 1, 1)
        )
# Save
dir.create(file.path(data_path, "final"))
saveRDS(ps, file.path(data_path, "final", 
        paste0(Sys.Date(), "_phyloseq_from_supplement.Rds")))
