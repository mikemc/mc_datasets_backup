library(tidyverse)
library(magrittr)
library(phyloseq)

dotenv::load_dot_env("../.env")
script_path <- getwd()
data_path <- file.path(Sys.getenv("DATA_DIR"),
    "bakker2018")

out_path <- file.path(data_path, "from_NB", "Amplicon_sequencing_library_1")
list.files(out_path)

# Load the metadata from the SRA
run_info <- readxl::read_excel(file.path(data_path, "from_NB",
        "Run_details.xlsx"))
names(run_info) %>% sort


# Load the sequence table
list.files(out_path, pattern = "rds")
# sts <- list.files(out_path, pattern = "rds", full.names = T) %>%
#     map(readRDS)
st <- file.path(out_path, "seqtab.nochim.library1.rds") %>%
    readRDS

sam <- run_info %>%
    select(Sample_Name, SRA_run = Run, Mock,
        PCR_conditions, Replicate, Center_name = Center_Name,
        Library_prep_method = Library_Prep_Method)

OTU <- otu_table(st, taxa_are_rows = FALSE)
SAM <- sam %>%
    data.frame(row.names = sam$SRA_run) %>%
    sample_data
ps <- phyloseq(OTU, SAM)
sample_names(ps) <- sample_data(ps)$Sample_Name

# Add taxonomy
cns <- c("SV", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus",
    "Species")
tax <- list.files(out_path, pattern = "taxa", full.names = T) %>%
    map_dfr(readr::read_tsv, skip = 1, col_names = cns) %>%
    distinct()
# One SV appears twice - for some reason, assigment didn't work well one time
tax %>%
    group_by(SV) %>%
    filter(n() > 1)
# Let's get rid of the one without genus assignment
sv <- tax %>%
    group_by(SV) %>%
    filter(n() > 1) %$%
    SV[1]
tax <- tax %>%
    filter((SV != sv) | !is.na(Genus))
# Merge phyloseq
mat <- tax %>%
    select(-SV) %>%
    as("matrix")
rownames(mat) = tax$SV
TAX <- tax_table(mat)
ps <- merge_phyloseq(ps, TAX)

# Merge sequences
sequences <- Biostrings::DNAStringSet(taxa_names(ps))
names(sequences) <- taxa_names(ps)
ps <- merge_phyloseq(ps, sequences)

# Save
dir.create(file.path(data_path, "final"))
saveRDS(ps, file.path(data_path, "final", 
        paste0(Sys.Date(), "_phyloseq_from_NB.Rds")))


