library(tidyverse)
library(rentrez)
library(xml2)

dotenv::load_dot_env("../.env")
script_path <- getwd()
data_path <- file.path(Sys.getenv("DATA_DIR"),
    "costea2017")

## Get the sample metadata from the supplemental file provided by the authors

# Supplementary Data 1 : Protocol descriptors
# https://media.nature.com/original/nature-assets/nbt/journal/v35/n11/extref/nbt.3960-S4.xlsx
#
# Supplementary Data 2 : Members and composition of mock community
# https://media.nature.com/original/nature-assets/nbt/journal/v35/n11/extref/nbt.3960-S5.xlsx
#
# Supplementary Data 3 : Sample description
# https://media.nature.com/original/nature-assets/nbt/journal/v35/n11/extref/nbt.3960-S6.xlsx

urls <- c("https://media.nature.com/original/nature-assets/nbt/journal/v35/n11/extref/nbt.3960-S4.xlsx",
    "https://media.nature.com/original/nature-assets/nbt/journal/v35/n11/extref/nbt.3960-S5.xlsx",
    "https://media.nature.com/original/nature-assets/nbt/journal/v35/n11/extref/nbt.3960-S6.xlsx")
names(urls) <- c("protocol_descriptors.xlsx", "mock_composition.xlsx",
    "sample_description.xlsx")

fns <- file.path(data_path, "supplementary_files", names(urls))
dir.create(file.path(data_path, "supplementary_files"))
walk2(urls, fns, download.file)

## Import the metadata for the Phase 3 experiment

# The metadata for the three experimental phases of the Costea2017 experiment
# (Phases I, II, and III, which I'll call 1, 2, and 3) is given on a single
# Excel sheet in separate rectagular ranges
cell_ranges <- c("A3:B192", "D3:F77", "H3:J32")
names(cell_ranges) <- 1:3
sam_ls <- cell_ranges %>%
    map(~ readxl::read_xlsx(fns[3], range = .))
names(sam_ls) <- names(cell_ranges)
sam <- sam_ls %>%
    bind_rows(.id = "Phase")

## Get the NCBI/ENA run accessions to connect samples to their sequence data
# https://www.ebi.ac.uk/ena/data/view/PRJEB14847
# https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJEB14847
# In NCBI: Accession: PRJEB14847 ID: 391591
sra_id <- entrez_link(dbfrom = "bioproject", id = "391591", db="sra") %>% 
    {.$links$bioproject_sra}
# Fetching the full records is slow (takes a min or two); the full sra records
# are large, about 100Mb when stored as text.
xml <- entrez_fetch(db = "sra", id = sra_id, rettype = "xml") %>%
    read_xml
# write_lines(xml, "/tmp/sra_records.xml")

# The Library name can be used to connect back to the Sample name in the study
# metadata, and the run accession can be used to download the sequence data
# using NCBI's SRA toolkit or to find the ENA download links
tb <- tibble(
    Library = xml %>%
        xml_find_all(".//LIBRARY_NAME") %>%
        xml_text,
    Run_accession = xml %>%
        xml_find_all(".//RUN/IDENTIFIERS/PRIMARY_ID") %>%
        xml_text,
    )

# Getting just the shorter xml strings in the summaries is much faster, but
# more annoying to parse because the xml string is mangled (escape sequences
# are inserted for all tags). I may switch to using these in the future

## SCRATCH
# sra <- entrez_summary(db = "sra", id = sra_id)
# xml0 <- sra %>% 
#     map_chr("expxml")
# xml1 <- xml0 %>%
#     {glue::collapse(.)}
# xml0[1] %>%
#     str_replace_all("&lt;", "<") %>%
#     str_replace_all("&gt;", ">") %>%
#     str_sub(3, -3)
## \SCRATCH

# Next, I want to add the run accessions to the sample metadata table. But the
# sample names in the nucleotide servers and in the sample metadata can include
# some extra characters, causing them not to match up, and these differ for the
# three Phases.

## SCRATCH
# # To get the library names from the sample names, we need to remove the middle
# # portion of the Phase 1 names, and add "_DA" to the end
# phase1_template <- "([A-B]1_)[0-9]+(?:_2)?_(\\w+)"
# sam <- sam %>%
#     mutate(Library = case_when(
#             Phase == 1 ~ str_match(Sample, phase1_template) %>% 
#                 {paste0(.[,2], .[,3])},
#             Phase == 2 ~ Sample,
#             Phase == 3 ~ Sample,
#             )) %>%
#     mutate(Library = paste0(Library, "_DA"))
# # Actually, this doesn't work
## /SCRATCH

# I'm just interested in the Phase 3 data for now, for which its easier to
# specify the relationship between the Library and Sample names---the Library
# name is just the Sample name with "_DA" tacked on. So let's just save the
# phase 3 metadata and SRA info.
ph3 <- sam %>%
    filter(Phase == 3) %>%
    mutate(Library = paste0(Sample, "_DA")) %>%
    left_join(tb, by = "Library")
readr::write_csv(ph3, file.path(script_path, "sample_metadata.csv"))
