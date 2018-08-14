library(tidyverse)
library(rentrez)
library(xml2)

dotenv::load_dot_env("../.env")
script_path <- getwd()
data_path <- file.path(Sys.getenv("DATA_DIR"),
    "fouhy2016")

# https://www.ncbi.nlm.nih.gov/sra/SRP071776
# https://www.ncbi.nlm.nih.gov/bioproject/PRJNA315115

sra_id <- entrez_link(dbfrom = "bioproject", id = "315115", db="sra") %>% 
    {.$links$bioproject_sra}
xml <- entrez_fetch(db = "sra", id = sra_id, rettype = "xml") %>%
    read_xml
write_lines(xml, "/tmp/ncbi_records.xml")

# s <- entrez_search(db = "sra", term = "SRP071776[ACCN]", retmax=1e3)
# xml <- entrez_fetch(db = "sra", id = s$ids, rettype = "xml") %>%
#     read_xml

tb <- tibble(
    Library = xml %>%
        xml_find_all(".//LIBRARY_NAME") %>%
        xml_text,
    SRA_run = xml %>%
        xml_find_all(".//RUN/IDENTIFIERS/PRIMARY_ID") %>%
        xml_text,
    Mock = xml %>% 
        xml_find_all(".//SAMPLE/DESCRIPTION") %>%
        xml_text,
    Description = xml %>% 
        xml_find_all(".//DESIGN_DESCRIPTION") %>%
        xml_text,
    )
print(tb, n=Inf)

tb <- tb %>%
    mutate(Type = case_when(
            Mock %in% c("HM-280", "HM-281") ~ "Cells",
            Mock == "HM-782D" ~ "DNA"),
        Solvent = case_when(
            Mock == "HM-280" ~ "PBS",
            Mock == "HM-281" ~ "Glycerol")
        )

# Next, parse everything from the Library name, make new sample names. Also
# verify that the Type and Solvent are consistent with the library name. Note
# that "Deg" is short for degenerate and indicates the V1V2 degenerate primers.
#
# should also consider if there is a safer way to parse that strictly checks
# the hierarchical structure, making sure everything within an
# EXPERIMENT_PACKAGE goes together.

# Format is <16s Sequencer Primer Extraction Solvent> OR
# <16s Sequencer Primer 'Mock DNA'>

template <- paste0("16s (PGM|MiSeq) (v[1|4]v[2|5](?: Deg)?+) ",
    "(Q|RBB|Mock) (PBS|Glycerol|DNA)"
    )
sam <- tb$Library %>%
    str_match(template) %>%
    as_tibble
names(sam) <- c("Library", "Sequencer", "Primer", "Extraction", "Solvent")
sam <- sam %>%
    mutate(Extraction = case_when(
            Extraction == "Q" ~ "Qiagen",
            Extraction == "RBB" ~ "RBB"
            ),
        Solvent = ifelse(Solvent == "DNA", NA, Solvent),
        Primer = case_when(
            Primer == "v1v2" ~ "V12",
            Primer == "v1v2 Deg" ~ "V12d",
            Primer == "v4v5" ~ "V45",
            )
        )
sam <- sam %>%
    select(Library, Primer, Sequencer, Extraction, Solvent) %>%
    left_join(tb, by = "Library",
        suffix = c("", ".tb"))
all.equal(sam$Solvent, sam$Solvent.tb)
# Create abbreviated sample names
sam <- sam %>%
    mutate(Sample = paste0(
            str_sub(Sequencer, 1, 1), 
            Primer,
            str_sub(Extraction, 1, 1),
            str_sub(Solvent, 1, 1)
            ) %>% str_replace("NANA", "MD")
        )
# Reorganize for saving
sam <- sam %>%
    select(Sample, Type, Primer, Sequencer, Extraction, Solvent, Mock,
        Description = Library, SRA_run)

readr::write_csv(sam, file.path(script_path, "sample_metadata.csv"))
