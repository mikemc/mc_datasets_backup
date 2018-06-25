library(tidyverse)
library(rentrez)
library(xml2)

dotenv::load_dot_env("../.env")
script_path <- getwd()
data_path <- file.path(Sys.getenv("DATA_DIR"),
    "brooks2015")

# The Brooks2015 dataset is filed under the BioProject
# https://www.ncbi.nlm.nih.gov/bioproject/267701. We want to get the SRA run
# accessions for all samples in the BioProject. One way to do this is to get
# the full record for all entries in the SRA database associated with the
# BioProject, which is stored in xml format.
sra_id <- entrez_link(dbfrom = "bioproject", id = "267701", db="sra") %>% 
    {.$links$bioproject_sra}
xml <- entrez_fetch(db = "sra", id = sra_id, rettype = "xml") %>%
    read_xml
# To view, run `write_lines(xml, "/tmp/xml.txt")` and open in a text editor.
# Doing so, we can see that the sample names given by the authors are stored in
# an element named "LIBRARY_NAME" (the only element with such a name), the SRA
# run (needed for downloading the reads) is stored in an element matching
# ".//RUN/IDENTIFIERS/PRIMARY_ID", and the rest of the info we need for that
# sample is contained in an element matching ".//DESIGN_DESCRIPTION".
tb <- tibble(
    Sample = xml %>%
        xml_find_all(".//LIBRARY_NAME") %>%
        xml_text,
    Run = xml %>%
        xml_find_all(".//RUN/IDENTIFIERS/PRIMARY_ID") %>%
        xml_text,
    Description = xml %>% 
        xml_find_all(".//DESIGN_DESCRIPTION") %>%
        xml_text %>%
        str_sub(2, -2),
    )
# Note that the `Sample` names have the format of `TRUTHa_b_c_d` where `a` is
# the plate, `b` is the barcode ID, and I'm not yet sure what `c` and `d` are;
# perhaps these relate to the experimental design.

# We can now parse the Description for the needed sample metadata:
template <- paste0(
    "Sequencing results for a mock community created by mixing",
    " equal amounts of (Cells|DNA|PCR Product) of",
    " \\{(.+) \\}", 
    " with barcode ID: ([0-9]+), barcode: ([A-Z]*), on plate: ([1-6])"
    )
desc_vars <- tb$Description %>% str_match(template) %>%
    as_tibble %>%
    select(-V1, Type = V2, Species = V3, BarcodeID = V4, Barcode = V5, 
        Plate = V6)
# The `Species` field has the form of a list separated by spaces; 
# e.g., "L. iners G. vaginalis P. bivia". We can use this form to calculate the
# the number of species in each mixed culture
desc_vars <- desc_vars %>%
    mutate(Num_species = str_count(Species, "[A-Z]\\. [a-z]+"))
# We'll save these variables from the description, along with the original
# description and the SRA run, in the same path as this script.
tb0 <- bind_cols(tb, desc_vars) %>%
    select(Sample, Type, Species, Num_species, 
        Plate, BarcodeID, Barcode,
        Description, 
        SRA_run = Run)
readr::write_csv(tb0, file.path(script_path, "sample_metadata.csv"))
