# Load packages
library(phyloseq)
library(dada2); packageVersion("dada2")

## Paths
# Path for raw sequencing data and pipeline output
data.path <- "~/data/mbqc/blinded_sequence_data/dada_out/"
# Path for silva training data
silva.path <- '~/data/silva/dada2_format'

## Load sequence table
st <- readRDS(file.path(data.path, "seqtab_all_nochim.rds"))

## Assign species
# Need at least 17.1 Gb of ram
spec <- assignSpecies(st, file.path(silva.path,
        "silva_species_assignment_v128.fa.gz"), allowMultiple=TRUE)
saveRDS(spec, file.path(data.path, "species.rds"))

