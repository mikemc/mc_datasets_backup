# Assign taxonomy down to species level.

# Load packages
library(dada2); packageVersion("dada2")

## Paths
# Path for raw sequencing data and pipeline output
data.path <- "~/data/mbqc/blinded_sequence_data/dada_out/"
# Path for silva training data
silva.path <- '~/data/silva/dada2_format'

## Taxonomy databases
tax.db <- file.path(silva.path, "silva_nr_v132_train_set.fa.gz")
species.db <- file.path(silva.path, "silva_species_assignment_v132.fa.gz")

## Load sequence table
st <- readRDS(file.path(data.path, "seqtab_all_nochim.Rds"))

## Assign taxonomy
# Assign tax up to Genus
taxa <- assignTaxonomy(st, tax.db, tryRC = TRUE, multithread=TRUE)
# Add species
taxa <- addSpecies(taxa, species.db, tryRC = TRUE, allowMultiple=TRUE)
saveRDS(taxa, file.path(data.path, 
        paste0(Sys.Date(), "taxonomy.Rds", sep="_")))
# Also save species assignment for later inspection
spec <- assignSpecies(st, species.db, tryRC = TRUE, allowMultiple=TRUE)
saveRDS(spec, file.path(data.path, 
        paste0(Sys.Date(), "species.Rds", sep="_")))
