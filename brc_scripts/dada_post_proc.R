# Combine sequence tables from all labs output by merging step in the DADA2
# pipeline, and finish processing according to the remaining steps shown in the
# DADA2 tutorial.

# Load packages
library(phyloseq)
library(dada2); packageVersion("dada2")

## Paths
# Path for raw sequencing data and pipeline output
data.path <- "~/data/mbqc/blinded_sequence_data/dada_out/"
# Path for silva training data
silva.path <- '~/data/silva/dada2_format'

#### Build combined sequence table and remove chimeras
## Build sequence table with chimeras present
labs <- c('B', 'C', 'E', 'F', 'H', 'J', 'K', 'N')
seqtab.paths <- file.path(data.path, paste0('seqtab_', labs, '.rds'))
seqtabs <- lapply(seqtab.paths, readRDS)
st.all <- do.call(mergeSequenceTables, seqtabs) 
saveRDS(st.all, file.path(data.path, "seqtab_all.rds"))
print(paste(sum(st.all), 'reads across', nrow(st.all), 'samples and', ncol(st.all), 'ASVs'))
## Remove chimeras
st <- removeBimeraDenovo(st.all, multithread=TRUE, verbose=TRUE)
saveRDS(st, file.path(data.path, "seqtab_all_nochim.rds"))
print(paste(sum(st), 'reads across', nrow(st), 'samples and', ncol(st), 'ASVs'))
remove(seqtabs, st.all)
gc()

#### Assign taxonomy
# Assign tax up to Genus
tax <- assignTaxonomy(st, file.path(silva.path,
        "silva_nr_v128_train_set.fa.gz"), multithread=TRUE)
saveRDS(tax, file.path(data.path, "taxonomy.rds"))
# Assign species
spec <- assignSpecies(st, file.path(silva.path,
        "silva_species_assignment_v128.fa.gz"), allowMultiple=TRUE)
saveRDS(spec, file.path(data.path, "species.rds"))
