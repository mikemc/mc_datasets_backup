# Combine sequence tables from all labs output by merging step in the DADA2
# pipeline, and finish processing according to the remaining steps shown in the
# DADA2 tutorial.

# Load packages
# library(phyloseq)
library(dada2); packageVersion("dada2")

## Paths
# Path for raw sequencing data and pipeline output
data.path <- "~/data/mbqc/blinded_sequence_data/dada_out/"
# Path where this and other BRC DADA2 scripts are
script.path <- "~/metagenomics_calibration/mbqc/dada2_pipeline/brc_scripts"

#### Build combined sequence table and remove chimeras
## Build sequence table with chimeras present
labs <- c('B', 'C', 'E', 'F', 'H', 'J', 'K', 'N')
seqtab.paths <- file.path(data.path, paste0('seqtab_', labs, '.Rds'))
seqtabs <- lapply(seqtab.paths, readRDS)
st.all <- do.call(mergeSequenceTables, seqtabs) 
saveRDS(st.all, file.path(data.path, "seqtab_all.Rds"))
print(paste(sum(st.all), 'reads across', nrow(st.all), 'samples and', ncol(st.all), 'ASVs'))
## Remove chimeras
st <- removeBimeraDenovo(st.all, multithread=TRUE, verbose=TRUE)
saveRDS(st, file.path(data.path, "seqtab_all_nochim.Rds"))
print(paste(sum(st), 'reads across', nrow(st), 'samples and', ncol(st), 'ASVs'))
remove(seqtabs, st.all)
gc()

## Assign taxonomy
source(file.path(script.path, "assign_taxonomy.R"))
