# Load packages
library(phyloseq)
library(dada2); packageVersion("dada2")

## Paths
# Path for raw sequencing data and pipeline output
data.path <- "~/active_research/metagenomics_calibration/mbqc/data"
# Path for silva training data
silva.path <- '~/data/silva/dada2_format'

#### Get the sample data for incorporation into a phyloseq object
# Will include samples from all labs except HL-A. Samples from other labs not
# being used will be filtered automatically when merging into a phyloseq object 
sd.all <- readRDS(file.path(data.path, 'mbqc_sample_data.rds'))

## Determine which of the original sample data fields to keep
# Helper functions for classifying sample variables
is.drylab.var <- function (field) {
    a = table(sd.all[, c('dry_lab', field)], useNA='no')
    all(rowSums(a>0) <= 1)
}
is.extlab.var <- function (field) {
    a = table(sd.all[, c('extraction_wetlab', field)], useNA='no')
    all(rowSums(a>0) <= 1)
}
is.seqlab.var <- function (field) {
    a = table(sd.all[, c('sequencing_wetlab', field)], useNA='no')
    all(rowSums(a>0) <= 1)
}
# Classify variables as dry lab, sequencing, or extraction. Seq and Ext may
# overlap b/c of the crude classification and fact that there is usually only
# local and central extraction. We also miss some extraction vars like kit
# maker and model
fields <- colnames(sd.all)
drylab.vars <- fields[sapply(fields, is.drylab.var)]
seqlab.vars <- fields[sapply(fields, is.seqlab.var)]
extlab.vars <- fields[sapply(fields, is.extlab.var)]
wetlab.vars <- union(seqlab.vars, extlab.vars)
# Make sure dry and wetlab vars don't overlap
intersect(drylab.vars, wetlab.vars)
# See what vars we didn't classify, and which ones we should get rid of
setdiff(fields, union(drylab.vars, wetlab.vars))
# Some of these are measures of diversity in the final processed community data
diversity.vars <- c("observed_species", "simpson_reciprocal", "chao1", "PD_whole_tree")
# Keep all except the drylab and diversity variables
final.vars <- setdiff(fields, union(drylab.vars, diversity.vars))
length(final.vars) # 49 variables remaining

## Build new sample data table
# Restrict to relevant vars and get rid of duplicate rows
sd1 <- unique(sd.all[,final.vars])
# Check that each Bioinformatics.ID now appears exactly oncd
a <- table(sd1$Bioinformatics.ID, useNA='ifany')
# 36 IDs appear twice
a[a>1]
# I think these are the samples that were extracted in lab A (or assigned to A
# but centrally extracted) and sequenced in both A and E. Let's check
# TODO: Ask MBQC about why these samples have duplicated Bioinf IDS
problem.ids <- names(a[a>1])
subset(sd1, Bioinformatics.ID %in% problem.ids, 
       select=c(extraction_wetlab, sequencing_wetlab, blinded_lab))
# That seems to be the case. So we should be ok if we get rid of all the
# samples sequenced in A
sd2 <- subset_samples(sd1, sequencing_wetlab != 'HL-A')
a2 <- table(sd2$Bioinformatics.ID, useNA='ifany')
all(a2==1) # TRUE
# Pad the IDs to 10 characters and set as sample names
target.length <- 10
sd2$Bioinformatics.ID <- stringr::str_pad(sd2$Bioinformatics.ID,
                                          target.length, side='left',
                                          pad='0')
sample_names(sd2) <- sd2$Bioinformatics.ID
sampledata <- sd2
remove(sd.all, sd1, sd2)
saveRDS(sampledata, file.path(data.path, 'brc_dada_out', "sample_data.rds"))

#### Build sequence table for all specimens
## Build sequence table with chimeras present
labs <- c('B', 'C', 'E', 'F', 'H', 'J', 'K', 'N')
seqtab.paths <- file.path(data.path, 'brc_dada_out', paste0('seqtab_', labs, '.rds'))
seqtabs <- lapply(seqtab.paths, readRDS)
st.all <- do.call(mergeSequenceTables, seqtabs) 
saveRDS(st.all, file.path(data.path, 'brc_dada_out', "seqtab_all.rds"))
print(paste(sum(st.all), 'reads across', nrow(st.all), 'samples and', ncol(st.all), 'ASVs'))
# [1] "77705284 reads across 1796 samples and 196025 ASVs"
## Remove chimeras
st <- removeBimeraDenovo(st.all, multithread=TRUE, verbose=TRUE)
saveRDS(st, file.path(data.path, 'brc_dada_out', "seqtab_all_nochim.rds"))
# st <- readRDS(file.path(data.path, 'brc_dada_out', "seqtab_all_nochim.rds"))
print(paste(sum(st), 'reads across', nrow(st), 'samples and', ncol(st), 'ASVs'))
# [1] "68881850 reads across 1796 samples and 31311 ASVs"
remove(seqtabs, st.all)

#### Filter samples and SVs
# Doing this now will speed up taxonomy assignment, since many ASVs have very
# low prevalence and/or total abundance. So let's just do something somewhat
# arbitrary but probably safe to reduce the number of ASVs
asv.prev <- colSums(st>0)
# Over half of ASVs have prevalance = 1. Filtering to prevalence >= 20 reduces
# the number of SVs from ~30,000 to ~5,000
qplot(seq(30), sapply(seq(30), function(min.prev) sum(asv.prev>=min.prev)), 
     xlab="Minimum Prevalence", ylab="Number of ASVs")
# qplot(seq(20), sapply(seq(20), function(min.prev) sum(st[,asv.prev>=min.prev])), 
#      xlab="Minimum Prevalence", ylab="Total Reads")
# Filtering to prev >= 20 keeps ~99% of reads
sum(st[,asv.prev>=20]) / sum(st) # 0.9892816
st.filt <- st[,asv.prev>20]
saveRDS(st.filt, file.path(data.path, 'brc_dada_out', "seqtab_filt.rds"))

#### Assign taxonomy
system.time(tax <- assignTaxonomy(st.filt, file.path(silva.path,
            "silva_nr_v128_train_set.fa.gz"), multithread=TRUE))
# Took 13 minutes to assign tax for ~5000 ASVs
saveRDS(tax, file.path(data.path, 'brc_dada_out', "taxonomy_filt.rds"))
# Assign species (this just takes one to a few minutes)
spec <- assignSpecies(st.filt, file.path(silva.path,
        "silva_species_assignment_v128.fa.gz"), allowMultiple=TRUE)
saveRDS(spec, file.path(data.path, 'brc_dada_out', "species_filt.rds"))
## Also do for the full set of ASVs
tax <- assignTaxonomy(st, file.path(silva.path,
        "silva_nr_v128_train_set.fa.gz"), multithread=TRUE)
saveRDS(tax, file.path(data.path, 'brc_dada_out', "taxonomy.rds"))
# species assignment has to be run on the cluster: On my pc, fails with "Error:
# cannot allocate vector of size 17.1 Gb". 
spec <- assignSpecies(st, file.path(silva.path,
        "silva_species_assignment_v128.fa.gz"), allowMultiple=TRUE)
saveRDS(spec, file.path(data.path, 'brc_dada_out', "species.rds"))
