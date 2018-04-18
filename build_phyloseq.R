# Load packages
library(phyloseq)
library(stringr)

#### Build a core-ASV phyloseq object
data_path <- "~/active_research/metagenomics_calibration/mbqc/data"

# Load data and build new phyloseq object
sampledata <- readRDS(file.path(data_path, 'brc_dada_out', "sample_data.Rds"))
seqtab <- readRDS(file.path(data_path, 'brc_dada_out', "seqtab_all_nochim.Rds"))
tax <- readRDS(file.path(data_path, 'brc_dada_out', "2018-04-12_taxonomy.Rds"))

ps <- phyloseq(otu_table(seqtab, taxa_are_rows = FALSE),
    sample_data(sampledata), tax_table(tax))
sequences <- Biostrings::DNAStringSet(taxa_names(ps))
names(sequences) <- taxa_names(ps)
ps <- merge_phyloseq(ps, sequences)
remove(sampledata, seqtab, tax, sequences)
gc()

## Filter to just the core SVs seen in every sequencing lab
# First, get a 0-1 array specifiying if the SV was found in a seqlab
st <- merge_samples(otu_table(ps), sample_data(ps)$sequencing_wetlab)
st <- transform_sample_counts(st, function (x) ifelse(x > 0, 1, 0))
# Prune to just taxa that show up in all labs
keep.taxa <- apply(st, 2, all)
ps0 <- prune_taxa(keep.taxa, ps)
# The number of taxa is reduced to 843 but 94.6% of reads are kept
sum(otu_table(ps0)) / sum(otu_table(ps))

# Compare the tax table before and after: One Archaeal SV (of 8) is maintained
# and all Eukaryotes are dropped
table(tax_table(ps)[,'Kingdom'])
table(tax_table(ps0)[,'Kingdom'])
# All remaining taxa have a characterized Phylum, which is nice
all(with(data.frame(tax_table(ps0)),
        !is.na(Phylum) & !Phylum %in% c("", "uncharacterized")))

## Sample data processing
ps <- ps0
remove(ps0)
# Add prefix to sample names to avoid phyloseq bug
sample_names(ps) <- paste0('s', sample_data(ps)$Bioinformatics.ID)
## Tweak sample data
# Remap extraction_wetlab=NA to 'central' for pre-extracted samples
sam <- sample_data(ps)
sam$extraction_wetlab[sam$pre.extracted_bool=="TRUE"] <- "Central"
# Add the extraction + sequencing pair as a field. Consists of a two-letter
# string, with 'pJ' for example denoting pre-extraction and sequencing in lab J.
# Negative controls are set to NA
sam$lab_pair <- with(sam, paste0(
        ifelse(pre.extracted_bool, 'p', str_sub(extraction_wetlab, -1)),
        str_sub(sequencing_wetlab, -1)))
sam$lab_pair[is.na(sam$extraction_wetlab)] <- NA
# Simplified versions of the lab names to use when plotting
sam$extr <- str_sub(sam$lab_pair, 1, 1)
sam$seq <- str_sub(sam$lab_pair, -1)
sample_data(ps) <- sam

## Save the filtered ps object, as this is what we will work with for now
saveRDS(ps, file.path(data_path, 'mbqc_core_phyloseq.Rds'))
