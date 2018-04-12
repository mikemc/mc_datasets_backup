#!/usr/bin Rscript
# Full dada2 pipeline for paired-end labs. Argument must be the name (letter
# only) of a sequencing lab besides A and I 

## Initialize 
# Lab should be supplied by the first (and only) command line argument
args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
  stop("Sequencing lab must be supplied\n", call.=FALSE)
} else if (length(args) > 1) {
  stop("Too many args\n", call.=FALSE)
}
lab <- args[1] 
seqlabs <- c('B', 'C', 'D', 'E', 'F', 'H', 'J', 'K', 'L', 'M', 'N')
if (!(lab %in% seqlabs)) {
  stop("Invalid sequencing lab\n", call.=FALSE)
}

print(paste0('Run Dada2 pipeline for HL-', lab))

# Load packages
library(phyloseq)
library(dada2); packageVersion("dada2")

# Path for raw sequencing data and pipeline output
data.path <- "~/data/mbqc/blinded_sequence_data"
# Path for MC git
mc.path <- '~/metagenomics_calibration'

# Primer seqs
F515 <- "GTGCCAGCMGCCGCGGTAA"
R806 <- "GGACTACHVGGGTWTCTAAT"
F318 <- "ACTCCTACGGGAGGCAGCAG"
F515_HLC <- "GTGTGCCAGCMGCCGCGGTAA"
            
## Get the sample names (Bioinformatics.ID) into a dataframe
sd.all <- readRDS(file.path(mc.path, 'mbqc/data/', 'mbqc_sample_data.Rds'))
df <- subset_samples(sd.all, sequencing_wetlab != 'HL-A')
df <- unique(df[, c('Bioinformatics.ID', 'sequencing_wetlab')])
target.length = 10
df$Bioinformatics.ID <- stringr::str_pad(df$Bioinformatics.ID,
                                         target.length, side='left',
                                         pad='0')
sample_names(df) <- df$Bioinformatics.ID
remove(sd.all)

## Trim and truncation parameters
trim.params <- data.frame(matrix(0, ncol = 4, nrow = length(seqlabs)))
rownames(trim.params) <- seqlabs
colnames(trim.params) <- c('trimF', 'trimR', 'truncF', 'truncR')
trim.params['B', c('truncF', 'truncR')] <- c(230, 210)
trim.params['C', c('trimF', 'trimR')] <- nchar(c(F515_HLC, R806))
trim.params['C', c('truncF', 'truncR')] <- c(240, 245)
trim.params['D', c('truncF', 'truncR')] <- c(245, 245)
trim.params['E', c('truncF', 'truncR')] <- c(148, 148)
trim.params['F', c('trimF', 'trimR', 'truncF', 'truncR')] <- c(1, 1, 148, 148)
trim.params['H', c('truncF', 'truncR')] <- c(245, 220)
trim.params['J', c('truncF', 'truncR')] <- c(245, 228)
trim.params['K', c('truncF', 'truncR')] <- c(148, 148)
trim.params['L', c('truncF', 'truncR')] <- c(127, 148)
trim.params['M', c('truncF', 'truncR')] <- c(160, 120)
trim.params['N', c('trimF', 'trimR', 'truncF', 'truncR')] <- c(1, 0, 225, 222)
# trim.params['I', 'truncF'] = 245

## Set filenames
bio.ids <- df$Bioinformatics.ID[df$sequencing_wetlab==paste0('HL-', lab)]
# Raw reads
fnFs <- file.path(data.path, 'raw', paste0('HL-', lab), paste0(bio.ids, '_R1.fastq.gz'))
fnRs <- file.path(data.path, 'raw', paste0('HL-', lab), paste0(bio.ids, '_R2.fastq.gz'))
names(fnFs) <- names(fnRs) <- bio.ids
# Trimmed and filtered reads
filtFs <- file.path(data.path, 'filtered', paste0(names(fnFs), "_R1_filt.fastq.gz"))
filtRs <- file.path(data.path, 'filtered', paste0(names(fnRs), "_R2_filt.fastq.gz"))
names(filtFs) <- names(filtRs) <- bio.ids

## Trim and filter
# Need to use matchIDs for lab N
cat('Trimming and filtering', length(fnFs), 'samples\n')
out <- filterAndTrim(fnFs, filtFs,
                     fnRs, filtRs,
                     trimLeft=trim.params[lab, c('trimF', 'trimR')],
                     truncLen=trim.params[lab, c('truncF', 'truncR')],
                     matchIDs=ifelse(lab=='N', TRUE, FALSE),
                     maxN=0, maxEE=c(2,2), truncQ=2,
                     rm.phix=TRUE, compress=TRUE, multithread=TRUE)
print(head(out))
# Future steps only use samples where >0 reads passed the filter
bad.ids <- bio.ids[out[,'reads.out'] == 0]
ok.ids <- bio.ids[out[,'reads.out'] > 0]
if (length(bad.ids) > 0) {
    cat('Omitting', bad.ids, '\n')
}
## Derep
print('Dereplicating')
derepFs <- derepFastq(filtFs[ok.ids], verbose=TRUE)
derepRs <- derepFastq(filtRs[ok.ids], verbose=TRUE)
names(derepFs) <- names(derepRs) <- ok.ids
## Learn errors
print('Learning errors')
set.seed(1)
print(system.time(errF <- learnErrors(derepFs, multithread=TRUE, randomize=TRUE)))
set.seed(1)
print(system.time(errR <- learnErrors(derepRs, multithread=TRUE, randomize=TRUE)))
## Save the error profiles
print(paste('Saving error profiles to', file.path(data.path, 'dada_out')))
saveRDS(errF, file.path(data.path, 'dada_out', paste0('errF_', lab, '.Rds')))
saveRDS(errR, file.path(data.path, 'dada_out', paste0('errR_', lab, '.Rds')))
## Run dada2
print('Running dada')
print(system.time(dadaFs <- dada(derepFs, err=errF, pool=TRUE, multithread=TRUE)))
print(system.time(dadaRs <- dada(derepRs, err=errR, pool=TRUE, multithread=TRUE)))
print(dadaFs[[1]])
## Merge the denoised forward and reverse reads:
print('Merging')
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
print(head(mergers[[1]]))
## Construct sequence table
print('Constructing sequence table')
seqtab <- makeSequenceTable(mergers)
print(dim(seqtab))
print(table(nchar(getSequences(seqtab))))
## Remove chimeras
print('Removing chimeras')
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus",
                                    multithread=TRUE, verbose=TRUE)
print(dim(seqtab.nochim))
print(sum(seqtab.nochim)/sum(seqtab))
## Track reads through pipeline
getN <- function(x) sum(getUniques(x))
out.ok <- out[out[,'reads.out'] > 0,]
track <- cbind(out.ok, sapply(dadaFs, getN), sapply(mergers, getN),
               rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled",
                     "nonchim")
rownames(track) <- ok.ids
# Add on the samples that failed to pass filter
out.bad <- out[out[,'reads.out'] == 0,]
track.bad <- cbind(out.bad, 0, 0, 0, 0)
rownames(track.bad) <- bad.ids
track <- rbind(track, track.bad)
print('Track reads through pipeline (6 samples shown)')
print(head(track))
## Save objects
print(paste('Saving objects to', file.path(data.path, 'dada_out')))
saveRDS(seqtab, file.path(data.path, 'dada_out', 
                          paste0('seqtab_', lab, '.Rds')))
saveRDS(seqtab.nochim, file.path(data.path, 'dada_out',
                                 paste0('seqtab_nochim_', lab, '.Rds')))
saveRDS(track, file.path(data.path, 'dada_out', paste0('track_', lab, '.Rds')))
