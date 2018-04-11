library(stringr)
# library(phyloseq); 
library(dada2);packageVersion("dada2")
library(Biostrings); library(ShortRead)
library(ggplot2)

df = readRDS(file.path('~/active_research/metagenomics_calibration/mbqc/data',
             'blinded_sequence_table.rds'))
path = "~/data/mbqc/blinded_sequence_data"
setwd(path)
            
# M matches A or C
# H matches not G
# V matches not T
# W matches A or T
F515 <- "GTGCCAGCMGCCGCGGTAA"
R806 <- "GGACTACHVGGGTWTCTAAT"
F318 <- "ACTCCTACGGGAGGCAGCAG"
# HL-C used an F515 primer with an extra 2bp inserted in after the first 2bp of
# the F515 primer sequence quoted in the paper (see
# mbqc_handling_protocols.xlsx)
F515_HLC <- "GTGTGCCAGCMGCCGCGGTAA"

# All labs have paired end reads except for HL-I
fns.all = list.files(path, recursive=TRUE)
fnIs = grep(pattern="raw/HL-I/.*_R1.fastq.gz$", fns.all, value=TRUE)
fnFs = grep(pattern="raw/HL-[^I]/.*_R1.fastq.gz$", fns.all, value=TRUE)
fnRs = grep(pattern="raw/HL-[^I]/.*_R2.fastq.gz$", fns.all, value=TRUE)
names(fnIs) <- sapply(strsplit(fnIs, "[_/]"), "[", 3)
names(fnFs) <- sapply(strsplit(fnFs, "[_/]"), "[", 3)
names(fnRs) <- sapply(strsplit(fnRs, "[_/]"), "[", 3)
identical(names(fnRs), names(fnFs))

#### Choose parameters 

# Load a few samples of F and R sequences for each lab to check quality
# profiles. I'll choose three from each lab except for the one lab (M) for
# which only there are only two.
set.seed(1)
split.samples = split(df$Bioinformatics.ID, df$sequencing_wetlab)
names(split.samples) = str_sub(names(split.samples), -1)
seqlabs = names(split.samples)
num.samples = sapply(split.samples,length)
test.samples = mapply(sample, split.samples, sapply(num.samples, min, 3))

trim.params =  data.frame(matrix(0, ncol = 4, nrow = length(seqlabs)))
rownames(trim.params) = seqlabs
colnames(trim.params) = c('trimF', 'trimR', 'truncF', 'truncR')

## read lengths example
# fq = readFastq(fnFs[1])
# table(width(fq))

# Function for getting the read length distribution from a list of fastq files
read_lengths = function(fns) {
    lapply(lapply(fns, readFastq), function (x) width(x))
}
read_lengths_summary = function(fns) {
    lapply(read_lengths(fns), summary)
}
read_lengths_table = function(fns) {
    lapply(read_lengths(fns), table)
}

## Check quality profiles and read lengths; pick params
# I will aim for trimming at least 2 bp from the ends of the shortest reads.
# Quality profiles for all labs
#
# qpFs = lapply(test.samples[-7], function(x) plotQualityProfile(fnFs[x]))
# qpRs = lapply(test.samples[-7], function(x) plotQualityProfile(fnRs[x]))
# qpIs = plotQualityProfile(fnIs[test.samples$I])
# saveRDS(qpFs, file.path(path, 'qc', paste0('quality_profiles_Fs', '.rds')))
# saveRDS(qpRs, file.path(path, 'qc', paste0('quality_profiles_Rs', '.rds')))
# saveRDS(qpIs, file.path(path, 'qc', paste0('quality_profiles_labI', '.rds')))
#
# Load previously saved quality profiles
qpFs = readRDS(file.path(path, 'qc', paste0('quality_profiles_Fs', '.rds')))
qpRs = readRDS(file.path(path, 'qc', paste0('quality_profiles_Rs', '.rds')))
qpIs = readRDS(file.path(path, 'qc', paste0('quality_profiles_labI', '.rds')))

## Check read lengths
# rlFs = lapply(test.samples[-7], function(x) read_lengths(fnFs[x])) 
# rlRs = lapply(test.samples[-7], function(x) read_lengths(fnRs[x])) 
# rlIs = read_lengths(fnIs[test.samples$I])
# saveRDS(rlFs, file.path(path, 'qc', paste0('read_lengths_Fs', '.rds')))
# saveRDS(rlRs, file.path(path, 'qc', paste0('read_lengths_Rs', '.rds')))
# saveRDS(rlIs, file.path(path, 'qc', paste0('read_lengths_labI', '.rds')))
rlFs = readRDS(file.path(path, 'qc', paste0('read_lengths_Fs', '.rds')))
rlRs = readRDS(file.path(path, 'qc', paste0('read_lengths_Rs', '.rds')))
rlIs = readRDS(file.path(path, 'qc', paste0('read_lengths_labI', '.rds')))

## HL-B
# Forward
read_lengths(fnFs[test.samples$B]) # Min 246; vast majority 251 or 250
qpFs$B + geom_vline(xintercept=230) # Moderate quality drop after 230
# Reverse
read_lengths(fnRs[test.samples$B]) # Min 246; vast majority 251 or 250
qpRs$B + geom_vline(xintercept=210) # Spike around 160, drop at 210
# Set params
trim.params['B', c('truncF', 'truncR')] = c(230, 210)

## HL-C
# Forward
read_lengths(fnFs[test.samples$C]) # All 243
qpFs$C + geom_vline(xintercept=240)
# Reverse
read_lengths(fnRs[test.samples$C]) # All 251
qpRs$C + geom_vline(xintercept=245)
# Set params
trim.params['C', c('trimF', 'trimR')] = nchar(c(F515_HLC, R806))
trim.params['C', c('truncF', 'truncR')] = c(240, 245)
# Reads will have length truncX - trimX

## HL-D
# Forward
read_lengths(fnFs[test.samples$D]) # Wide range, Most are 251. Suggests trunc 246
qpFs$D + geom_vline(xintercept=245) # Stays high, so let's use 245
# Reverse
read_lengths(fnRs[test.samples$D]) # Wide range, Most are 251. Suggests trunc 245
qpRs$D + geom_vline(xintercept=245) # Mod drop at 245
# Set params
trim.params['D', c('truncF', 'truncR')] = c(245, 245)

## HL-E
# Forward
read_lengths(fnFs[test.samples$E]) # vast majority 151bp; suggests 149
qpFs$E + geom_vline(xintercept=148) # trunc 148 looks good
# Reverse
read_lengths(fnRs[test.samples$E]) # vast majority 148 to 151 bp
qpRs$E + geom_vline(xintercept=148) # Q stays good; truc at 148
# Set params
trim.params['E', c('truncF', 'truncR')] = c(148, 148)

## HL-F
# Forward
read_lengths_table(fnFs[test.samples$F]) # F reads all 152 bp
qpFs$F + geom_vline(xintercept=148)
# R
read_lengths_table(fnRs[test.samples$F]) # R reads all 152 bp 
qpRs$F + geom_vline(xintercept=148)
# Set params
trim.params['F', c('truncF', 'truncR')] = c(148, 148)
# 2018-01-19: Results from running the dada2 pipeline suggest we need to trim
# the first bp from the forward and reverse reads, which perhaps resulted from
# incompeltely trimmed primer sequences. On second look at the QPs, we can see
# that the first bp in Fs and Rs has signiticantly higher Q score than the
# following bps, which might be an indicator of such a problem.

## HL-H
# Forward
read_lengths(fnFs[test.samples$H]) # all 251
qpFs$H + geom_vline(xintercept=245) # slow decay; stays decent
# Reverse
read_lengths(fnRs[test.samples$H]) # all 251
qpRs$H + geom_vline(xintercept=220) # Mod decay, avg <30 at 220
# Set params
trim.params['H', c('truncF', 'truncR')] = c(245, 220)

## HL-J
# Forward
read_lengths(fnFs[test.samples$J]) # all 250
qpFs$J + geom_vline(xintercept=245) # Small drop ~235, but ok
# Reverse
read_lengths(fnRs[test.samples$J]) # all 250
qpRs$J + geom_vline(xintercept=228) # Mod drop at ~230
# Set params
trim.params['J', c('truncF', 'truncR')] = c(245, 228)

## HL-K
# Forward
read_lengths(fnFs[test.samples$K]) # Vast majority 151
qpFs$K + geom_vline(xintercept=148) # stays good
# Reverse
read_lengths(fnRs[test.samples$K]) # Vast majority 150 or 151
qpRs$K + geom_vline(xintercept=148) # stays good
# Set params
trim.params['K', c('truncF', 'truncR')] = c(148, 148)

## HL-L
# Samples have lots of reads, but terrible quality. Only 150bp paired end, so
# getting 20bp overlap is going to be a problem if we truncate heavily.
# Forward reads are worse, so will trim as much as possible from them
# Lengths:
read_lengths(fnFs[test.samples$L]) # all 150
read_lengths(fnRs[test.samples$L]) # all 150
# Forward: Big drop at 81 but also terrible after 94
# Trim as much as possible to get 20 bp overlap (assuming max 255 length reads)
qpFs$L + geom_vline(xintercept=127)
# Reverse reads: Mode drops to Q~2 at ~100bp, but rebounds
qpRs$L + geom_vline(xintercept=148)
# Set params
trim.params['L', c('truncF', 'truncR')] = c(127, 148)

## HL-M
# Forward
read_lengths(fnFs[test.samples$M]) # all 175
# Ok, with slow drop at end. Much better than reverse, so keep as much as we
# can
qpFs$M + geom_vline(xintercept=173)
# Reverse
read_lengths(fnRs[test.samples$M]) # all 175
# Generally terrible quality throughout, with drop at end
# 25-percentile drops to Q2 at ~150
# Mean drops to ~20 at 135
qpRs$M + geom_vline(xintercept=133)
# Set params
trim.params['M', c('truncF', 'truncR')] = c(173, 133)

## HL-N
# Quality decent, but a drop off we can cut since long reads
# Forward
read_lengths_table(fnFs[test.samples$N]) # vast maj 251-253
qpFs$N + geom_vline(xintercept=225) # avg ~Q30 at ~225
# Reverse
read_lengths_table(fnRs[test.samples$N]) # vast maj 247-249
qpRs$N + geom_vline(xintercept=222) # avg ~Q30 at ~222
# Set params
trim.params['N', c('truncF', 'truncR')] = c(225, 222)
# 2018-01-19: Results of the dada2 pipeline suggest that we need to trim the
# first bp from the forward reads. Unlike for lab F, the QPs in this case don't
# seem to indicate anything special about the first bp.

## HL-I, single-end reads
# vast maj 252-253, with one ~315
# Not sure if I should truncate these reads. There is a big drop to Q=2 in many
# reads near 248, but I can just let dada truncate at these
read_lengths(fnIs[test.samples$I])
qpIs
# Try not trimming?
# trim.params['I', 'truncF'] = 

# Lines copied from above to make sure all params set
trim.params['B', c('truncF', 'truncR')] = c(230, 210)
trim.params['C', c('trimF', 'trimR')] = nchar(c(F515_HLC, R806))
trim.params['C', c('truncF', 'truncR')] = c(240, 245)
trim.params['D', c('truncF', 'truncR')] = c(245, 245)
trim.params['E', c('truncF', 'truncR')] = c(148, 148)
trim.params['F', c('truncF', 'truncR')] = c(148, 148)
trim.params['H', c('truncF', 'truncR')] = c(245, 220)
trim.params['J', c('truncF', 'truncR')] = c(245, 228)
trim.params['K', c('truncF', 'truncR')] = c(148, 148)
trim.params['L', c('truncF', 'truncR')] = c(127, 148)
# trim.params['M', c('truncF', 'truncR')] = c(173, 133)
# trim.params['M', c('truncF', 'truncR')] = c(160, 130)
trim.params['M', c('truncF', 'truncR')] = c(160, 120)
trim.params['N', c('truncF', 'truncR')] = c(225, 222)
trim.params['I', 'truncF'] = 245
trim.params
## Check overlap lengths
# Want at least 20 when true sequence length is 255
with(trim.params, truncF + truncR - trimF - trimR - 255)

# Save trim.params
saveRDS(trim.params, file.path(path, 'dada_out', paste0('trim_params', '.rds')))


#### Trim and filter; learn errors; run dada2

filt_path = file.path(path, "filtered")
filtFs = file.path(filt_path, paste0(names(fnFs), "_R1_filt.fastq.gz"))
names(filtFs) = names(fnFs)
filtRs = file.path(filt_path, paste0(names(fnRs), "_R2_filt.fastq.gz"))
names(filtRs) = names(fnRs)
filtIs = file.path(filt_path, paste0(names(fnIs), "_R1_filt.fastq.gz"))
names(filtIs) = names(fnIs)

# Make sure folders exist: filtered, dada_output

# run.times =  data.frame(matrix(0, ncol = 4, nrow = length(seqlabs)))
# rownames(run.times) = seqlabs
# colnames(run.times) = c('errF', 'errR', 'dadaF', 'dadaR')

## Labs with paired-end reads
# for (lab in seqlabs[seqlabs!='I']) {
for (lab in seqlabs[-c(7)]) {
    # Delete old objects to make sure they don't get incorrectly saved in case
    # of an error, and clear memory
    remove(out, derepFs, derepRs, errF, errR, dadaFs, dadaRs, mergers, seqtab,
           seqtab.nochim, track)
    gc()
    # Get the samples we will be processing
    sample.names = split.samples[[lab]]
    print(paste('Starting HL', lab, 'with', length(sample.names), 'samples'))
    ## Trim and filter
    print('Trimming and filtering')
    out = filterAndTrim(fnFs[sample.names], filtFs[sample.names],
                        fnRs[sample.names], filtRs[sample.names],
                        trimLeft=trim.params[lab, c('trimF', 'trimR')],
                        truncLen=trim.params[lab, c('truncF', 'truncR')],
                        maxN=0, maxEE=c(2,2), truncQ=2,
                        rm.phix=TRUE, compress=TRUE, multithread=TRUE)
    print(head(out))
    ## Derep
    print('Dereplicating')
    derepFs <- derepFastq(filtFs[sample.names], verbose=TRUE)
    derepRs <- derepFastq(filtRs[sample.names], verbose=TRUE)
    # Name the derep-class objects by the sample names
    names(derepFs) <- sample.names
    names(derepRs) <- sample.names
    ## Learn errors
    print('Learning errors')
    print(system.time(errF <- learnErrors(derepFs, multithread=TRUE)))
    print(system.time(errR <- learnErrors(derepRs, multithread=TRUE)))
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
    track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN),
                   rowSums(seqtab), rowSums(seqtab.nochim))
    colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled",
                         "nonchim")
    rownames(track) <- sample.names
    print('Track reads through pipeline (6 samples shown)')
    print(head(track))
    ## Save objects
    print(paste('Saving objects to', file.path(path, 'dada_out')))
    saveRDS(errF, file.path(path, 'dada_out', paste0('errF_', lab, '.rds')))
    saveRDS(errR, file.path(path, 'dada_out', paste0('errR_', lab, '.rds')))
    saveRDS(seqtab, file.path(path, 'dada_out', 
                              paste0('seqtab_', lab, '.rds')))
    saveRDS(seqtab.nochim, file.path(path, 'dada_out',
                                     paste0('seqtab_nochim_', lab, '.rds')))
    saveRDS(track, file.path(path, 'dada_out', paste0('track_', lab, '.rds')))
}

## HL-B run, ~1million reads used
# errF: 7403.465    3.044 2066.953 user/system/elapsed
# errR: 6779.642    1.930 1894.020
# dadaFs: 6468.391    1.731 1807.510
# dadaRs: 5803.288    1.326 1626.241
## HL-C run, only ~380,000 reads
# errF: 1433.130    0.891  403.826
# errR: 2485.705    1.306  691.231
## HL-E run, 
# errF: 3424.731    2.039 1053.212
# errR: 11435.859     2.719  3197.149
# 2151.107    0.650  679.815
# Seems like a lot of chimeric reads
## HL-F
# Total reads used:  314019 
# 378.272   0.783 116.992
# 1064.410    0.784  302.560
# 226.748   0.110  72.569
# 668.331   0.106 191.691
# Some samples had a very high chimera fraction
# Should check if any association of G versus F
## HL-H
# 453068 reads
# 2201.154    2.169  640.072
# 1641.277    1.226  471.087
# Seems to have gone well
## HL-J
# 639046 reads
# 6732.671    4.379 1870.422
# 5176.936    3.802 1462.322
# 5836.274    1.006 1607.726
# 5215.681    2.090 1487.472
# high chimera fraction
## HL-K
# 375106 reads
# 233.277   0.346  76.719
# 326.343   0.300  98.818
# 91.202   0.064  30.881
# 135.041   0.053  41.517
# Good quality? Very high % of reads made it through
## HL-L
# Total reads used:  287509 
# 332.975   0.436  98.966 
# 892.226   0.473 249.778
# 100.221   0.057  30.109
# 326.306   0.116  91.762
# High chimera fraction
# Large majority of reads filtered; should increase EE
## HL-M
# 200963 reads
# 607.396   0.310 181.550
# 748.953   0.220 235.755
# Huge amount filtered during merging; perhaps reads not matched?
## HL-N
# Total reads used:  101225
# 293.714   0.293  84.542
# 399.820   0.320 114.287
# 157.988   0.043  44.592
# 160.413   0.040  45.429
# Table looks ok

# Done: all except D

## Lab N note
# The paired files had the same lengths for this specimen, but not for all
# specimens, and so mergeIDs=TRUE is needed on the cluster run. It may be worth
# looking into the files that don't match up

## HL-D
# Running filterAndTrim with default, matchIDs=FALSE, fails b/c of unequal
# numbers of reads in the R1 and R2 files. Running with matchIDs=TRUE fixes the
# error, but filters a majority of reads (my guess is b/c of id matching, as
# opposed to quality and length issues)
lab='D'
sample.names = split.samples[[lab]]
out = filterAndTrim(fnFs[sample.names], filtFs[sample.names],
                    fnRs[sample.names], filtRs[sample.names],
                    trimLeft=trim.params[lab, c('trimF', 'trimR')],
                    truncLen=trim.params[lab, c('truncF', 'truncR')],
                    maxN=0, maxEE=c(2,2), truncQ=2,
                    matchIDs=TRUE,
                    rm.phix=TRUE, compress=TRUE, multithread=TRUE)
# None will merge :(
# Troubleshooting: 
# mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, returnRejects=TRUE)
# mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, justConcatenate=TRUE)
errF = readRDS(file.path(path, 'dada_out', paste0('errF_', lab, '.rds')))
errR = readRDS(file.path(path, 'dada_out', paste0('errR_', lab, '.rds')))
## Track reads through pipeline---nothing after merging
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), 0, 0)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled",
                     "nonchim")
rownames(track) <- sample.names
saveRDS(track, file.path(path, 'dada_out', paste0('track_', lab, '.rds')))

## HL-M
# trim.params['M', c('truncF', 'truncR')] = c(160, 130)
out = filterAndTrim(fnFs[sample.names], filtFs[sample.names],
                    fnRs[sample.names], filtRs[sample.names],
                    trimLeft=trim.params[lab, c('trimF', 'trimR')],
                    truncLen=trim.params[lab, c('truncF', 'truncR')],
                    maxN=0, maxEE=c(10,10), truncQ=2,
                    rm.phix=TRUE, compress=TRUE, multithread=TRUE)


#### HL-I: single-end reads

# For choosing error tolerance, can use vsearch command: 
# vsearch --fastq_eestats2 ~/data/mbqc/blinded_sequence_data/raw/HL-I/4375574631_R1.fastq.gz --ee_cutoffs 1,2,4,6,8 --length_cutoffs 240,250,5 --output -
# Suggests would keep ~90% at maxEE=5 and 245bp but doesn't account for the
# truncation at Q=2.
# Increasing maxEE from 2 to 5 seems to make almost no difference b/c of truncQ = 2

lab = 'I'
# Free memory
gc()
# Get the samples we will be processing
sample.names = split.samples[[lab]]
print(paste('Starting HL', lab, 'with', length(sample.names), 'samples'))
## Trim and filter
print('Trimming and filtering')
out = filterAndTrim(fnIs[sample.names], filtIs[sample.names],
                    trimLeft=0,
                    truncLen=trim.params['I', 'truncF'],
                    maxN=0, maxEE=2, truncQ=2,
                    rm.phix=TRUE, compress=TRUE, multithread=TRUE)
# Truncation length 247 was chosen after seeing that truncLen=0 and truncQ=2
# led to a mode of read lengths 247. Makes sense, since median Q drops to 2 at
# 248. BJC suggested better to truncate at 245, to stay enough away from the
# crappy Q=2 bases, so reran with that on 11/03. Looking at vsearch profile on
# one sample shows that 245 is the point where the median rises to above 30.
print(head(out))
## Derep
print('Dereplicating')
derepIs <- derepFastq(filtIs[sample.names], verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepIs) <- sample.names
## Learn errors
print('Learning errors')
print(system.time(errI <- learnErrors(derepIs, multithread=TRUE)))
## Run dada2
print('Running dada')
print(system.time(dadaIs <- dada(derepIs, err=errI, pool=TRUE, multithread=TRUE)))
print(dadaIs[[1]])
## Merge the denoised forward and reverse reads:
print('Single-end reads; skipping merging')
## Construct sequence table
print('Constructing sequence table')
seqtab <- makeSequenceTable(dadaIs)
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
track <- cbind(out, sapply(dadaIs, getN), rowSums(seqtab),
               rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "tabled", "nonchim")
rownames(track) <- sample.names
print('Track reads through pipeline (6 samples shown)')
print(head(track))
## Save objects
print(paste('Saving objects to', file.path(path, 'dada_out')))
saveRDS(errI, file.path(path, 'dada_out', paste0('err_', lab, '.rds')))
saveRDS(seqtab, file.path(path, 'dada_out', 
                          paste0('seqtab_', lab, '.rds')))
saveRDS(seqtab.nochim, file.path(path, 'dada_out',
                                 paste0('seqtab_nochim_', lab, '.rds')))
saveRDS(track, file.path(path, 'dada_out', paste0('track_', lab, '.rds')))


#### Inspect error profiles
errFs = lapply(file.path(path, 'dada_out', paste0('errF_', seqlabs[-7], '.rds')), readRDS)
errRs = lapply(file.path(path, 'dada_out', paste0('errF_', seqlabs[-7], '.rds')), readRDS)
errI = readRDS(file.path(path, 'dada_out', 'err_I.rds'))
names(errFs) = names(errRs) = seqlabs[-7]

for (lab in seqlabs[-7]) {
    print(paste('HL', lab, ': Forward'))
    print(plotErrors(errFs[[lab]]))
    readline(prompt="Press [enter] to continue")
    print(paste('HL', lab, ': Reverse'))
    print(plotErrors(errRs[[lab]]))
    readline(prompt="Press [enter] to continue")
}
plotErrors(errFs$D)
plotErrors(errRs$D)
plotErrors(errI)
# Look ok except for HL-L, and some noisiness in some labs perhaps from low
# read counts
