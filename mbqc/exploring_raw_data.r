library(phyloseq)

ps.all = readRDS(file.path('~/active_research/metagenomics_calibration/mbqc/data/',
                           'mbqc_phyloseq.rds'))
meta = sample_data(ps.all)
meta = subset_samples(meta, (MBQC.ID == 'DZ15296') & (sequencing_wetlab != 'HL-A'))
# saveRDS(meta,
#         file.path('~/active_research/metagenomics_calibration/mbqc/data/',
#                   "mbqc_sample_data.rds"))
df = unique(meta[, c('Bioinformatics.ID', 
                     'extraction_wetlab', 'sequencing_wetlab',
                     "X16S_primer",
                     "seq_machine", 'paired_end_reads', "read_length",
                     "frac_quality_bases", "phiX_frac", 'read_count'
                     # "seq_chem_version",
                     # "extraction_kit_maker", "extraction_kit_model",
                     )])
df$mbqc_name = sample_names(df)
target.length = 10
df$Bioinformatics.ID = stringr::str_pad(df$Bioinformatics.ID,
                                        target.length, side='left',
                                        pad='0')
sample_names(df) = df$Bioinformatics.ID
#saveRDS(df, file.path(path, "blinded_sequence_table.rds"))

# Init
library(phyloseq)
library(dada2);packageVersion("dada2")
library(stringr)
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
names(fnIs) <- sapply(strsplit(fnIs, "[_/]"), `[`, 3)
names(fnFs) <- sapply(strsplit(fnFs, "[_/]"), `[`, 3)
names(fnRs) <- sapply(strsplit(fnRs, "[_/]"), `[`, 3)
identical(names(fnRs), names(fnFs))

# Get data in right order (Not sure if this is needed)
all(rownames(df) %in% c(names(fnFs), names(fnIs))) # TRUE
all(rownames(df)[df$sequencing_wetlab!='HL-I'] %in% names(fnRs)) # TRUE
fnFs <- fnFs[df$Bioinformatics.ID]
fnRs <- fnRs[df$Bioinformatics.ID]

# Inspect HL-C and remove primer sequences
# test.samples = sample(rownames(df)[df$sequencing_wetlab=='HL-C'], 3)
test.samples = c('3762855612', '4963000700', '6905611088')

## Look for primer sequences
# Work with dereplicated sequences
drpF <- derepFastq(fnFs[test.samples], verbose=TRUE)
drpR <- derepFastq(fnRs[test.samples], verbose=TRUE)
# Check lengths
seqsF = names(drpF[[1]]$uniques)
seqsR = names(drpR[[1]]$uniques)
table(sapply(seqsF, nchar)) # 243
table(sapply(seqsR, nchar)) # 251
# Look for primer sequences, starting with the forward reads
str_sub(seqsF[1:10], 1, 22)
# GTGTGCCAGCAGCCGCGGTAAT (observed)
# GTG--CCAGCMGCCGCGGTAA  (F515 primer with 2bp insertion)
# So we'll need to trim the first 21bp
# How many sequences match?
matches = str_sub(seqsF, 1, 21) %in% c("GTGTGCCAGCAGCCGCGGTAA", "GTGTGCCAGCCGCCGCGGTAA")
sum(matches); length(seqsF) # 10856 out of 12819 match
# What's going on with the ones that don't match?
str_sub(seqsF[!matches][1:10], 1, 21)
# Some don't match because of indels or SNPS. Others seem to be totally
# different sequences

seqsF = DNAStringSet(seqsF)
seqsR = DNAStringSet(seqsR)
# mean(vcountPattern(DNAString(F515_HLC), seqsF, max.mismatch=4, min.mismatch=0,
#                    with.indels=TRUE, fixed=FALSE, algorithm="auto"))
qplot(c(neditStartingAt(DNAString(F515), seqsF, with.indels=TRUE, fixed=FALSE)))
qplot(c(neditStartingAt(DNAString(F515_HLC), seqsF, with.indels=TRUE, fixed=FALSE)))
qplot(c(neditStartingAt(DNAString(R806), seqsR, with.indels=TRUE, fixed=FALSE)))

# Remove primers with cutadapt
# cutadapt -g ^GTGCCAGCMGCCGCGGTAA -o output.fastq input.fastq
# Might want to use the sequence with the two inserted bases that actually
# appears most often in the reads

for (fn in fnFs[test.samples]) {
    outfile = file.path(path, "primers_removed", basename(fn))
    command = paste('cutadapt -g', paste0('^', F515_HLC), '-e 0.2', 
                    '-o', outfile, fn)
    system(command)
}

# Reverse reads
str_sub(seqsR[1:10], 1, 20) # Matches R806

# Cranking the error rate up to 60% lets us match 99%+ of sequences, with still
# a very low expected number of matches by chance

for (fn in fnRs[test.samples]) {
    outfile = file.path(path, "primers_trimmed", basename(fn))
    command = paste('cutadapt -g', paste0('^', R806), '-e 0.6', 
                    '-o', outfile, fn)
    system(command)
}

# Trim paired reads together
for (sn in test.samples) {
    fnF = fnFs[sn]
    fnR = fnRs[sn]
    outfileF = file.path("primers_trimmed", basename(fnF))
    outfileR = file.path("primers_trimmed", basename(fnR))
    untrimmedF = file.path("untrimmed", basename(fnF))
    untrimmedR = file.path("untrimmed", basename(fnR))
    command = paste('cutadapt', 
                    '-g', paste0('^', F515_HLC), '-e 0.2', 
                    '-G', paste0('^', R806),
                    '-o', outfileF, 
                    '-p', outfileR, 
                    paste0('--untrimmed-output=', untrimmedF),
                    paste0('--untrimmed-paired-output=', untrimmedR),
                    fnF, fnR)
    system(command)
}

# primer-removed sequence fastq files
prFs <- file.path(path, "primers_trimmed", basename(fnFs[test.samples]))
prRs <- file.path(path, "primers_trimmed", basename(fnRs[test.samples]))
names(prFs) = test.samples
names(prRs) = test.samples

# Check lengths
drpF <- derepFastq(prFs[test.samples], verbose=TRUE)
drpR <- derepFastq(prRs[test.samples], verbose=TRUE)
# Check lengths
seqsF = names(drpF[[1]]$uniques)
seqsR = names(drpR[[1]]$uniques)
table(sapply(seqsF, nchar)) # 219-226
table(sapply(seqsR, nchar)) # 227-235

# Check quality profiles
# Quality stays high, so I'll choose truncation length to remove 5bp from the ends of
# the shortest trimmed reads.
qpf = plotQualityProfile(prFs)
qpr = plotQualityProfile(prRs)
qpf + geom_vline(xintercept=214)
qpr + geom_vline(xintercept=222)

# Filter and trim
filtFs <- file.path(path, "filtered", basename(prFs))
names(filtFs) = names(prFs)
filtRs <- file.path(path, "filtered", basename(prRs))
names(filtRs) = names(prRs)

outC <- filterAndTrim(prFs[test.samples], filtFs[test.samples],
                      prRs[test.samples], filtRs[test.samples],
                      truncLen=c(214, 222), maxEE=1, rm.phix=TRUE, multithread=TRUE)
outC; median(outC[,2]/outC[,1]) # 83% of reads kept

### Run dada pipeline
# Dereplicate
drpF <- derepFastq(filtFs[test.samples], verbose=TRUE)
drpR <- derepFastq(filtRs[test.samples], verbose=TRUE)
# Learn Errors
errFC <- learnErrors(drpF, multithread=TRUE)
errRC <- learnErrors(drpR, multithread=TRUE)

# Run dada
ddFC <- dada(drpF, err=errFC, pool=TRUE, multithread=TRUE)
ddRC <- dada(drpR, err=errRC, pool=TRUE, multithread=TRUE)

# Merge
mmC <- mergePairs(ddFC, drpF, ddRC, drpR, verbose=TRUE)
# >98% successfully merged for the 3 test samples

# Make table
stC <- makeSequenceTable(mmC)
table(nchar(getSequences(stC)))
# Almost all sequences are 253 +-1, with 4 sequences that are ~350
# I wonder if being longer than 252 suggests an extra base from the primers
# made it through?

# Remove chimeras
st <- removeBimeraDenovo(stC, multithread=TRUE, verbose=TRUE)
# Identified 51 bimeras out of 440 input sequences.
sum(st)/sum(stC) # 0.9727925
table(nchar(getSequences(st)))
# Still three unsually long sequences (~350bp). These should be removed.

# Remove usually long sequences; only a tiny fraction of reads
st1 = st[, nchar(getSequences(st)) < 260]
sum(st1)/sum(st) # 0.9999068
table(nchar(getSequences(st1)))

# Check that the differences in length are not due to errors in primer trimming
library(DECIPHER)
seqs <- getSequences(st1)
names(seqs) <- seqs
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

##############################################################################
path = "~/data/mbqc/blinded_sequence_data"
setwd(path)

## Check if the two sequencing runs in lab E have similar QPs

## Check if the G and H samples sequenced in lab F have similar QPs
g.samples = rownames(subset(df, extraction_wetlab=='HL-G'))
f.samples = rownames(subset(df, extraction_wetlab=='HL-F'))
# g.qps = plotQualityProfile(fnFs[g.samples])
# f.qps = plotQualityProfile(fnFs[f.samples])
qps = plotQualityProfile(fnFs[c(g.samples, f.samples)])
# Look very similar
