# Init
library(stringr)
library(phyloseq); library(dada2);packageVersion("dada2")
library(Biostrings); 
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

# # Get data in right order (Not sure if this is needed)
# all(rownames(df) %in% c(names(fnFs), names(fnIs))) # TRUE
# all(rownames(df)[df$sequencing_wetlab!='HL-I'] %in% names(fnRs)) # TRUE
# fnFs <- fnFs[df$Bioinformatics.ID]
# fnRs <- fnRs[df$Bioinformatics.ID]

#### Look for primer seqs

# Load an example set of F and R sequences for each lab
seqlabs = unique(df$sequencing_wetlab)
set.seed(1)
split.samples = split(df$Bioinformatics.ID, df$sequencing_wetlab)
test.samples = sapply(split.samples, sample, 1)

# Dereplicate
drpF <- derepFastq(fnFs[test.samples[names(test.samples)!='HL-I']],
                   verbose=TRUE)
drpR <- derepFastq(fnRs[test.samples[names(test.samples)!='HL-I']],
                   verbose=TRUE)
drpI <- derepFastq(fnIs[test.samples[names(test.samples)=='HL-I']],
                   verbose=TRUE)

## Check edit distances for primers
# Only HL-C has primer sequences, which makes sense since I think it used an
# in-house primer set while the others used EMP or Schloess 2013
table(df[,c('sequencing_wetlab', 'X16S_primer')])

# Check all labs with paired reads
for (lab in seqlabs[seqlabs!='HL-I']) {
    seqsF = DNAStringSet(names(drpF[[test.samples[lab]]]$uniques))
    seqsR = DNAStringSet(names(drpR[[test.samples[lab]]]$uniques))
    print(lab)
    print(summary(c(neditStartingAt(DNAString(F515), seqsF, with.indels=TRUE,
                                    fixed=FALSE))))
    print(summary(c(neditStartingAt(DNAString(R806), seqsR, with.indels=TRUE,
                                    fixed=FALSE))))
}

# HL-C: primers present; HL-C specific forward primer
seqsF = DNAStringSet(names(drpF[[test.samples['HL-C']]]$uniques))
seqsR = DNAStringSet(names(drpR[[test.samples['HL-C']]]$uniques))
summary(c(neditStartingAt(DNAString(F515), seqsF, with.indels=TRUE, fixed=FALSE)))
summary(c(neditStartingAt(DNAString(F515_HLC), seqsF, with.indels=TRUE, fixed=FALSE)))
summary(c(neditStartingAt(DNAString(R806), seqsR, with.indels=TRUE, fixed=FALSE)))
# qplot(c(neditStartingAt(DNAString(F515), seqsF, with.indels=TRUE, fixed=FALSE)))
# qplot(c(neditStartingAt(DNAString(F515_HLC), seqsF, with.indels=TRUE, fixed=FALSE)))
# qplot(c(neditStartingAt(DNAString(R806), seqsR, with.indels=TRUE, fixed=FALSE)))

# HL-I: single-end reads, no primers
seqsI = DNAStringSet(names(drpI$uniques))
summary(c(neditStartingAt(DNAString(F515), seqsI, with.indels=TRUE, fixed=FALSE)))
# qplot(c(neditStartingAt(DNAString(F515), seqsI, with.indels=TRUE, fixed=FALSE)))
# Min distance is 3, but almost all are >=8

#### Trim primers from HL-C samples

# Trim paired reads together. Only reads where both ends match the primer
# sequence are written into the primers_trimmed folder
for (sn in split.samples[['HL-C']]) {
    fnF = fnFs[sn]
    fnR = fnRs[sn]
    outfileF = file.path("primers_trimmed", basename(fnF))
    outfileR = file.path("primers_trimmed", basename(fnR))
    untrimmedF = file.path("untrimmed", basename(fnF))
    untrimmedR = file.path("untrimmed", basename(fnR))
    command = paste('cutadapt', 
                    '-g', paste0('^', F515_HLC), '-e 0.1',
                    '-G', paste0('^', R806),
                    '-o', outfileF, 
                    '-p', outfileR, 
                    paste0('--untrimmed-output=', untrimmedF),
                    paste0('--untrimmed-paired-output=', untrimmedR),
                    fnF, fnR)
    system(command)
}

# ~99% of R1 reads match the primer sequence, but only ~90% of R2 reads match
# the primer sequence at this low error threshold.  The lowest percentage was
# was 87.6% for sample 7826871407
