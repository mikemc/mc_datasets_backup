# Read in metadata and get the sequencing data
df <- read.table(file.path(path, "male_stool_36y.tsv"), header=TRUE, colClasses="character")
path <- "~/MBQC"
setwd(path)
urlFs <- paste0("http://downloads.ihmpdcc.org/data/MBQCBD2/", df$sequencing_wetlab, "/R1/", df$Bioinformatics.ID, "_R1.fastq.gz")
urlRs <- paste0("http://downloads.ihmpdcc.org/data/MBQCBD2/", df$sequencing_wetlab, "/R2/", df$Bioinformatics.ID, "_R2.fastq.gz")
###for(url in c(urlFs, urlRs)) {
###  system(paste("wget", url))
###}

# Init
library(dada2);packageVersion("dada2")
F515 <- "GTGCCAGCMGCCGCGGTAA"
R806 <- "GGACTACHVGGGTWTCTAAT"
fnFs <- list.files(path, pattern="_R1.fastq.gz$")
fnRs <- list.files(path, pattern="_R2.fastq.gz$")
names(fnFs) <- sapply(strsplit(fnFs, "_"), `[`, 1)
names(fnRs) <- sapply(strsplit(fnRs, "_"), `[`, 1)
identical(names(fnRs), names(fnFs))

# Get data in right order
rownames(df) <- df$Bioinformatics.ID
all(rownames(df) %in% names(fnFs)) # TRUE
fnFs <- fnFs[df$Bioinformatics.ID]
fnRs <- fnRs[df$Bioinformatics.ID]
isB <- df$sequencing_wetlab == "HL-B"
isE <- df$sequencing_wetlab == "HL-E"
#saveRDS(df, file.path(path, "RDS", "df.rds"))

### Inspect and Filter
plotQualityProfile(sample(fnFs[isB], 3)) # 220 truncation length
plotQualityProfile(sample(fnRs[isB], 3)) # 160
# Miseq 2x250 looks like

plotQualityProfile(sample(fnFs[isE], 3)) # 145
plotQualityProfile(sample(fnRs[isE], 3)) # 140
# Hiseq(?) 2x150 looks like

# Filter and Trim as appropriate for each sequencing site
filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))
outB <- filterAndTrim(fnFs[isB], filtFs[isB], fnRs[isB], filtRs[isB],
                      truncLen=c(220, 160), maxEE=1, rm.phix=TRUE, multithread=TRUE)
outB; median(outB[,2]/outB[,1]) # 84% of reads kept

outE <- filterAndTrim(fnFs[isE], filtFs[isE], fnRs[isE], filtRs[isE],
                      truncLen=c(145, 140), maxEE=1, rm.phix=TRUE, multithread=TRUE)
outE; median(outE[,2]/outE[,1]) # 92% of reads kept

### Run dada pipeline
# Dereplicate
drpF <- derepFastq(filtFs, verbose=TRUE)
drpR <- derepFastq(filtRs, verbose=TRUE)
# Learn Errors
errFB <- learnErrors(drpF[isB], multithread=TRUE) # 750 elapsed, in seconds, by system.time()
errRB <- learnErrors(drpR[isB], multithread=TRUE) # 611 elapsed
errFE <- learnErrors(drpF[isE], multithread=TRUE) # 494 elapsed
errRE <- learnErrors(drpR[isE], multithread=TRUE) # 1181 elapsed
#saveRDS(errFB, file.path(path, "RDS", "errFB.rds"))
#saveRDS(errRB, file.path(path, "RDS", "errRB.rds"))
#saveRDS(errFE, file.path(path, "RDS", "errFE.rds"))
#saveRDS(errRE, file.path(path, "RDS", "errRE.rds"))

# Run dada
ddFB <- dada(drpF[isB], err=errFB, pool=TRUE, multithread=TRUE) # 450 elapsed
ddRB <- dada(drpR[isB], err=errRB, pool=TRUE, multithread=TRUE) # 335 elapsed
ddFE <- dada(drpF[isE], err=errFE, pool=TRUE, multithread=TRUE) # 337 elapsed
ddRE <- dada(drpR[isE], err=errRE, pool=TRUE, multithread=TRUE) # 800 elapsed
# pool=TRUE here for maximum comparability across samples
#saveRDS(ddFB, file.path(path, "RDS", "ddFB.rds"))
#saveRDS(ddRB, file.path(path, "RDS", "ddRB.rds"))
#saveRDS(ddFE, file.path(path, "RDS", "ddFE.rds"))
#saveRDS(ddRE, file.path(path, "RDS", "ddRE.rds"))

# Merge
mmB <- mergePairs(ddFB, drpF[isB], ddRB, drpR[isB], verbose=TRUE)
# >95% of reads successfully merged
mmE <- mergePairs(ddFE, drpF[isE], ddRE, drpR[isE], verbose=TRUE)
# >90% of reads successfully merged

# Make table
stB <- makeSequenceTable(mmB)
stE <- makeSequenceTable(mmE)
table(nchar(getSequences(stB)))
table(nchar(getSequences(stE)))
# Both almost all in 251-254 expected range
# stB a few longer, stE a few shorter, but very few
sta <- mergeSequenceTables(stB, stE)
rownames(sta) <- sapply(strsplit(rownames(sta), "_"), `[`, 1)
identical(rownames(sta), df$Bioinformatics.ID) # TRUE
saveRDS(sta, file.path(path, "RDS", "sta.rds"))

# Remove chimeras
st <- removeBimeraDenovo(sta, multithread=TRUE, verbose=TRUE)
sum(st[isB,])/sum(sta[isB,]) # 93.8%
sum(st[isE,])/sum(sta[isE,]) # 73.5%
saveRDS(st, file.path(path, "RDS", "st.rds"))
# That is a high chimera fraction in lab E!!

### Filter SVs based on prevalence
plot(colSums(st), log="y")
plot(colSums(st>0))
table(colSums(st[isB,]>0), colSums(st[isE,]>0))
# May be a "bad" labE sample, as a big prevalance spiked in 10B/9E
# Also a lot more single-labE-sample SVs than single-labB-sample
rowSums(st) # 4775275923 has 30k reads, well below any others (min 74k)
plot(seq(20), sapply(seq(20), function(min.prev) sum(st[,colSums(st>0)>=min.prev])), 
     xlab="Minimum Prevalence", ylab="Total Reads")
# From this and table above, choosing a minimum prevalence of 7
ss <- st[,colSums(st>0) >= 7]
ncol(ss);ncol(st);sum(ss)/sum(st) # ~99% of the reads kept after prevalence filter
#saveRDS(ss, file.path(path, "RDS", "ss.rds"))

# Assign taxonomy
tax <- assignTaxonomy(ss, "~/tax/silva_nr_v128_train_set.fa.gz", multithread=TRUE)
table(tax[,"Phylum"]) # ASV numbers dominated by Firmicutes (437/525)
rowsum(colSums(ss), tax[,"Phylum"])/sum(ss)
# But reads not nearly as much: 56% Firmicutes, 40% Bacteroidetes
spec <- assignSpecies(ss, "~/tax/silva_species_assignment_v128.fa.gz", allowMultiple=TRUE)
unname(spec[1:20,])
# Quite a lot of the abundant stuff can be assigned to species (or at least multiple)
#saveRDS(tax, file.path(path, "RDS", "tax.rds"))
#saveRDS(spec, file.path(path, "RDS", "spec.rds"))

### Analysis (can use readRDS to load in ss, tax, df)
path <- "~/MBQC"
ss <- readRDS(file.path(path, "RDS", "ss.rds"))
df <- readRDS(file.path(path, "RDS", "df.rds"))
tax <- readRDS(file.path(path, "RDS", "tax.rds"))
ff <- sweep(ss, 1, rowSums(ss), "/") # Convert to frequency
medians <- apply(ff, 2, median) # Median frequency for each SV
bb <- log(sweep(ff, 2, medians, "/")) # Bias = log(freq/median)
plot(factor(df$sequencing_wetlab), bb[,1])
plot(factor(df$extraction_wetlab), bb[,1])
# Wow, huge effect driven by extraction for that ASV

# Make a plot for the top 20 ASVs
library(ggplot2); library(reshape2)
pdf <- cbind(df, bb[,1:20])
pdf <- melt(pdf, id.vars=colnames(df), variable.name="sequence", value.name="Bias")
pdf <- cbind(pdf, tax[pdf$sequence,])
pdf$ASV <- paste0("ASV", match(pdf$sequence, colnames(ss)))
ggplot(data=pdf, aes(x=extraction_wetlab, y=Bias, color=Genus)) + 
  geom_boxplot() + geom_jitter() + facet_wrap(~ASV) + theme_bw()
# Looks quite promising!

# Make a phylogenetic tree (from F1000 workflow)
seqs <- colnames(ss)
names(seqs) <- seqs # names propagate to the tip labels of the tree
library(DECIPHER); library(phangorn)
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
system.time(fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, 
            rearrangement="stochastic", control=pml.control(trace=0))) #874 elapsed
saveRDS(fitGTR, file.path(path, "RDS", "fitGTR.rds"))

library(phyloseq)
psa <- phyloseq(tax_table(tax), sample_data(df), otu_table(ss, taxa_are_rows = FALSE),
               phy_tree(fitGTR$tree))
psa # Got it

plot_ordination(psa, ordinate(psa, method="MDS"), color="extraction_wetlab") + geom_text(aes(label=Bioinformatics.ID))
# One outlier HL-E sample, perhaps should be removed, 4775275923
sample_sums(psa) # Yeah that's the small one, down to 13k from 30k before removing <7 prevalance
ps <- subset_samples(psa, !Bioinformatics.ID %in% "4775275923"); ps
plot_ordination(ps, ordinate(ps, method="MDS", distance="bray"), color="extraction_wetlab", shape="sequencing_wetlab")
# Clear separation based on central or other extraction
plot_ordination(ps, ordinate(ps, method="MDS", distance="w-unifrac"), color="extraction_wetlab", shape="sequencing_wetlab")
# Same story, except central extractions even more tightly clustered now
plot_ordination(ps, ordinate(ps, method="MDS", distance="unifrac"), color="extraction_wetlab", shape="sequencing_wetlab")
# Can see clustering of sequencing lab w/in NA extraction most clearly here - contaminants?



