---
  title: "DADA2 for NAAMES + EXPORTS"
  # Sasha Kramer
  # based on https://benjjneb.github.io/dada2/tutorial.html and https://github.com/dcat4/amplicon_bioinformatics/
  # for my ref: https://docs.google.com/document/d/1biwZQDe_oPKZukS95RLqvApjsRlfPGvntYuT29kgqr0/edit
---

## Clear your workspace:
rm(list=ls())

## Set working directory
setwd("~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/")

## Load packages and libraries
.cran_packages <- c("gridExtra", "knitr", "data.table", "ggplot2")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst], repos='http://cran.us.r-project.org')
}

sapply(c(.cran_packages), require, character.only = TRUE)
library("dada2")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("DECIPHER", "Biostrings", "dada2"))

library(dada2); packageVersion("dada2")

## Set up new directories to raw data and for storing results at each step of data processing:

# write the name of your library here:
Lib1 <- "July2020"
Lib2 <- "Nov2020"
Lib3 <- "Dec2020"

# establish paths to raw miseq files
path1 <- paste0("~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/",Lib1)
path2 <- paste0("~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/",Lib2)
path3 <- paste0("~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/",Lib3)

# designate a path to store all your plots of processing results.
plots_path1 <- paste0("~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/",Lib1,"/Processing/")
plots_path2 <- paste0("~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/",Lib2,"/Processing/")
plots_path3 <- paste0("~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/",Lib3,"/Processing/")

# create directories for storing plots of each library:
dir.create(plots_path1, recursive = TRUE)
dir.create(plots_path2, recursive = TRUE)
dir.create(plots_path3, recursive = TRUE)

#### Workflow for Miseq data
## Step 1:
# Filter and Trim sequences with dada2_ML_filtering_Nov2018_dc.R

## Paths for filtering script (where your raw reads are):
pathF1 <- path1
pathR1 <- path1
pathF2 <- path2
pathR2 <- path2
pathF3 <- path3
pathR3 <- path3

# Check quality scores first from: https://github.com/dcat4/amplicon_bioinformatics/blob/master/dada2_pipe/dada2_ML_inspect_Qscores_Mar20.R

## Library 1: July 2020 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## grab fastq's
fastqFs <- sort(list.files(pathF1, pattern="_R1_001.fastq"))
fastqRs <- sort(list.files(pathR1, pattern="_R2_001.fastq"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

# plot Quality scores; may need to subset if you're library is heavily multiplexed
# saves to a pdf file for review
plot.Rquals <- plotQualityProfile(paste0(pathR1,"/",fastqRs[1:6]))
ggsave(paste0(plots_path1,"Rqualplot1.pdf"), plot.Rquals, device="pdf")
plot.Rquals <- plotQualityProfile(paste0(pathR1,"/",fastqRs[7:12]))
ggsave(paste0(plots_path1,"Rqualplot2.pdf"), plot.Rquals, device="pdf")
plot.Rquals <- plotQualityProfile(paste0(pathR1,"/",fastqRs[13:18]))
ggsave(paste0(plots_path1,"Rqualplot3.pdf"), plot.Rquals, device="pdf")
plot.Rquals <- plotQualityProfile(paste0(pathR1,"/",fastqRs[19:23]))
ggsave(paste0(plots_path1,"Rqualplot4.pdf"), plot.Rquals, device="pdf")

plot.Fquals <- plotQualityProfile(paste0(pathF1,"/",fastqFs[1:6]))
ggsave(paste0(plots_path1,"Fqualplot1.pdf"), plot.Fquals, device="pdf")
plot.Fquals <- plotQualityProfile(paste0(pathF1,"/",fastqFs[7:12]))
ggsave(paste0(plots_path1,"Fqualplot2.pdf"), plot.Fquals, device="pdf")
plot.Fquals <- plotQualityProfile(paste0(pathF1,"/",fastqFs[13:18]))
ggsave(paste0(plots_path1,"Fqualplot3.pdf"), plot.Fquals, device="pdf")
plot.Fquals <- plotQualityProfile(paste0(pathF1,"/",fastqFs[19:23]))
ggsave(paste0(plots_path1,"Fqualplot4.pdf"), plot.Fquals, device="pdf")

## Library 2: Nov 2020 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## grab fastq's
fastqFs <- sort(list.files(pathF2, pattern="_R1_001.fastq"))
fastqRs <- sort(list.files(pathR2, pattern="_R2_001.fastq"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

# plot Quality scores; may need to subset if your library is heavily multiplexed
# saves to a pdf file for review
plot.Rquals <- plotQualityProfile(paste0(pathR2,"/",fastqRs))
ggsave(paste0(plots_path2,"Rqualplot.pdf"), plot.Rquals, device="pdf")

plot.Fquals <- plotQualityProfile(paste0(pathF2,"/",fastqFs))
ggsave(paste0(plots_path2,"Fqualplot.pdf"), plot.Fquals, device="pdf")

## Library 3: Dec 2020 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## grab fastq's
fastqFs <- sort(list.files(pathF3, pattern="_R1_001.fastq"))
fastqRs <- sort(list.files(pathR3, pattern="_R2_001.fastq"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

# plot Quality scores; may need to subset if you're library is heavily multiplexed
# saves to a pdf file for review
plot.Rquals <- plotQualityProfile(paste0(pathR3,"/",fastqRs[1:9]))
ggsave(paste0(plots_path3,"Rqualplot1.pdf"), plot.Rquals, device="pdf")
plot.Rquals <- plotQualityProfile(paste0(pathR3,"/",fastqRs[10:18]))
ggsave(paste0(plots_path3,"Rqualplot2.pdf"), plot.Rquals, device="pdf")
plot.Rquals <- plotQualityProfile(paste0(pathR3,"/",fastqRs[19:27]))
ggsave(paste0(plots_path3,"Rqualplot3.pdf"), plot.Rquals, device="pdf")

plot.Fquals <- plotQualityProfile(paste0(pathF3,"/",fastqFs[1:9]))
ggsave(paste0(plots_path3,"Fqualplot1.pdf"), plot.Fquals, device="pdf")
plot.Fquals <- plotQualityProfile(paste0(pathF3,"/",fastqFs[10:18]))
ggsave(paste0(plots_path3,"Fqualplot2.pdf"), plot.Fquals, device="pdf")
plot.Fquals <- plotQualityProfile(paste0(pathF3,"/",fastqFs[19:27]))
ggsave(paste0(plots_path3,"Fqualplot3.pdf"), plot.Fquals, device="pdf")

# Now filter: https://github.com/dcat4/amplicon_bioinformatics/blob/master/dada2_pipe/dada2_ML_filtering_Mar20_dc.R
## Where you want the filtered reads to go:
filtMethod <- "" # this can be empty if you want, but if not make it an indication of filter parameters you used.
filtpathF1 <- file.path("~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/",Lib1,"filtered_reads",filtMethod) # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR1 <- file.path("~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/",Lib1,"filtered_reads",filtMethod) # ...
filtpathF2 <- file.path("~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/",Lib2,"filtered_reads",filtMethod) # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR2 <- file.path("~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/",Lib2,"filtered_reads",filtMethod) # ...
filtpathF3 <- file.path("~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/",Lib3,"filtered_reads",filtMethod) # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR3 <- file.path("~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/",Lib3,"filtered_reads",filtMethod) # ...

## Where you want the readsin/out .csv file to go:
pathout1 <- paste0("~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/",Lib1,"/Processing/filtering/",filtMethod)
pathout2 <- paste0("~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/",Lib2,"/Processing/filtering/",filtMethod)
pathout3 <- paste0("~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/",Lib3,"/Processing/filtering/",filtMethod)

## Library 1: July 2020 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# File parsing
fastqFs <- sort(list.files(pathF1, pattern="_R1_001.fastq"))
fastqRs <- sort(list.files(pathR1, pattern="_R2_001.fastq"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

#
out <- filterAndTrim(fwd=file.path(pathF1, fastqFs), filt=file.path(filtpathF1, fastqFs),
                     rev=file.path(pathR1, fastqRs), filt.rev=file.path(filtpathR1, fastqRs),
                     truncLen=c(140,120), maxEE=2, truncQ=2, maxN=0, rm.phix=TRUE,
                     compress=TRUE, verbose=TRUE, multithread=TRUE)

dir.create(pathout1)
write.csv(out, paste0(pathout1,"/readsin_readsout.csv"))

## Library 2: Nov 2020 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# File parsing
fastqFs <- sort(list.files(pathF2, pattern="_R1_001.fastq"))
fastqRs <- sort(list.files(pathR2, pattern="_R2_001.fastq"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

#
out <- filterAndTrim(fwd=file.path(pathF2, fastqFs), filt=file.path(filtpathF2, fastqFs),
                     rev=file.path(pathR2, fastqRs), filt.rev=file.path(filtpathR2, fastqRs),
                     truncLen=c(140,120), maxEE=2, truncQ=2, maxN=0, rm.phix=TRUE,
                     compress=TRUE, verbose=TRUE, multithread=TRUE)

dir.create(pathout2)
write.csv(out, paste0(pathout2,"/readsin_readsout.csv"))

## Library 3: Dec 2020 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# File parsing
fastqFs <- sort(list.files(pathF3, pattern="_R1_001.fastq"))
fastqRs <- sort(list.files(pathR3, pattern="_R2_001.fastq"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

#
out <- filterAndTrim(fwd=file.path(pathF3, fastqFs), filt=file.path(filtpathF3, fastqFs),
                     rev=file.path(pathR3, fastqRs), filt.rev=file.path(filtpathR3, fastqRs),
                     truncLen=c(140,120), maxEE=2, truncQ=2, maxN=0, rm.phix=TRUE,
                     compress=TRUE, verbose=TRUE, multithread=TRUE)

dir.create(pathout3)
write.csv(out, paste0(pathout3,"/readsin_readsout.csv"))

## Step 2:
## Infer Sequence Variants with: https://github.com/dcat4/amplicon_bioinformatics/blob/master/dada2_pipe/dada2_ML_infer_seq_variants_Mar20_dc.R

### Paths for this script are the same - using the dataset I filtered above
# same thing as filt method, just keeping track of parameter changes. This gets added onto the filtMethod
# in the directory name:
ISVmethod <- ""

pathout_ISV1 <- paste0("~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/",Lib1,"/sample_inference/",filtMethod,ISVmethod)
pathout_ISV2 <- paste0("~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/",Lib2,"/sample_inference/",filtMethod,ISVmethod)
pathout_ISV3 <- paste0("~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/",Lib3,"/sample_inference/",filtMethod,ISVmethod)

dir.create(pathout_ISV1)
dir.create(pathout_ISV2)
dir.create(pathout_ISV3)

### Library 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# File parsing
filtFs <- list.files(filtpathF1, pattern="_R1_001.fastq", full.names = TRUE)
filtRs <- list.files(filtpathR1, pattern="_R2_001.fastq", full.names = TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)

## learn error rates:
errF <- learnErrors(filtFs, nbases=1e8, multithread=FALSE,randomize=TRUE)

# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=FALSE,randomize=TRUE)

# # visualize and save error rates
# plotErrors(errF, nominalQ=TRUE)
plot.FerrMdl <- plotErrors(errF, nominalQ=TRUE)
ggsave(paste0(plots_path1,"F_errMdl.pdf"), plot.FerrMdl, device="pdf")

# plotErrors(errR, nominalQ=TRUE)
plot.RerrMdl <- plotErrors(errR, nominalQ=TRUE)
ggsave(paste0(plots_path1,"R_errMdl.pdf"), plot.RerrMdl, device="pdf")

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR, trimOverhang = TRUE)
  mergers[[sam]] <- merger
}

rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)

# keep track of reads merged and tabled:
Nreads.mergedAndTabled <- rowSums(seqtab)
write.csv(Nreads.mergedAndTabled, paste0(pathout_ISV1,"readsMergedAndTabled.csv"))

# Inspect distribution of sequence lengths
ampSizes <- table(nchar(getSequences(seqtab)))

saveRDS(seqtab, paste0(pathout_ISV1,"seqtab.rds")) # CHANGE ME to where you want sequence table saved

### Library 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# File parsing
filtFs <- list.files(filtpathF2, pattern="_R1_001.fastq", full.names = TRUE)
filtRs <- list.files(filtpathR2, pattern="_R2_001.fastq", full.names = TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)

## learn error rates:
errF <- learnErrors(filtFs, nbases=1e8, multithread=FALSE,randomize=TRUE)

# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=FALSE,randomize=TRUE)

# # visualize and save error rates
plot.FerrMdl <- plotErrors(errF, nominalQ=TRUE)
ggsave(paste0(plots_path2,"F_errMdl.pdf"), plot.FerrMdl, device="pdf")

plot.RerrMdl <- plotErrors(errR, nominalQ=TRUE)
ggsave(paste0(plots_path2,"R_errMdl.pdf"), plot.RerrMdl, device="pdf")

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR, trimOverhang = TRUE)
  mergers[[sam]] <- merger
}

rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)

# keep track of reads merged and tabled:
Nreads.mergedAndTabled <- rowSums(seqtab)
write.csv(Nreads.mergedAndTabled, paste0(pathout_ISV2,"readsMergedAndTabled.csv"))

# Inspect distribution of sequence lengths
ampSizes <- table(nchar(getSequences(seqtab)))

saveRDS(seqtab, paste0(pathout_ISV2,"seqtab.rds")) # CHANGE ME to where you want sequence table saved


### Library 3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# File parsing
filtFs <- list.files(filtpathF3, pattern="_R1_001.fastq", full.names = TRUE)
filtRs <- list.files(filtpathR3, pattern="_R2_001.fastq", full.names = TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)

## learn error rates:
errF <- learnErrors(filtFs, nbases=1e8, multithread=FALSE,randomize=TRUE)

# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=FALSE,randomize=TRUE)

# # visualize and save error rates
# plotErrors(errF, nominalQ=TRUE)
plot.FerrMdl <- plotErrors(errF, nominalQ=TRUE)
ggsave(paste0(plots_path3,"F_errMdl.pdf"), plot.FerrMdl, device="pdf")

# plotErrors(errR, nominalQ=TRUE)
plot.RerrMdl <- plotErrors(errR, nominalQ=TRUE)
ggsave(paste0(plots_path3,"R_errMdl.pdf"), plot.RerrMdl, device="pdf")

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR, trimOverhang = TRUE)
  mergers[[sam]] <- merger
}

rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)

# keep track of reads merged and tabled:
Nreads.mergedAndTabled <- rowSums(seqtab)
write.csv(Nreads.mergedAndTabled, paste0(pathout_ISV3,"readsMergedAndTabled.csv"))

# Inspect distribution of sequence lengths
ampSizes <- table(nchar(getSequences(seqtab)))

saveRDS(seqtab, paste0(pathout_ISV3,"seqtab.rds")) # CHANGE ME to where you want sequence table saved

## Step 3: Combine libraries and remove chimeras with: https://github.com/dcat4/amplicon_bioinformatics/blob/master/dada2_pipe/dada2_ML_remove_chimera_Mar20.R

projectName <- "all18s_naames_exp"

# Where you want the output sequence table saved:
RCmethod <- "" # same as above, there's only 2 options here though really.
pathout_RC <- paste0("~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/",projectName,"/chimera_removal/",filtMethod,ISVmethod)

# Merge multiple runs (if necessary)
st1 <- readRDS(paste0(pathout_ISV1,"seqtab.rds"))
st2 <- readRDS(paste0(pathout_ISV2,"seqtab.rds"))
st3 <- readRDS(paste0(pathout_ISV3,"seqtab.rds"))

st.all <- mergeSequenceTables(st1, st2, st3) # merge all runs

seqtab.nochime <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)

Nreads.nochime <- rowSums(seqtab.nochime)

dir.create(pathout_RC,recursive = TRUE)
saveRDS(seqtab.nochime, paste0(pathout_RC,"seqtab.rds")) # CHANGE ME to where you want sequence table saved
write.csv(Nreads.nochime,paste0(pathout_RC,"readsChimeraFiltered.csv"))
