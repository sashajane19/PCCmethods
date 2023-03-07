---
  title: "Run Catlett ensembleTax on NAAMES and EXPORTS"
  # Sasha Kramer, based on code by Dylan Catlett et al.
  # https://github.com/dcat4/ensembleTax
  # results of DADA2 filtering/trimming are saved: seqtab.rds in ~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/all18S_naames_exp/chimera_removal
---

## Picking up at end of DADA2_NAAMES_EXPORTS.R so wd and directories should be the same
setwd("~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/")

## Load libraries and data
devtools::install_github("pr2database/pr2database")
install.packages(c("dplyr"))
install.packages(c("stringr"))
install.packages(c("usethis"))
install.packages(c("devtools"))
install.packages(c("knitr"))
install.packages(c("rmarkdown"))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(devtools)
devtools::install_github("dcat4/ensembleTax", build_manual = TRUE, build_vignettes = TRUE)
packageVersion("ensembleTax")

library("ensembleTax")
library("Biostrings")
library("stringr")

# DADA2 assignTaxonomy with silva:
bayes.silva <- assignTaxonomy(seqtab.nochime, "~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/tax/silva_nr_v138_train_set.fa", multithread=TRUE, minBoot = 0, outputBootstraps=TRUE)

# DECIPHER idtaxa with silva:
library(DECIPHER); packageVersion("DECIPHER")
dna <- DNAStringSet(getSequences(seqtab.nochime)) # Create a DNAStringSet from the ASVs
load("~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/tax/SILVA_SSU_r138_2019.RData")
idtax.silva <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors

# DADA2 assignTaxonomy with PR2:
bayes.pr2 <- assignTaxonomy(seqtab.nochime, "~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/tax/pr2_version_4.14.0_SSU_dada2.fasta", multithread = TRUE, taxLevels
                            = c("Kingdom","Supergroup","Division","Class","Order","Family","Genus","Species"), minBoot = 0, outputBootstraps=TRUE)

# DECIPHER idtaxa with PR2:
trainset <- readRDS("~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/tax/pr2_version_4.14.0_SSU.decipher.trained.rds")
idtax.pr2 <- IdTaxa(dna, trainset, strand="top", processors=NULL, verbose=FALSE) # use all processors

# Save files in case future work is lost...
saveRDS(bayes.silva,"~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/all18s_naames_exp/initial_tax_tabs/bayes_silva_Sept21.rds")
saveRDS(bayes.pr2,"~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/all18s_naames_exp/initial_tax_tabs/bayes_pr2_Sept21.rds")
saveRDS(idtax.silva,"~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/all18s_naames_exp/initial_tax_tabs/idtax_silva_Sept21.rds")
saveRDS(idtax.pr2,"~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/all18s_naames_exp/initial_tax_tabs/idtax_pr2_Sept21.rds")

# Create an ASV "rubric" to track ASVs through the ensembleTax pipeline
asv.rubric <- DNAStringSet(getSequences(seqtab.nochime))
# this creates names (sv1, sv2, ..., svX) for each ASV
snam <- vector(mode = "character", length = length(asv.rubric))
for (i in 1:length(asv.rubric)) {
  snam[i] <- paste0("sv", as.character(i))
}
names(asv.rubric) <- snam

# Put data into the same formats for next steps
bayes.pr2 <- bayestax2df(bayes.pr2, db = "pr2", boot = 60, rubric = asv.rubric, return.conf = FALSE)
bayes.silva <- bayestax2df(bayes.silva, db = "silva", boot = 60, rubric = asv.rubric, return.conf = FALSE)
# go to Sasha_idtaxpr2_fix.R --> idtax.pr2 <- idtax2df(idtax.pr2, db = "pr2", boot = 50, rubric = asv.rubric, return.conf = FALSE)
idtax.silva <- idtax2df(idtax.silva, db = "silva", boot = 50, rubric = asv.rubric, return.conf = FALSE)

# sanity check
identical(bayes.pr2[, 1:2], bayes.silva[, 1:2])
identical(bayes.pr2[, 1:2], idtax.silva[, 1:2])
identical(bayes.pr2[, 1:2], idtax.pr2[, 1:2])

# Switch here to etax_fullpipe_furreal.R from Dylan

# Now continue with Dylan's pipeline from github at line 172
# Remove proks, metazoa, fungi, streptophyta
rm.i <- unique(c(which(etax.pr2$kingdom %in% c("Bacteria", "Archaea")),
                 which(etax.pr2$kingdom %in% c("Bacteria", "Archaea")),
                 which(etax.pr2$division %in% c("Metazoa", "Fungi", "Streptophyta", "Rhodophyta")),
                 which(etax.pr2$class %in% c("Phaeophyceae" , "Ulvophyceae"))))
bayes.silva.map <- bayes.silva.map[-rm.i , ]
idtax.silva.map <- idtax.silva.map[-rm.i , ]
bayes.pr2 <- bayes.pr2[-rm.i , ]
idtax.pr2 <- idtax.pr2[-rm.i , ]
etax.pr2 <- etax.pr2[-rm.i , ]

bayes.silva <- bayes.silva[-rm.i , ]
idtax.silva <- idtax.silva[-rm.i , ]
bayes.pr2.map <- bayes.pr2.map[-rm.i , ]
idtax.pr2.map <- idtax.pr2.map[-rm.i , ]
etax.silva <- etax.silva[-rm.i , ]

# Compute and plot % ASVs unassigned at each rank:
library("ggplot2")
library("reshape2")
nasum <- function(taxdf){
  notuz <- nrow(taxdf)
  x <- is.na(taxdf)
  ii <- colSums(x) / notuz
  return(ii)
}
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
# with the individuals and 4-way ensembles:
xx.pr2 <- list(idtax.pr2, bayes.pr2,
               idtax.silva.map, bayes.silva.map,
               etax.pr2)
names(xx.pr2) <- c("idtax-pr2", "bayes-pr2",
                   "idtax-silva","bayes-silva",
                   "ensemble-res")
xx.pr2 <- lapply(xx.pr2, function(x) x[, -c(1,2)])
tina <- lapply(xx.pr2, nasum)
yaboi <- matrix(unlist(tina), nrow=length(xx.pr2), byrow=TRUE)
rownames(yaboi) <- names(xx.pr2)
colnames(yaboi) <- colnames(xx.pr2[[1]])
yaboi <- as.data.frame(t(yaboi))
yaboi$rankz <- rownames(yaboi)
yaboi <- melt(yaboi, id.vars = "rankz")
p.pr2.2 <- ggplot(yaboi, aes(fill = variable, x = rankz, y = value)) +
  geom_bar(stat="identity", color = "black", position=position_dodge(width=0.8)) +
  labs(x = "Taxonomic Rank", y = "Proportion of ASVs Unassigned") +
  scale_x_discrete(limits = colnames(xx.pr2[[1]])) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_text(size = 12, face="bold"),
        axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12, face="bold"),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "white"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "white"),
        axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_fill_manual(name = "Taxonomy table",
                    breaks = names(xx.pr2),
                    values = cbPalette)
ggsave(paste0("~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/Figures/","ASVs_unassigned_PR2.pdf"), p.pr2.2, device="pdf",width=11,height=8.5,units="in")

xx.silva <- list(idtax.pr2.map, bayes.pr2.map,
                 idtax.silva, bayes.silva,
                 etax.silva)
names(xx.silva) <- c("idtax-pr2", "bayes-pr2",
                     "idtax-silva","bayes-silva",
                     "ensemble-res")
xx.silva <- lapply(xx.silva, function(x) x[, -c(1,2)])
tina <- lapply(xx.silva, nasum)
yaboi <- matrix(unlist(tina), nrow=length(xx.silva), byrow=TRUE)
rownames(yaboi) <- names(xx.silva)
colnames(yaboi) <- colnames(xx.silva[[1]])
yaboi <- as.data.frame(t(yaboi))
yaboi$rankz <- rownames(yaboi)
yaboi <- melt(yaboi, id.vars = "rankz")
p.silva.2 <- ggplot(yaboi, aes(fill = variable, x = rankz, y = value)) +
  geom_bar(stat="identity", color = "black", position=position_dodge(width=0.8)) +
  labs(x = "Taxonomic Rank", y = "Proportion of ASVs Unassigned") +
  scale_x_discrete(limits = colnames(xx.silva[[1]])) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_text(size = 12, face="bold"),
        axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12, face="bold"),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "white"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "white"),
        axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_fill_manual(name = "Taxonomy table",
                    breaks = names(xx.pr2),
                    values = cbPalette)
ggsave(paste0("~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/Figures/","ASVs_unassigned_silva.pdf"), p.silva.2, device="pdf",width=11,height=8.5,units="in")

# Compare assignments between the ensembles and the individual taxonomy tables:
library("dplyr")
# fcn that compares two taxonomy tables:
tblcomper <- function(y,x) {
  perf <- dplyr::intersect(x, y) # perfectly-matching rows (assignments)
  tmp.x <- dplyr::setdiff(x, perf)
  tmp.y <- dplyr::setdiff(y, perf)
  if (!identical(tmp.x[, 1], tmp.y[, 1])) {
    tmp.x <- sort_my_taxtab(tmp.x, ranknames = colnames(tmp.x)[3:ncol(tmp.x)])
    tmp.y <- sort_my_taxtab(tmp.y, ranknames = colnames(tmp.y)[3:ncol(tmp.x)])
  }
  xna <- is.na(tmp.x)
  yna <- is.na(tmp.y)
  # This warning happens 1 million times when there's no NA's. It doesn't matter b/c inf is always bigger so suppressing.
  ## Warning in min(which(z)): no non-missing arguments to min; returning Inf
  x.minna <- apply(xna, MARGIN = 1, function(z) suppressWarnings(min(which(z))))
  y.minna <- apply(yna, MARGIN = 1, function(z) suppressWarnings(min(which(z))))
  x.mo <- which(x.minna > y.minna)
  y.mo <- which(y.minna > x.minna)

  yunder.i <- c()
  yover.i <- c()
  mis.i <- c()

  # subset where x has more resolved assignments and only to cols where
  # both have assignments, then see if names match
  yunder <- 0
  ymis <- 0
  if (length(x.mo) > 0){
    for (i in 1:length(x.mo)){
      if ((y.minna[x.mo[i]]-1) < 3){ # if kingdom is unassigned just add to under-classification
        yunder <- yunder+1
        # yunder.i <- append(yunder.i, tmp.x$svN)
      } else {
        tmp.tmpx <- tmp.x[x.mo[i] , 3:(y.minna[x.mo[i]]-1)]
        tmp.tmpy <- tmp.y[x.mo[i] , 3:(y.minna[x.mo[i]]-1)]

        if (all(tmp.tmpx == tmp.tmpy)) {
          yunder <- yunder+1
          # yunder.i <- append(yunder.i, tmp.x$svN)
        } else {
          ymis <- ymis+1
          # mis.i <- append(mis.i, tmp.x$svN)
        }
      }
    }
  }

  # repeat above where y is more resolved than x:
  yover <- 0
  xmis <- 0
  if (length(y.mo) > 0){
    for (i in 1:length(y.mo)){
      if ((x.minna[y.mo[i]]-1) < 3){ # if kingdom is unassigned just add to under-classification
        yover <- yover+1
        # yover.i <- append(yover.i, tmp.x$svN)
      } else {
        tmp.tmpx <- tmp.x[y.mo[i] , 3:(x.minna[y.mo[i]]-1)]
        tmp.tmpy <- tmp.y[y.mo[i] , 3:(x.minna[y.mo[i]]-1)]

        if (all(tmp.tmpx == tmp.tmpy)) {
          yover <- yover+1
          # yover.i <- append(yover.i, tmp.x$svN)
        } else {
          xmis <- xmis+1
          # mis.i <- append(mis.i, tmp.x$svN)
        }
      }
    }
  }

  if (yover+xmis != length(y.mo) || yunder+ymis != length(x.mo)) {
    stop("somethings wrong i think")
  }

  perf <- nrow(perf)
  # whatever's left is where both tables have the same ranks named, but different names. this is a misclassification:
  moremis <- nrow(x) - (xmis+ymis+yover+yunder+perf)
  result <- data.frame(all.match = c(perf), mis = c(ymis+xmis+moremis), over = c(yover), under = c(yunder))
  if (sum(result) == nrow(x)) {
    return(result)
    # return(list(result, mis.i, yover.i, yunder.i, perf.i))
  } else {
    stop("noooooooooooooo")
  }
}

### this is the pr2 comparisons... scroll past plot saving for silva
# make a list of all taxonomy tables:
tbl.list <- list(idtax.pr2, bayes.pr2,
                 idtax.silva.map, bayes.silva.map)
all.comp <- lapply(tbl.list, FUN = tblcomper, x = etax.pr2)
all.comp <- base::do.call(base::rbind.data.frame, all.comp)
row.names(all.comp) <- c("idtax-pr2", "bayes-pr2", "idtax-silva", "bayes-silva")
all.comp <- all.comp / rowSums(all.comp)
all.comp$tbl <- row.names(all.comp)
plt.all.comp <- melt(all.comp)

# plot the results:
p.comp2.pr2 <- ggplot(plt.all.comp, aes(fill = tbl, x = variable, y = value)) +
  geom_bar(stat="identity", color = "black", position=position_dodge(width=0.8)) +
  labs(x = "Relative to ensemble-res", y = "Proportion of ASVs") +
  scale_x_discrete(breaks = c("all.match","mis","over","under"),
                   labels = c("Agree", "Misclassified", "Overclassified","Underclassified")) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0),
                     limits = c(0, 1), expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_text(size = 12, face="bold"),
        axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12, face="bold"),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "white"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "white"),
        axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_fill_manual(name = "Taxonomy",
                    breaks = c("idtax-pr2", "bayes-pr2", "idtax-silva", "bayes-silva", "ensemble-acc"),
                    values = cbPalette[c(1:5)])

print(p.comp2.pr2)
ggsave(paste0("~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/Figures/","compare_assignments_PR2.pdf"), p.comp2.pr2, device="pdf",width=11,height=8.5,units="in")

write.csv(seqtab.nochime,"~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/all18s_naames_exp/seqtab_nochime.csv",row.names=TRUE)
write.csv(etax.pr2,"~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/all18s_naames_exp/etax_pr2.csv",row.names=TRUE)
write.csv(etax.silva,"~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/all18s_naames_exp/etax_silva.csv",row.names=TRUE)
