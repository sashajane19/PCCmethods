# Sasha adapting Dylan code
# this script adds maps the cleaned up taxonomy table onto the compilation
# of functional assignments for particular lineages of protists.
setwd("~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/tax2fcn/")

tax <- read.csv("naames_exp_etax_20211026_noNA.csv",header = TRUE)
trophs <- read.csv("merged_trophic_assignments_20210810.csv",
                   header = TRUE, stringsAsFactors = FALSE)

# extract col of most-resolved ASV:
lowest <- apply(tax, MARGIN = 1, FUN = function(x) max(which(!is.na(x))))
tax.troph <- cbind(tax, data.frame(matrix(NA, nrow = nrow(tax), ncol = 2)))
colnames(tax.troph) <- c(colnames(tax), c("phyto", "troph"))
no.troph <- vector(mode = "character")
tmp.no.troph <- vector(mode = "character")
for (i in 1:length(lowest)) {
  col <- lowest[i]
  taxnam <- tax[i, col]
  while (col > 3 && !(taxnam %in% trophs$Tax_name)) {
    col <- col-1
    tmp.no.troph <- c(tmp.no.troph, taxnam)
    taxnam <- tax[i, col]
  }
  # so now if col == 3 the lineage will be unassigned
  if (col == 3) {
    # lineage is unassigned...
    # that means do nothing since they're already NA...
    # but do save the names:
    no.troph <- c(no.troph, tmp.no.troph)
  } else if (taxnam %in% trophs$Tax_name) {
    # taxnam was found in trophs. so extract the trophs and assign
    row.match <- which(trophs$Tax_name %in% taxnam)
    troph.match <- trophs[row.match , ]

    # assign pigmented or not:
    tax.troph$phyto[i] <- troph.match$pigmented

    if (all(is.na(troph.match[, c("Troph2", "Troph3", "Troph4", "Exceptions1", "Exceptions2",
                                  "Exceptions3", "Exceptions4", "Exceptions5", "Exceptions6",
                                  "Exceptions7", "Exceptions8", "Exceptions9", "Exceptions10")]))) {
      # assign the Troph column as it's not ambiguous
      tax.troph$troph[i] <- troph.match$Troph

    } else {
      # assign NA to Troph column (so just do nothing...)
    }
  }

  tmp.no.troph <- vector(mode = "character")
}

no.troph <- unique(no.troph)

# proportion of ASVs known phyto or known troph
sum(is.na(tax.troph$phyto))/nrow(tax.troph)
sum(is.na(tax.troph$troph))/nrow(tax.troph)
saveRDS(tax.troph, file = "naames_exp_eTax_w_troph_20211026.rds")
write.csv(tax.troph, file = "naames_exp_eTax_w_troph_20211026.csv")
write.csv(no.troph, file = "naames_exp_no_troph_20211026.csv")
utax.troph <- unique(tax.troph[, 3:ncol(tax.troph)], MARGIN = 1)
write.csv(utax.troph, file = "naames_exp_unique_taxNtroph_20211026.csv")

# read in your otu table and sample data...
otu <- read.csv("naames_exp_etax_otu_noNA.csv",header = TRUE)
# get ASVs w/ known "phyto" or not
sv.known.phyto <- tax.troph$svN[tax.troph$phyto %in% c(1) | tax.troph$phyto %in% c(0)]
# hist of rel abundance of known phyto or not across all samples:
hist(rowSums(otu[, sv.known.phyto]))

sv.known.troph <- tax.troph$sv[!(is.na(tax.troph$troph))]
hist(rowSums(otu[, sv.known.troph]))

unktax.troph1 <- unique(tax.troph[is.na(tax.troph$phyto) , 3:ncol(tax.troph)], MARGIN = 1) # where phyto is unknown
unktax.troph2 <- unique(tax.troph[is.na(tax.troph$troph) , 3:ncol(tax.troph)], MARGIN = 1) # where troph is unknown
write.csv(unktax.troph1, file = "unknown_phyto_lin_20211026.csv")
write.csv(unktax.troph2, file = "unknown_troph_lin_20211026.csv")

