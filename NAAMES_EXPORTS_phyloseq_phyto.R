# Sasha Kramer
# UCSB IGPMS
# 20211130

# NAAMES + EXPORTS 18S with phyloseq (just phyto samples)
# Tutorial from: https://vaulot.github.io/tutorials/Phyloseq_tutorial.html

# Install packages
if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
BiocManager::install("phyloseq")

install.packages("dplyr")     # To manipulate dataframes
install.packages("readxl")    # To read Excel files into R
install.packages("ggplot2")   # for high quality graphics

library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names

otu_phy <- read_excel("~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/phyloseq/NA_EXP_phylo_phyto.xlsx", sheet = "OTU")
tax_phy <- read_excel("~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/phyloseq/NA_EXP_phylo_phyto.xlsx", sheet = "taxo")
samples_df <- read_excel("~/Documents/UCSB/Research/Data/NAAMES/Genomics/18S/phyloseq/NA_EXP_phylo_phyto.xlsx", sheet = "samples")

otu_phy <- otu_phy %>%
  tibble::column_to_rownames("otu")

tax_phy <- tax_phy %>%
  tibble::column_to_rownames("otu")

samples_df <- samples_df %>%
  tibble::column_to_rownames("sample")

otu_phy <- as.matrix(otu_phy)
tax_phy <- as.matrix(tax_phy)

OTUphy = otu_table(otu_phy, taxa_are_rows = TRUE)
TAXphy = tax_table(tax_phy)
samples = sample_data(samples_df)

phylo_phy <- phyloseq(OTUphy, TAXphy, samples)
phylo_phy

total = median(sample_sums(phylo_phy))
standf = function(x, t=total) round(t * (x / sum(x)))
phylo_phy = transform_sample_counts(phylo_phy, standf)

# plot all samples based on division
plot_div <- plot_bar(phylo_phy, fill = "division") +
geom_bar(aes(color=division, fill=division), stat="identity", position="stack")
ggsave(paste0("~/Documents/UCSB/Research/Figures/HPLC_18S_16S/20211130/","plot_div_phyto.pdf"), plot_div, device="pdf",width=13,height=7,units="in")

# regroup NAAMES + EXPORTS samples
phylo_cruise <- merge_samples(phylo_phy, "cruise")
plot_bar(phylo_cruise, fill = "division") +
  geom_bar(aes(color=division, fill=division), stat="identity", position="stack")

# try a heatmap
plot_heatmap(phylo_phy, method = "NMDS", distance = "bray")

# too crowded
phylo_abund <- filter_taxa(phylo_phy, function(x) sum(x > total*0.05) > 0, TRUE)
phylo_abund

plot_heatmap(phylo_abund, method = "NMDS", distance = "bray")

heatm <- plot_heatmap(phylo_abund, method = "MDS", distance = "(A+B-2*J)/(A+B-J)",
             taxa.label = "class", taxa.order = "class",
             trans=NULL, low="beige", high="red", na.value="beige")
ggsave(paste0("~/Documents/UCSB/Research/Figures/HPLC_18S_16S/20211130/","heatmap_phyto.pdf"), heatm, device="pdf",width=13,height=7,units="in")

# look at diversity
adiv <- plot_richness(phylo_phy, measures=c("Chao1", "Shannon"))
ggsave(paste0("~/Documents/UCSB/Research/Figures/HPLC_18S_16S/20211130/","chao_shannon.pdf"), adiv, device="pdf",width=13,height=7,units="in")
adiv_c <- plot_richness(phylo_phy, measures=c("Chao1", "Shannon"), x="cruise", color = "lat")
ggsave(paste0("~/Documents/UCSB/Research/Figures/HPLC_18S_16S/20211130/","chao_shannon_cruise.pdf"), adiv_c, device="pdf",width=13,height=7,units="in")

# ordination
phylo.ord <- ordinate(phylo_phy, "NMDS", "bray")
ord_all <- plot_ordination(phylo_phy, phylo.ord, type="taxa", color="class",
                title="ASVs", label="class") +
  facet_wrap(~division, 2)
ggsave(paste0("~/Documents/UCSB/Research/Figures/HPLC_18S_16S/20211130/","ordination_all.pdf"), ord_all, device="pdf",width=13,height=7,units="in")

ord_c <- plot_ordination(phylo_phy, phylo.ord, type="samples", color="lat",
                shape="cruise", title="Samples") + geom_point(size=3)
ggsave(paste0("~/Documents/UCSB/Research/Figures/HPLC_18S_16S/20211130/","ordination_cruise.pdf"), ord_c, device="pdf",width=9,height=7,units="in")

ord_ct <- plot_ordination(phylo_phy, phylo.ord, type="split", color="class",
                shape="cruise", title="biplot") +
  geom_point(size=3)
ggsave(paste0("~/Documents/UCSB/Research/Figures/HPLC_18S_16S/20211130/","ordination_cruise_tax.pdf"), ord_ct, device="pdf",width=13,height=7,units="in")

# networks
plot_net(phylo_phy, distance = "(A+B-2*J)/(A+B)", type = "taxa",
         maxdist = 0.7, color="class", point_label="genus")

plot_net(phylo_abund, distance = "(A+B-2*J)/(A+B)", type = "taxa",
         maxdist = 1, color="class", point_label="genus")

