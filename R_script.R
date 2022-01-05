
library(microbiome) 
library(phyloseq)
library(ggplot2)
library(RColorBrewer) 
library(ggpubr) 
library(DT)
library(data.table)
library(dplyr) 
library(ape)
library(vegan)
library(pheatmap)

### upload data 
Phylo <- read_phyloseq(otu.file = "Biom_file.biom1", 
                            taxonomy.file = NULL, 
                            metadata.file = "Metadata.csv", 
                            type = "biom")
treefile <- read.tree("Tree_file.tree")
Phylo <- merge_phyloseq(Phylo, treefile)

### subset project data
Phylo_1 <- subset_samples(Phylo,  ProjectName == "Fe_SO4_Andrea_V4_new") 
Phylo_1 <- prune_taxa(taxa_sums(Phylo_1) > 0, Phylo_1)
Phylo_1  ## 2894 taxa 62 samples 

### data check
datatable(tax_table(Phylo_1))

Phylo_2 <- subset_taxa(Phylo_1, Domain !="NA")
Phylo_2 <- subset_taxa(Phylo_2, Family !="f__Mitochondria")
Phylo_2 <- subset_taxa(Phylo_2, Order != "o__Chloroplast")
Phylo_2 <- subset_taxa(Phylo_2, Domain !="k__Archaea")
Phylo_2  ## 2826; 62

rarecurve(t(otu_table(Phylo_2)), step=50, cex=0.5)

summarize_phyloseq(Phylo_2)
nrRead <- sort(sample_sums(Phylo_2))
nrRead

### subset controls 
ps3 <- subset_samples(Phylo_2,  Seq_ID == "MOCK3" | Seq_ID == "MOCK4" | Seq_ID == "water") 
ps3 <- prune_taxa(taxa_sums(ps3) > 0, ps3)
ps3

### Figure PcOA
### subset per experiment 
ps4 <- subset_samples(Phylo_2,  Experiment  == "4") 
ps4<- prune_taxa(taxa_sums(ps4) > 0, ps4)
ps4

ord <- ordinate(ps4, method = "PCoA", distance = "unifrac")  ### "jaccard", binary = TRUE
p <- plot_ordination(ps4, ord, 
                     color = "Material",
                     shape = "Day")
p1 <- p + theme_bw() +
  scale_colour_manual(values = c("salmon1", "deepskyblue", "gray30"))
p1

### Figure Heatmap
ps12 <- subset_samples(Phylo_2,  Experiment    == "1"  | Experiment    == "2") 
ps12 <- prune_taxa(taxa_sums(ps12) > 0, ps12)
ps12    # 1499 - 25

Trans1 <- aggregate_taxa(ps12, 'Genus')  ### 519 taxa
Cl1_Com <- microbiome::transform(Trans1, "compositional")
Cl1Cor <- core(Cl1_Com, detection = 0.01, prevalence = 10/100)
Cl1Cor

otu <- abundances(Cl1Cor)
metadata <- meta(Cl1Cor)

all.samples <- intersect(rownames(metadata),colnames(otu))
my.metadata <- metadata

data <- metadata

otu <- otu[,all.samples]

my.metadata <- my.metadata[all.samples,]

myd <- otu[,rownames(my.metadata)]
myd <- myd[rowSums(myd)!=0,]

G1 = rownames(data[data$Experiment == "1",])
G2 = rownames(data[data$Experiment == "2",]) 

fcs <- c()
for (tax in rownames(myd)) {
  fcs[[tax]] <- mean(myd[tax, as.character(G1)]) - mean(myd[tax, G2])
}


pvals <- check_wilcoxon(myd,G1 ,G2 , sort = T, p.adjust.method = "BH", paired=FALSE)
r <- names(which(pvals < 0.05))
r <- names(which(pvals < 2))

str(r)
taxa <- r 

############################### FIGURE
meta <- meta(Cl1)
head(meta)
str(meta)
col_groups <- dplyr::select(meta, Experiment, Day, Material)  

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 80)
d <- distance(Cl1, method = "unifrac")   ## select distance
hpws <- hclust(d, method = "ward.D2")

# Reorder Intervention  levels
meta$Experiment = factor(meta$Experiment, levels = c("1", "2")) 
IntCol <- c("steelblue", "chartreuse3")
names(IntCol) <- levels(meta$Experiment)

# Reorder Day COLOM 1
meta$Day <- factor(meta$Day, levels = c("311", "218", "219", "583", "582","269", "599", "775"))
AgeCol <- c("burlywood4", "cadetblue2","aquamarine3" ,"slategray" ,"pink2","thistle3", "darkseagreen2", "tan1")
names(AgeCol) <- levels(meta$Day)

# Reorder Day COLUM 2
#meta$Day <- factor(meta$Day, levels = c("322", "219", "220", "592", "277","615", "795"))
#AgeCol <- c("burlywood4", "cadetblue2","aquamarine3" ,"slategray" ,"pink2","thistle3", "darkseagreen2", "tan1")
#names(AgeCol) <- levels(meta$Day)

# Reorder Material 
meta$Material <- factor(meta$Material, levels = c("sediment", "inoculum", "liquid"))
HouseCol <- c("pink", "blue", "green")
names(HouseCol) <- levels(meta$Material)

AnnColour <- list(
  Experiment = IntCol,
  Day = AgeCol,
  Marterial = HouseCol) 

p <- pheatmap(myd[taxa,], scale = 'row',color = my_palette,cluster_cols = hpws, cutree_cols = 6, 
              cluster_rows = TRUE, annotation_col = col_groups, annotation_colors = AnnColour)



