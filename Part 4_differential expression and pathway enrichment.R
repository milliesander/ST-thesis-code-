# Spatial transcriptomic analysis
# Author: Millie Sander
# Part 4: differential expression analysis  

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Seurat")
BiocManager::install("ggplot2")
BiocManager::install("patchwork")
BiocManager::install("dplyr")
BiocManager::install("sctransform")
BiocManager::install("harmony")
BiocManager::install("RColorBrewer")
BiocManager::install("writexl")
BiocManager::install("clusterProfiler")
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("Matrix", version = 1.5)
BiocManager::install("DESeq2")
BiocManager::install("glmGamPoi")
BiocManager::install("ggvenn")
BiocManager::install("MetBrewer")
BiocManager::install("enrichR")
BiocManager::install("GeneOverlap")
BiocManager::install("devtools")
BiocManager::install("tidyverse")
BiocManager::install("rrvgo")

devtools::install_github('YingMa0107/CARD')
devtools::install_github('xuranw/MuSiC')

# load required packages
library(Seurat)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(ggplot2)
library(patchwork)
library(dplyr)
library(sctransform)
library(harmony)
library(RColorBrewer)
library(glmGamPoi)
library(ggvenn)
library(MetBrewer)
library(enrichR)
library(GeneOverlap)
library(devtools)
library(CARD)
library(WGCNA)
library(hdWGCNA)
library(cowplot)
library(patchwork)
library(igraph)

Idents(hippocampus) <- "type"

# first upon the entire hippocampus
hippocampus_DEGs <- FindMarkers(hippocampus, assay = "SCT", ident.1 = "HET", ident.2 = "WT", verbose = FALSE, recorrect_umi = FALSE, logfc.threshold = 0) 
write.xlsx(hippocampus_DEGs, file = "/lustre/home/ms739/spatial_transcriptomics/thesis/hippocampus_DEGs.xlsx")

# filtering for significant (p < 0.05) DEGs
hippocampus_DEGs <- hippocampus_DEGs[which(hippocampus_DEGs$p_val < 0.05),]
write.xlsx(hippocampus_DEGs, file = "/lustre/home/ms739/spatial_transcriptomics/thesis/p_hippocampus_DEGs.xlsx")

# next, individual clusters 
# example shown on cluster 41 
CA1_pyrmaidal_DEGs <- FindMarkers(hippocampus, assay = "SCT", ident.1 = "41_HET", ident.2 = "41_WT", verbose = FALSE, recorrect_umi = FALSE, logfc.threshold = 0)
write.xlsx(CA1_pyrmaidal_DEGs, file = "/lustre/home/ms739/spatial_transcriptomics/thesis/CA1_pyramidal_DEGs.xlsx")

CA1_pyrmaidal_DEGs <- CA1_pyrmaidal_DEGs[which(CA1_pyrmaidal_DEGs$p_val < 0.05),]
write.xlsx(CA1_pyrmaidal_DEGs, file = "/lustre/home/ms739/spatial_transcriptomics/thesis/p_CA1_pyramidal_DEGs.xlsx")

# repeat for each other cluster (6, 23, 26, 35, 36)

# volcano plots 
# entire hippocampus 

hippocampus_col <- hippocampus_DEGs

hippocampus_col <- hippocampus_col %>% 
  mutate(
    Expression = case_when(avg_log2FC >= 0 & p_val_adj <= 0.05 ~ "p adj upregulated",
              avg_log2FC <= 0 & p_val_adj <= 0.05 ~ "p adj downregulated", 
              avg_log2FC >= 0 & p_val <= 0.05 ~ "upregulated",
                           avg_log2FC <= 0 & p_val <= 0.05 ~ "downregulated",
                           TRUE ~ "Unchanged"))

hipp_p1 <- ggplot(hippocampus_col, aes(avg_log2FC, -log(p_val,10))) +
  geom_point(aes(color = Expression), size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"p")) +
  scale_color_manual(values = c("gray50", "#6495ed", "dodgerblue3", "firebrick3", "#ff69b4")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 

# adding App gene label 

hipp_app <- hippocampus_p_adj[1,]
hipp_app$symbol <- rownames(hipp_app)

hipp_p2 <- hipp_p1 + geom_label_repel(data = hipp_app,
                   mapping = aes(avg_log2FC, -log(p_val,10), label = symbol),
                   size = 2)

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/hipp_gg_volcano.pdf", height = 4, width = 5)
hipp_p2
dev.off() 

# repeated for each cluster independently
# example shown using cluster 41 (CA1 pyramidal cell layer): 
Clust41_col <- CA1_pyramidal_DEGs
Clust41_col <- Clust41_col %>% 
  mutate(
    Expression = case_when(avg_log2FC >= 0 & p_val_adj <= 0.05 ~ "p adj upregulated",
              avg_log2FC <= 0 & p_val_adj <= 0.05 ~ "p adj downregulated", 
              avg_log2FC >= 0 & p_val <= 0.05 ~ "upregulated",
                           avg_log2FC <= 0 & p_val <= 0.05 ~ "downregulated",
                           TRUE ~ "Unchanged"))

Clust41_p1 <- ggplot(Clust41_col, aes(avg_log2FC, -log(p_val,10))) +
  geom_point(aes(color = Expression), size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"p")) +
  scale_color_manual(values = c("gray50", "#6495ed", "firebrick3", "#ff69b4", "dodgerblue3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 

# adding FDR significant gene labels 
Clust41_FDR <- Clust41_col[1:2,]
Clust41_FDR$symbol <- rownames(Clust41_FDR)

Clust41_p2 <- Clust41_p1 + geom_label_repel(data = Clust41_FDR,
                   mapping = aes(avg_log2FC, -log(p_val,10), label = symbol),
                   size = 2)

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/Clust41_gg_volcano.pdf", height = 4, width = 5)
Clust41_p2
dev.off() 

# plotting expression of APP per cluster, between genotypes
# extracting gene counts per barcode, with genotype, cluster, and sample info 
App <- FetchData(hippocampus, vars = c("App", "ident", "orig.ident"), layer = "counts", clean = "none") 

# seperate based on cluster 
App_cluster41 <- subset(App, ident == "41")
App_cluster26 <- subset(App, ident == "26")
App_cluster36 <- subset(App, ident == "36")
App_cluster23 <- subset(App, ident == "23")
App_cluster35 <- subset(App, ident == "35")
App_cluster6 <- subset(App, ident == "6")

# remove names of clusters so we can average 
App_cluster41 <- subset(App_cluster41, select = -c(ident))
App_cluster26 <- subset(App_cluster26, select = -c(ident))
App_cluster36 <- subset(App_cluster36, select = -c(ident))
App_cluster23 <- subset(App_cluster23, select = -c(ident))
App_cluster35 <- subset(App_cluster35, select = -c(ident))
App_cluster6 <- subset(App_cluster6, select = -c(ident))

# average expression by sample
App_avg_41 <- aggregate(.~orig.ident, data = App_cluster41, mean)
App_avg_26 <- aggregate(.~orig.ident, data = App_cluster26, mean)
App_avg_36 <- aggregate(.~orig.ident, data = App_cluster36, mean)
App_avg_23 <- aggregate(.~orig.ident, data = App_cluster23, mean)
App_avg_35 <- aggregate(.~orig.ident, data = App_cluster35, mean)
App_avg_6 <- aggregate(.~orig.ident, data = App_cluster6, mean)

# export and combine with genotype info 
write.xlsx(App_avg_41, file = "/lustre/home/ms739/spatial_transcriptomics/thesis/counts/App_avg_41.xlsx")
write.xlsx(App_avg_26, file = "/lustre/home/ms739/spatial_transcriptomics/thesis/counts/App_avg_26.xlsx")
write.xlsx(App_avg_36, file = "/lustre/home/ms739/spatial_transcriptomics/thesis/counts/App_avg_36.xlsx")
write.xlsx(App_avg_23, file = "/lustre/home/ms739/spatial_transcriptomics/thesis/counts/App_avg_23.xlsx")
write.xlsx(App_avg_35, file = "/lustre/home/ms739/spatial_transcriptomics/thesis/counts/App_avg_35.xlsx")
write.xlsx(App_avg_6, file = "/lustre/home/ms739/spatial_transcriptomics/thesis/counts/App_avg_6.xlsx")

# plotted using external software (OriginPro)
# creating a venn diagram of overlapping/cluster specific DEGs (nominally significant)
DEGs <- list("Cluster 41" = rownames(CA1_pyrmaidal_DEGs),
    "Cluster 26" = rownames(CA3_pyrmaidal_DEGs), 
    "Cluster 36" = rownames(DG_DEGs), 
    "Cluster 23" = rownames(cluster35_DEGs), 
    "Cluster 35" = rownames(cluster23_DEGs), 
    "Cluster 6" = rownames(cluster6_DEGs)
    )

custom_colour <- c("#0073C2FF", "#EFC000FF", "#CD534CFF", "#CD135CFF","#00C2FF","#009E73")

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/venn.pdf", height = 16, width = 16)
venn(DEGs, ilabels = F, zcol = c("#0073C2FF", "#EFC000FF", "#CD534CFF", "#CD135CFF","#00C2FF","#009E73"),
     box = F,ilcs=2, sncs = 3, ellipse=T)
dev.off()

# identifying FDR corrected DEGs 
FDR_CA1_pyramidal <- subset(CA1_pyrmaidal_DEGs, select = -c(p_val, avg_log2FC, pct.1, pct.2, p_val_adj))
FDR_CA3_pyramidal <- subset(CA3_pyrmaidal_DEGs, select = -c(p_val, avg_log2FC, pct.1, pct.2, p_val_adj))
FDR_DG_granule <- subset(DG_DEGs, select = -c(p_val, avg_log2FC, pct.1, pct.2, p_val_adj))
FDR_clust35 <- subset(cluster35_DEGs, select = -c(p_val, avg_log2FC, pct.1, pct.2, p_val_adj))
FDR_clust6 <- subset(cluster6_DEGs, select = -c(p_val, avg_log2FC, pct.1, pct.2, p_val_adj))
FDR_clust23 <- subset(cluster23_DEGs, select = -c(p_val, avg_log2FC, pct.1, pct.2, p_val_adj))

# creating a venn diagram of overlapping/cluster specific DEGs (nominally significant)
FDR_DEGs <- list("Cluster 41" = rownames(FDR_CA1_pyramidal),
    "Cluster 26" = rownames(FDR_CA3_pyramidal), 
    "Cluster 36" = rownames(FDR_DG_granule), 
    "Cluster 23" = rownames(FDR_clust35), 
    "Cluster 35" = rownames(FDR_clust23), 
    "Cluster 6" = rownames(FDR_clust6)
    )

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/FDR_venn.pdf", height = 16, width = 16)
venn(FDR_DEGs, ilabels = F, zcol = c("#0073C2FF", "#EFC000FF", "#CD534CFF", "#CD135CFF","#00C2FF","#009E73"),
     box = F,ilcs=2, sncs = 3, ellipse=T)
dev.off()

# to actually identify FDR corrected DEGs specific to clusters: 
colnames(FDR_CA1_pyramidal) <- "CA1_pyr"
colnames(FDR_CA3_pyramidal) <- "CA3_pyr"
colnames(FDR_DG_granule) <- "DG_granule"
colnames(FDR_clust35) <- "clust35"
colnames(FDR_clust6) <- "clust6"
colnames(FDR_clust23) <- "clust23"

FDR_CA1_pyramidal$gene <- rownames(FDR_CA1_pyramidal)
FDR_CA3_pyramidal$gene <- rownames(FDR_CA3_pyramidal)
FDR_DG_granule$gene <- rownames(FDR_DG_granule)
FDR_clust35$gene <- rownames(FDR_clust35)
FDR_clust6$gene <- rownames(FDR_clust6)
FDR_clust23$gene <- rownames(FDR_clust23)

FDR_merged <- merge(FDR_CA1_pyramidal, FDR_CA3_pyramidal, by = "gene", all.x = TRUE, all.y = TRUE)
FDR_merged <- merge(FDR_merged, FDR_DG_granule, by = "gene", all.x = TRUE, all.y = TRUE)
FDR_merged <- merge(FDR_merged, FDR_clust23, by = "gene", all.x = TRUE, all.y = TRUE)
FDR_merged <- merge(FDR_merged, FDR_clust35, by = "gene", all.x = TRUE, all.y = TRUE)
FDR_merged <- merge(FDR_merged, FDR_clust6, by = "gene", all.x = TRUE, all.y = TRUE)
# did not include hippocampus as a whole 

rownames(FDR_merged) <- FDR_merged$gene
FDR_merged <- subset(FDR_merged, select = -c(gene))

# example shown for CA1 pyramidal cell layer (cluster 41), repeated for each cluster: 
CA1pyr_specific <- subset(FDR_merged, CA1_pyr != "NA") 
CA1pyr_specific <- CA1pyr_specific[is.na(CA1pyr_specific$CA3_pyr),]
CA1pyr_specific <- CA1pyr_specific[is.na(CA1pyr_specific$DG_granule),]
CA1pyr_specific <- CA1pyr_specific[is.na(CA1pyr_specific$clust35),]
CA1pyr_specific <- CA1pyr_specific[is.na(CA1pyr_specific$clust23),]
CA1pyr_specific <- CA1pyr_specific[is.na(CA1pyr_specific$clust6),]
write.xlsx(CA1pyr_specific, file = "/lustre/home/ms739/spatial_transcriptomics/thesis/CA1pyr_specific.xlsx")

# pathway enrichment analysis 
# conducted upon the hippocampal data and each cluster independently
# example shown using entire hippocampal data

GO_hipp <- rownames(hippocampus_p_adj)
GO_hipp_entrez <- mapIds(org.Mm.eg.db, keys = GO_hipp, keytype = "SYMBOL", column = "ENTREZID")
GO_hipp_results <- enrichGO(gene = GO_hipp_entrez,
                       OrgDb = "org.Mm.eg.db",
                       ont = "BP")

# network construction of terms
edox_hipp <- setReadable(GO_hipp_results, 'org.Mm.eg.db', 'ENTREZID')

# checking how many terms are returned
GO_terms_hipp <- as.list(edox_hipp@result$ID)
GO_terms_hipp <- unique(GO_terms_hipp)
length(GO_terms_hipp)
# 6377

# reduce terms using rrvgo package

# to create the similarity matrix of GO terms: 
simMatrix_hipp <- calculateSimMatrix(edox_hipp$ID,
                                orgdb="org.Mm.eg.db",
                                ont="BP",
                                method="Rel")

scores_hipp <- setNames(-log10(edox_hipp$qvalue), edox_hipp$ID)

reducedTerms_hipp <- reduceSimMatrix(simMatrix_hipp,
                                scores_hipp,
                                threshold=0.9,
                                orgdb="org.Mm.eg.db")
# threshold 0.8 = 87
# threshold 0.9 = 34

# visualising parent terms 
# take the name and size of the parent terms and merge with the original enrichment dataframe
# from this, take the qvalue of the original term 
# then will plot the parent term, against it's original qvalue, and the size of the spots will be the size of the parent term group 

reducedTerms_hipp$parentTerm <- as.factor(reducedTerms_hipp$parentTerm)
parentSize <- list(summary(reducedTerms_hipp$parentTerm))
parentSize <- as.data.frame(parentSize)
colnames(parentSize) <- c("parentSize")
parentSize$term <- rownames(parentSize)

# next merge the dataframes so we have the size of each parent term contained within the reduced term dataframe 
termSig <- data.frame(edox_hipp$Description, edox_hipp$qvalue)
colnames(termSig) <- c("term", "qvalue")

hipp_parent <- merge(termSig, parentSize, by = "term", all = FALSE)

# ordering by qvalue and group size 
hipp_parent <- hipp_parent[order(hipp_parent$qvalue),]

# selecting the top 10
top_hipp <- head(hipp_parent, 10)

top10_hippocampus <- ggplot(top_hipp, aes(x = as.factor(term), y = -log10(qvalue) ,size = parentSize))+
  geom_point(fill = "#ffa500",colour = "#674200", shape = 21, )+
  coord_flip()+
  ggtitle("Hippocampus")+
  theme_bw()

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/FDR/top10_hipp_reduced.pdf", height = 3, width = 7)
top10_hippocampus
dev.off()

# repeat for each cluster individually 
