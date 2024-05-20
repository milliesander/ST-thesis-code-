# Spatial transcriptomic analysis
# Author: Millie Sander
# Part 2: clustering 

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
library(rrvgo)

# creating clusters using a spatial resolution of 2.5
clustered <- FindNeighbors(object = filtered_merged, assay = "SCT", reduction = "harmony", dims = 1:20,
 nn.method = "annoy", compute.SNN = TRUE)

clustered_2.5 <- FindClusters(object = clustered, resolution = 2.5, graph.name = "SCT_snn")

# UMAP of clusters split by genotype 
pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/DimPlot_2.5_genotype.pdf", height = 6, width = 12)
DimPlot(clustered_2.5, reduction = "umap", split.by = "type")
dev.off()

# subset hippocampal clusters 
Idents(clustered_2.5) <- "seurat_clusters"

# prepare data for differential expression analysis prior to subsetting 
clustered_2.5 <- PrepSCTFindMarkers(clustered_2.5)

hippocampus <- subset(clustered_2.5, idents = c(6, 23, 26, 35, 36, 41))

# plot hippocampal clusters on original images 
pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/hipp.pdf", height = 12, width = 20)
SpatialDimPlot(clustered, cells.highlight = CellsByIdentities(hippocampus), pt.size.factor = 8, ncol = 4,
  cols.highlight = c("#f08080", "#839e2e", "#21ccbb", "#d899e5", "#a9a9a9"))
dev.off()

# subsetting to remove any random outlier spots 

hippocampus <- subset(hippocampus, image4416_imagecol < 70, invert = TRUE)

hippocampus <- subset(hippocampus, image4420_imagecol < 95, invert = TRUE)
hippocampus <- subset(hippocampus, image4420_imagecol > 120, invert = TRUE)
hippocampus <- subset(hippocampus, image4420_imagerow < 360, invert = TRUE)
hippocampus <- subset(hippocampus, image4420_imagerow > 400, invert = TRUE)

hippocampus <- subset(hippocampus, image4417_imagerow < 75, invert = TRUE)
hippocampus <- subset(hippocampus, image4417_imagecol > 125, invert = TRUE)

hippocampus <- subset(hippocampus, image4419_imagerow < 220, invert = TRUE)
hippocampus <- subset(hippocampus, image4419_imagerow > 265, invert = TRUE)
hippocampus <- subset(hippocampus, image4419_imagecol < 90, invert = TRUE)
hippocampus <- subset(hippocampus, image4419_imagecol > 120, invert = TRUE)

hippocampus <- subset(hippocampus, image4466_imagerow < 490, invert = TRUE)
hippocampus <- subset(hippocampus, image4466_imagecol > 115, invert = TRUE)

hippocampus <- subset(hippocampus, image4467_imagecol < 90, invert = TRUE)

hippocampus <- subset(hippocampus, image4472_imagecol < 80, invert = TRUE)

hippocampus <- subset(hippocampus, image8221_imagecol < 100, invert = TRUE)
hippocampus <- subset(hippocampus, image8221_imagecol > 140, invert = TRUE)

# visualise final data subset
pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/hipp_subset.pdf", height = 8, width = 16)
SpatialDimPlot(hippocampus, pt.size.factor = 10, ncol = 3)
dev.off()

# plotting hippocampal clusters per sample, anatomically correct (Seurat compresses to a square)
ratio_4414 <- (max(hippocampus@images$image_4414@coordinates$imagerow) - min(hippocampus@images$image_4414@coordinates$imagerow)) / (max(hippocampus@images$image_4414@coordinates$imagecol) - min(hippocampus@images$image_4414@coordinates$imagecol))
ratio_4415 <- (max(hippocampus@images$image_4415@coordinates$imagerow) - min(hippocampus@images$image_4415@coordinates$imagerow)) / (max(hippocampus@images$image_4415@coordinates$imagecol) - min(hippocampus@images$image_4415@coordinates$imagecol))
ratio_4416 <- (max(hippocampus@images$image_4416@coordinates$imagerow) - min(hippocampus@images$image_4416@coordinates$imagerow)) / (max(hippocampus@images$image_4416@coordinates$imagecol) - min(hippocampus@images$image_4416@coordinates$imagecol))
ratio_4417 <- (max(hippocampus@images$image_4417@coordinates$imagerow) - min(hippocampus@images$image_4417@coordinates$imagerow)) / (max(hippocampus@images$image_4417@coordinates$imagecol) - min(hippocampus@images$image_4417@coordinates$imagecol))
ratio_4419 <- (max(hippocampus@images$image_4419@coordinates$imagerow) - min(hippocampus@images$image_4419@coordinates$imagerow)) / (max(hippocampus@images$image_4419@coordinates$imagecol) - min(hippocampus@images$image_4419@coordinates$imagecol))
ratio_4420 <- (max(hippocampus@images$image_4420@coordinates$imagerow) - min(hippocampus@images$image_4420@coordinates$imagerow)) / (max(hippocampus@images$image_4420@coordinates$imagecol) - min(hippocampus@images$image_4420@coordinates$imagecol))
ratio_4466 <- (max(hippocampus@images$image_4466@coordinates$imagerow) - min(hippocampus@images$image_4466@coordinates$imagerow)) / (max(hippocampus@images$image_4466@coordinates$imagecol) - min(hippocampus@images$image_4466@coordinates$imagecol))
ratio_4467 <- (max(hippocampus@images$image_4467@coordinates$imagerow) - min(hippocampus@images$image_4467@coordinates$imagerow)) / (max(hippocampus@images$image_4467@coordinates$imagecol) - min(hippocampus@images$image_4467@coordinates$imagecol))
ratio_4472 <- (max(hippocampus@images$image_4472@coordinates$imagerow) - min(hippocampus@images$image_4472@coordinates$imagerow)) / (max(hippocampus@images$image_4472@coordinates$imagecol) - min(hippocampus@images$image_4472@coordinates$imagecol))
ratio_8221 <- (max(hippocampus@images$image_8221@coordinates$imagerow) - min(hippocampus@images$image_8221@coordinates$imagerow)) / (max(hippocampus@images$image_8221@coordinates$imagecol) - min(hippocampus@images$image_8221@coordinates$imagecol))
ratio_4337 <- (max(hippocampus@images$image_4337@coordinates$imagerow) - min(hippocampus@images$image_4337@coordinates$imagerow)) / (max(hippocampus@images$image_4337@coordinates$imagecol) - min(hippocampus@images$image_4337@coordinates$imagecol))

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/hipp_4414.pdf")
SpatialDimPlot(hippocampus, images = "image_4414", pt.size.factor = 10) + theme(aspect.ratio = ratio_4414)
dev.off()

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/hipp_4415.pdf", height = 8, width = 16)
SpatialDimPlot(hippocampus, images = "image_4415", pt.size.factor = 20) + theme(aspect.ratio = ratio_4415)
dev.off()

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/hipp_4416.pdf", height = 8, width = 16)
SpatialDimPlot(hippocampus, images = "image_4416", pt.size.factor = 26) + theme(aspect.ratio = ratio_4416)
dev.off()

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/hipp_4417.pdf", height = 8, width = 16)
SpatialDimPlot(hippocampus, images = "image_4417", pt.size.factor = 26) + theme(aspect.ratio = ratio_4417)
dev.off()

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/hipp_4419.pdf", height = 8, width = 16)
SpatialDimPlot(hippocampus, images = "image_4419", pt.size.factor = 24) + theme(aspect.ratio = ratio_4419)
dev.off()

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/hipp_4420.pdf", height = 8, width = 16)
SpatialDimPlot(hippocampus, images = "image_4420", pt.size.factor = 25) + theme(aspect.ratio = ratio_4420)
dev.off()

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/hipp_4466.pdf", height = 8, width = 16)
SpatialDimPlot(hippocampus, images = "image_4466", pt.size.factor = 8) + theme(aspect.ratio = ratio_4466)
dev.off()

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/hipp_4467.pdf", height = 8, width = 16)
SpatialDimPlot(hippocampus, images = "image_4467", pt.size.factor = 26) + theme(aspect.ratio = ratio_4467)
dev.off()

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/hipp_4472.pdf", height = 8, width = 16)
SpatialDimPlot(hippocampus, images = "image_4472", pt.size.factor = 10) + theme(aspect.ratio = ratio_4472)
dev.off()

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/hipp_8221.pdf", height = 8, width = 16)
SpatialDimPlot(hippocampus, images = "image_8221", pt.size.factor = 30) + theme(aspect.ratio = ratio_8221)
dev.off()

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/hipp_4337.pdf", height = 8, width = 16)
SpatialDimPlot(hippocampus, images = "image_4337", pt.size.factor = 10) + theme(aspect.ratio = ratio_4337)
dev.off()
