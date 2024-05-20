# Spatial transcriptomic analysis
# Author: Millie Sander
# Part 1: data processing

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

#loading image data 
image_4337 <- Read10X_Image("/lustre/projects/Research_Project-191391/Project_10769/sample10769_4337/outs/spatial",image.name = "tissue_lowres_image.png",filter.matrix = TRUE,)
image_4414 <- Read10X_Image("/lustre/projects/Research_Project-191391/Project_10769/sample10769_4414/outs/spatial",image.name = "tissue_lowres_image.png",filter.matrix = TRUE,)
image_4415 <- Read10X_Image("/lustre/projects/Research_Project-191391/Project_10769/sample10769_4415/outs/spatial",image.name = "tissue_lowres_image.png",filter.matrix = TRUE,)
image_4416 <- Read10X_Image("/lustre/projects/Research_Project-191391/Project_10769/sample10769_4416/outs/spatial",image.name = "tissue_lowres_image.png",filter.matrix = TRUE,)
image_4417 <- Read10X_Image("/lustre/projects/Research_Project-191391/Project_10769/sample10769_4417/outs/spatial",image.name = "tissue_lowres_image.png",filter.matrix = TRUE,)
image_4419 <- Read10X_Image("/lustre/projects/Research_Project-191391/Project_10769/sample10769_4419/outs/spatial", image.name = "tissue_lowres_image.png", filter.matrix = TRUE,)
image_4420 <- Read10X_Image("/lustre/projects/Research_Project-191391/Project_10769/sample10769_4420/outs/spatial",image.name = "tissue_lowres_image.png",filter.matrix = TRUE,)
image_4466 <- Read10X_Image("/lustre/projects/Research_Project-191391/Project_10769/sample10769_4466/outs/spatial",image.name = "tissue_lowres_image.png",filter.matrix = TRUE,)
image_4467 <- Read10X_Image("/lustre/projects/Research_Project-191391/Project_10769/sample10769_4467/outs/spatial",image.name = "tissue_lowres_image.png",filter.matrix = TRUE,)
image_4472 <- Read10X_Image("/lustre/projects/Research_Project-191391/Project_10769/sample10769_4472/outs/spatial",image.name = "tissue_lowres_image.png",filter.matrix = TRUE,)
image_8221 <- Read10X_Image("/lustre/projects/Research_Project-191391/Project_10769/sample10769_8221/outs/spatial",image.name = "tissue_lowres_image.png",filter.matrix = TRUE,)

# loading sample data, example shown for sample 4416 and repeated per sample 

sample_4416 <- Load10X_Spatial(
  "/lustre/projects/Research_Project-191391/Project_10769/sample10769_4416/outs",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "image_4416",
  filter.matrix = TRUE,
  to.upper = FALSE,
  image = image_4416)


# setting sample info 
Idents(object = sample_4337) <- "sample_4337"
Idents(object = sample_4414) <- "sample_4414"
Idents(object = sample_4415) <- "sample_4415"
Idents(object = sample_4416) <- "sample_4416"
Idents(object = sample_4417) <- "sample_4417"
Idents(object = sample_4419) <- "sample_4419"
Idents(object = sample_4420) <- "sample_4420"
Idents(object = sample_4466) <- "sample_4466"
Idents(object = sample_4467) <- "sample_4467"
Idents(object = sample_4472) <- "sample_4472"
Idents(object = sample_8221) <- "sample_8221"

sample_4337$type = "HET"
sample_4414$type = "HET"
sample_4415$type = "HET"
sample_4416$type = "WT"
sample_4417$type = "HET"
sample_4419$type = "HET"
sample_4420$type = "WT"
sample_4466$type = "WT"
sample_4467$type = "WT"
sample_4472$type = "WT"
sample_8221$type = "HET" 

sample_4337$slide = "1"
sample_4414$slide = "1"
sample_4415$slide = "1"
sample_4416$slide = "1"
sample_4417$slide = "2"
sample_4419$slide = "2"
sample_4420$slide = "2"
sample_4466$slide = "2"
sample_4467$slide = "3"
sample_4472$slide = "3"
sample_8221$slide = "3" 

sample_4337$sex = "F"
sample_4414$sex = "M"
sample_4415$sex = "M"
sample_4416$sex = "M"
sample_4417$sex = "F"
sample_4419$sex = "F"
sample_4420$sex = "F"
sample_4466$sex = "M"
sample_4467$sex = "M"
sample_4472$sex = "M"
sample_8221$sex = "F" 

sample_4337$orig.ident = as.factor("sample_4337")
sample_4414$orig.ident = as.factor("sample_4414")
sample_4415$orig.ident = as.factor("sample_4415")
sample_4416$orig.ident = as.factor("sample_4416")
sample_4417$orig.ident = as.factor("sample_4417")
sample_4419$orig.ident = as.factor("sample_4419")
sample_4420$orig.ident = as.factor("sample_4420")
sample_4466$orig.ident = as.factor("sample_4466")
sample_4467$orig.ident = as.factor("sample_4467")
sample_4472$orig.ident = as.factor("sample_4472")
sample_8221$orig.ident = as.factor("sample_8221")


# merging and filtering datasets 
merged <- merge(sample_4337, sample_4414)
merged <- merge(merged, sample_4415)
merged <- merge(merged, sample_4416)
merged <- merge(merged, sample_4417)
merged <- merge(merged, sample_4419)
merged <- merge(merged, sample_4420)
merged <- merge(merged, sample_4466)
merged <- merge(merged, sample_4467)
merged <- merge(merged, sample_4472)
merged <- merge(merged, sample_8221)

merged <- PercentageFeatureSet(merged, "^mt-", col.name = "percent_mito")
merged <- PercentageFeatureSet(merged, "^Hb.*-", col.name = "percent_hb")
merged <- PercentageFeatureSet(merged, "^Rp.*-", col.name = "percent_ribo")


# violin plot of QC measures 
pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/nCounts_new.pdf", height = 8, width = 12)
VlnPlot(merged, features = c("nCount_Spatial"), pt.size = 0.01, ncol = 2) + NoLegend() 
dev.off()

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/nFeatures.pdf", height = 8, width = 12)
VlnPlot(merged, features = c("nFeature_Spatial"), pt.size = 0.1, ncol = 2) + NoLegend()
dev.off()

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/percent_mito.pdf", height = 8, width = 12)
VlnPlot(merged, features = c("percent_mito"), pt.size = 0.1, ncol = 2) + NoLegend()
dev.off()

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/percent_hb.pdf", height = 8, width = 12)
VlnPlot(merged, features = c("percent_hb"), pt.size = 0.1, ncol = 2) + NoLegend()
dev.off

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/percent_ribo.pdf", height = 8, width = 12)
VlnPlot(merged, features = c("percent_ribo"), pt.size = 0.1, ncol = 2) + NoLegend()
dev.off()


# spatial heatmap of QC measures
pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/nCount_spatial.pdf", height = 8, width = 14)
SpatialFeaturePlot(merged, features = c("nCount_Spatial"), keep.scale = "all", pt.size.factor = 8, ncol = 4)
dev.off()

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/nFeature_spatial.pdf", height = 12, width = 20)
SpatialFeaturePlot(merged, features = c("nFeature_Spatial"), keep.scale = "all", pt.size.factor = 8, ncol = 4)
dev.off()

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/percent_mito_spatial.pdf", height = 12, width = 20)
SpatialFeaturePlot(merged, features = c("percent_mito"), keep.scale = "all", pt.size.factor = 8, ncol = 4)
dev.off()

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/percent_hb_spatial.pdf", height = 12, width = 20)
SpatialFeaturePlot(merged, features = c("percent_hb"), keep.scale = "all", pt.size.factor = 8, ncol = 4)
dev.off()

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/percent_ribo_spatial.pdf", height = 12, width = 20)
SpatialFeaturePlot(merged, features = c("percent_ribo"), keep.scale = "all", pt.size.factor = 8, ncol = 4)
dev.off()


# visualising top expressed genes for potential contamination
A = merged@assays$Spatial@layers
A@x = A@x/rep.int(colSums(A), diff(A@p))
most_expressed <- order(Matrix::rowSums(A), decreasing = T)[20:1]
most_expressed <- as.matrix(most_expressed)
A <- as.matrix(A)

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/most_expressed_A.pdf", height = 6, width = 8)
boxplot(as.matrix(t(A[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell",
    col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)
dev.off() 


# UMAP plot of original data
pdf(file = "/lustre/home/ms739/spatial_transcriptomics/final_analysis/plots/cluster_dimplot.pdf", height = 6, width = 10)
DimPlot(merged, group.by = "orig.ident")
dev.off()

# data filtering 
data.filt <- merged[, merged$nFeature_Spatial > 200]
selected_f <- rownames(data.filt)[Matrix::rowSums(data.filt) > 3]
data.filt <- subset(data.filt, features = selected_f)
data.filt <- data.filt[!grepl("Gm42418", rownames(data.filt)), ]

data.filt <- ScaleData(data.filt) # vars.to.regress = "percent_mito"
data.filt <- FindVariableFeatures(data.filt, selection.method = "vst")
data.filt <- RunPCA(data.filt)
data.filt <- RunUMAP(data.filt, dims = 1:30, reduction = "pca")
data.filt <- FindNeighbors(data.filt, reduction = "pca", dims = 1:30)
data.filt <- FindVariableFeatures(data.filt, selection.method = "vst")

# UMAP after filtering
pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/dimPlot_filtered.pdf", height = 6, width = 8)
DimPlot(data.filt, group.by = "orig.ident")
dev.off() 

# splitting filtered data in order to normalise (SCTransform) individually 
split <- SplitObject(data.filt, split.by = "orig.ident")
# removed irrelevant images from each sample, e.g. to remove image 4414 from sample 4337 data:
# filt_4337@images[["image_4414"]] = NULL 
# repeat for all 

# normalising each sample individually using SCTransform
filtered <- list(filt_4337, filt_4414, filt_4415, filt_4416, filt_4417, filt_4419, filt_4420, filt_4466, filt_4467, filt_4472, filt_8221)

for (i in 1:length(filtered)) {
    filtered[[i]] <- SCTransform(filtered[[i]], vst.flavor = "v2", assay = "Spatial", verbose = FALSE)
}
features <- SelectIntegrationFeatures(object.list = filtered, nfeatures = 3000)

filtered_merged <- merge(filtered[[1]], 
y = filtered[2:length(filtered)],  
project = "j20", 
merge.data = TRUE)

DefaultAssay(filtered_merged) <- "SCT"
VariableFeatures(filtered_merged) <- features

filtered_merged <- ScaleData(filtered_merged)
filtered_merged <- RunPCA(object = filtered_merged, assay = "SCT", npcs = 20)

# determine the significance of PCs ("dims") (for downstream analysis)
elbow_plot <- ElbowPlot(filtered_merged)
pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/elbow_plot.pdf", height = 6, width = 8)
elbow_plot
dev.off()

# harmony batch correction on individually normalised samples
filtered_merged <- RunHarmony(object = filtered_merged,
                                    assay.use = "SCT",
                                    reduction = "pca",
                                    dims.use = 1:20,
                                    group.by.vars = c("slide", "sex"), 
                                    plot_convergence = TRUE)

filtered_merged <- RunUMAP(object = filtered_merged, assay = "SCT", reduction = "harmony", dims = 1:20)

# UMAP after batch correction 
pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/harmony_corrected.pdf", height = 6, width = 8)
DimPlot(filtered_merged, group.by = "orig.ident")
dev.off()


