# Spatial transcriptomic analysis
# Author: Millie Sander
# Part 3: cell type deconvolution  

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

# hippocampal SC reference: 
# set the directory
refdir = "/lustre/home/ms739/spatial_transcriptomics/thesis/CARD"
# set the path
reference_path <- file.path(refdir, "SCRef_hippocampus.RDS")
# load reference
hippocampus_reference <- readRDS(reference_path)

# create count matrix for SC data: rows = genes; columns = spot (barcode)
SC_counts <- hippocampus_reference@assays$RNA@counts

meta <- hippocampus_reference@meta.data
meta_data <- subset(meta, select = -c(nCount_RNA, nFeature_RNA, cluster, nUMI))
# left with "orig.ident" (sample info) and "liger_ident_coarse" (cell type)
# change order
meta_data <- meta_data[,c("liger_ident_coarse", "orig.ident")]
# rename
colnames(meta_data) <- c("cellType", "sampleInfo")

# extract counts and coordinates per sample
# example shown using sample 4416, repeated per sample 
counts_4416 <- sample_4416@assays$Spatial@counts
coords_4416 <- GetTissueCoordinates(sample_4416)
colnames(coords_4416) <- c("x", "y")

# create CARD object 
# example shown using sample 4416, repeat per sample 
CARD_obj_4416 = createCARDObject(
  sc_count = SC_counts,
  sc_meta = meta_data,
  spatial_count = counts_4416,
  spatial_location = coords_4416,
  ct.varname = "cellType",
  ct.select = NULL,
  sample.varname = "sampleInfo",
  minCountGene = 100,
  minCountSpot = 5) 

# deconvolute the spatial data using the new CARD object
CARD_obj_4337 = CARD_deconvolution(CARD_object = CARD_obj_4337)
CARD_obj_4414 = CARD_deconvolution(CARD_object = CARD_obj_4414)
CARD_obj_4415 = CARD_deconvolution(CARD_object = CARD_obj_4415)
CARD_obj_4416 = CARD_deconvolution(CARD_object = CARD_obj_4416)
CARD_obj_4417 = CARD_deconvolution(CARD_object = CARD_obj_4417)
CARD_obj_4419 = CARD_deconvolution(CARD_object = CARD_obj_4419)
CARD_obj_4420 = CARD_deconvolution(CARD_object = CARD_obj_4420)
CARD_obj_4466 = CARD_deconvolution(CARD_object = CARD_obj_4466)
CARD_obj_4467 = CARD_deconvolution(CARD_object = CARD_obj_4467)
CARD_obj_4472 = CARD_deconvolution(CARD_object = CARD_obj_4472)
CARD_obj_8221 = CARD_deconvolution(CARD_object = CARD_obj_8221)

# visualising the cell proportions for each spot
colors = c("#FFD92F","#4DAF4A","#FCCDE5","#D9D9D9","#377EB8","#7FC97F","#BEAED4",
    "#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02",
    "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")

# example shown using sample 4416, repeated per sample: 
pie_4416 <- CARD.visualize.pie(
  proportion = CARD_obj_4416@Proportion_CARD,
  spatial_location = CARD_obj_4416@spatial_location, 
  colors = colors, 
    radius = 0.7)

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/CARD/plots/pie_4416.pdf", height = 5, width = 7)
pie_4416
dev.off()

## select specific cell types that we are interested
ct.visualize = c("Astrocyte", "CA1", "CA3", "Denate", "Interneuron", "Microglia_Macrophages", "Oligodendrocyte", "Polydendrocyte")

# example using 4416, repeated per sample
cell_types_4416 <- CARD.visualize.prop(
  proportion = CARD_obj_4416@Proportion_CARD,        
  spatial_location = CARD_obj_4416@spatial_location, 
  ct.visualize = ct.visualize,
  colors = c("lightblue","lightyellow","red"), 
  NumCols = 4,                                 
        pointSize = 3.0)  

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/CARD/plots/cell_plots_4416.pdf", height = 10, width = 15)
cell_types_4416
dev.off()

# to increase spatial resolution, example using 4416, repeated per sample 
CARD_obj_4416 = CARD.imputation(CARD_obj_4416,NumGrids = 2000,ineibor = 10,exclude = NULL)

location_imputation_4416 = cbind.data.frame(x=as.numeric(sapply(strsplit(rownames(CARD_obj_4416@refined_prop),split="x"),"[",1)),
  y=as.numeric(sapply(strsplit(rownames(CARD_obj_4416@refined_prop),split="x"),"[",2)))
rownames(location_imputation_4416) = rownames(CARD_obj_4416@refined_prop)

p5_4416 <- CARD.visualize.prop(
  proportion = CARD_obj_4416@refined_prop,                         
  spatial_location = location_imputation_4416,            
  ct.visualize = ct.visualize,                    
  colors = c("lightblue","lightyellow","red"),    
  NumCols = 4)                                  

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/CARD/plots/high_res_4416.pdf", height = 10, width = 15)
p5_4416
dev.off()
