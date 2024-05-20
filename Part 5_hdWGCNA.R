# Spatial transcriptomic analysis
# Author: Millie Sander
# Part 5: hdWGCNA

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

# enable parallel processing for network analysis (optional)
enableWGCNAThreads(nThreads = 8)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# create meta data slot contain the coordinates for each spot 
DefaultAssay(hippocampus) <- "SCT"

image_df <- do.call(rbind, lapply(names(hippocampus@images), function(x){
  hippocampus@images[[x]]@coordinates
}))

# merge the image_df with the Seurat metadata
new_meta <- merge(hippocampus@meta.data, image_df, by='row.names')

# fix the row ordering to match the original seurat object
rownames(new_meta) <- new_meta$Row.names
ix <- match(as.character(colnames(hippocampus)), as.character(rownames(new_meta)))
new_meta <- new_meta[ix,]

# add the new metadata to the seurat object
hippocampus@meta.data <- new_meta
head(image_df)

hippocampus <- ScaleData(hippocampus)

hippocampus <- SetupForWGCNA(
  hippocampus,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "practise"
)

# after this point, can not subset data!

# creating metaspots, reduces sparcity of the matrix (zero entries)
# metaspots will be created for cells of the same group - important to specify to group by sample, and cluster 
hippocampus <- MetaspotsByGroups(
  hippocampus,
  group.by = c("orig.ident", "seurat_clusters"),
  ident.group = "seurat_clusters",
  assay = "SCT",
  slot = "counts",
  wgcna_name = "practise"
)

hippocampus  <- NormalizeMetacells(hippocampus)

# to extract the metaspot object: 
m_obj <- GetMetacellObject(hippocampus, wgcna_name = "practise")
m_obj

hippocampus  <- SetDatExpr(
  hippocampus,
  group.by=NULL, # if NULL, will use the seurat Idents as groups  
  group_name = NULL, 
  use_metacells=TRUE, 
  assay = "Spatial", 
  slot = "data"
)

hippocampus <- TestSoftPowers(hippocampus, networkType = "signed")
plot_list <- PlotSoftPowers(hippocampus)

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/hdWGCNA/hipp_soft_thresholding.pdf", height = 9, width = 15)
wrap_plots(plot_list, ncol=2)
dev.off() 

hippocampus <- ConstructNetwork(
  hippocampus,
  networkType = "signed", 
  TOMType = "signed",
  tom_name='test',
  overwrite_tom=TRUE, 
  soft_power = 12, 
  min.module.size = 50, 
  detectCutHeight = 0.995, 
  mergeCutHeight = 0.2 
)

# changing module colours (optional)

# get the module table
modules <- GetModules(hippocampus)
mods <- unique(modules$module)

# make a table of the module-color pairings
mod_colors_df <- dplyr::select(modules, c(module, color)) %>%
  distinct %>% arrange(module)
rownames(mod_colors_df) <- mod_colors_df$module

# get a table of just the module and it's unique color
mod_color_df <- GetModules(hippocampus) %>%
  dplyr::select(c(module, color)) %>%
  distinct %>% arrange(module)

# the number of unique modules (subtract 1 because the grey module stays grey):
n_mods <- nrow(mod_color_df) - 1

# using the "Signac" palette from metbrewer, selecting for the number of modules
new_colors <- c("#7ba5f0", "#b380df", "#ffc968", "#f08080", "#97eb7a")

# reset the module colors
hippocampus <- ResetModuleColors(hippocampus, new_colors)

hippocampus <- ModuleEigengenes(hippocampus, group.by.vars = "orig.ident")
hippocampus <- ModuleConnectivity(hippocampus, sparse = TRUE)

# renaming modules with prefix "SM" for "spatial modules"
hippocampus <- ResetModuleNames(
  hippocampus,
  new_name = "SM"
)
modules <- GetModules(hippocampus) %>% subset(module != 'grey')
head(modules[,1:3])

# get module eigengenes and gene-module assignment tables
MEs <- GetMEs(hippocampus)
modules <- GetModules(hippocampus)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# add the MEs to the seurat metadata so we can plot it with Seurat functions
hippocampus@meta.data <- cbind(hippocampus@meta.data, MEs)

# plot with Seurat's DotPlot function
p <- DotPlot(hippocampus, features=mods, group.by = 'seurat_clusters', dot.min=0.1)

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  xlab('') + ylab('')

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/hdWGCNA/module_dotplot.pdf", height = 6, width = 10)
p
dev.off()

# plotting modules spatially, upon sample 4416 as example 
sample_mods <- SpatialFeaturePlot(
  hippocampus,
  features = c(SM2_mod, SM3_mod, SM5_mod),
  alpha = c(0, 1),
  ncol = 5, 
  images = "image_4416",
  pt.size.factor = 16, keep.scale = "all", 
  crop = TRUE
) & theme(aspect.ratio = ratio_4416) & ggplot2::scale_fill_gradient2(limits = c(0, 15), mid = "#ffff80", high = "#e6004c")


pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/hdWGCNA/sample_mods_4416.pdf", height = 6, width = 10)
sample_mods
dev.off()  

# to plot modules with the correct aspect ratio
# upon image 4416 as an example: 

SM1_mod <- levels(modules$module); SM1_mod <- SM1_mod[SM1_mod == 'SM1']
SM2_mod <- levels(modules$module); SM2_mod <- SM2_mod[SM2_mod == 'SM2']
SM3_mod <- levels(modules$module); SM3_mod <- SM3_mod[SM3_mod == 'SM3']
SM4_mod <- levels(modules$module); SM4_mod <- SM4_mod[SM4_mod == 'SM4']
SM5_mod <- levels(modules$module); SM5_mod <- SM5_mod[SM5_mod == 'SM5']

ratio_4416 <- (max(hippocampus@images$image_4416@coordinates$imagerow) - min(hippocampus@images$image_4416@coordinates$imagerow)) / (max(hippocampus@images$image_4416@coordinates$imagecol) - min(hippocampus@images$image_4416@coordinates$imagecol))

SM1_4416 <- SpatialFeaturePlot(hippocampus, features = SM1_mod, alpha = c(0.1, 1), ncol = 5, images = "image_4416", pt.size.factor = 16, crop = TRUE) + theme(aspect.ratio = ratio_4416)
SM2_4416 <- SpatialFeaturePlot(hippocampus, features = SM2_mod, alpha = c(0.1, 1), ncol = 5, images = "image_4416", pt.size.factor = 16, crop = TRUE) + theme(aspect.ratio = ratio_4416)
SM3_4416 <- SpatialFeaturePlot(hippocampus, features = SM3_mod, alpha = c(0.1, 1), ncol = 5, images = "image_4416", pt.size.factor = 16, crop = TRUE) + theme(aspect.ratio = ratio_4416)
SM4_4416 <- SpatialFeaturePlot(hippocampus, features = SM4_mod, alpha = c(0.1, 1), ncol = 5, images = "image_4416", pt.size.factor = 16, crop = TRUE) + theme(aspect.ratio = ratio_4416)
SM5_4416 <- SpatialFeaturePlot(hippocampus, features = SM5_mod, alpha = c(0.1, 1), ncol = 5, images = "image_4416", pt.size.factor = 16, crop = TRUE) + theme(aspect.ratio = ratio_4416)

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/hdWGCNA/SM1_4416.pdf", height = 6, width = 10)
SM1_4416
dev.off()

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/hdWGCNA/SM2_4416.pdf", height = 6, width = 10)
SM2_4416
dev.off()

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/hdWGCNA/SM3_4416.pdf", height = 6, width = 10)
SM3_4416
dev.off()

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/hdWGCNA/SM4_4416.pdf", height = 6, width = 10)
SM4_4416
dev.off()

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/hdWGCNA/SM5_4416.pdf", height = 6, width = 10)
SM5_4416
dev.off()

# network plots within individual folders: 
ModuleNetworkPlot(hippocampus, 
  mods = "all", 
  outdir = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/hdWGCNA/ModuleNetworks") 
# this will create a new folder in the (WinSCP) path, with a PDF image file for each module 

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/hdWGCNA/HubGene_network_plot.pdf", height = 4, width = 4)
HubGeneNetworkPlot(
  hippocampus,
  n_hubs = 3, n_other=5,
  edge_prop = 0.75,
  mods = 'all'
)
dev.off()

# differential module eigengene (DME) expression analysis 
# hdWGCNA function "FindDMEs" is a special case of the Seurat function "FindMarkers"
# first need to extract barcodes for each group to be tested (HET and WT)
group1_HET <- hippocampus@meta.data %>% subset(type == "HET") %>% rownames
group2_WT <- hippocampus@meta.data %>% subset(type != "HET") %>% rownames

hippocampus_DMEs <- FindDMEs(
  hippocampus,
  barcodes1 = group1_HET,
  barcodes2 = group2_WT,
  test.use='wilcox',
  wgcna_name='practise'
)

lollyPlot <- PlotDMEsLollipop(
  hippocampus, 
  hippocampus_DMEs, 
  wgcna_name='practise', 
  pvalue = "p_val_adj", 
  cex = 4
)

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/hdWGCNA/lollyPlot.pdf", height = 3, width = 4)
lollyPlot
dev.off()

# pathway enrichment analysis 

dbs <- c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023')

hippocampus <- RunEnrichr(
  hippocampus,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 100, # number of genes per module to test
  wgcna_name = "practise", 
  wait_time = 60
)
# increased wait time to 60 (seconds per module) so that EnrichR website treats each module independently and doesnt overlap output 

# retrieve the output table
enrich_df <- GetEnrichrTable(hippocampus)

# EnrichrBarPlot will create an output file containing pdfs for each module and GO result 
EnrichrBarPlot(
  hippocampus,
  outdir = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/hdWGCNA/enrichr_plots", # name of output directory
  n_terms = 10, # number of enriched terms to show (sometimes more show if there are ties!)
  plot_size = c(6,8), # width, height of the output .pdfs
  logscale=TRUE # do you want to show the enrichment as a log scale?
)

## plotting the top term in each module and seeing overlap 
EnrichRDotPlot <- EnrichrDotPlot(
  hippocampus,
  mods = "all", # use all modules (this is the default behavior)
  database = "GO_Biological_Process_2023", # this has to be one of the lists we used above!!!
  n_terms=1 # number of terms for each module
)

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/plots/hdWGCNA/enrichrDotPlot.pdf", height = 3, width = 5.5)
EnrichRDotPlot
dev.off()

# module preservation analysis 

# hippocampal SC reference: 
refdir = "/lustre/home/ms739/spatial_transcriptomics/thesis/CARD"
reference_path <- file.path(refdir, "SCRef_hippocampus.RDS")
hippocampus_reference <- readRDS(reference_path)

hippocampus_reference <- FindVariableFeatures(hippocampus_reference, selection.method = "vst", nfeatures = 2000)
hippocampus_reference <- ScaleData(hippocampus_reference)
hippocampus_reference <- RunPCA(hippocampus_reference)
hippocampus_reference <- RunUMAP(hippocampus_reference, dims = 1:30, reduction = "pca")

hippocampus_reference <- ProjectModules(
  seurat_obj = hippocampus_reference,
  seurat_ref = hippocampus,
  modules = modules,
  wgcna_name_proj="hdWGCNA_projected", 
  wgcna_name = "tutorial", 
  overlap_proportion = 0.5)

hippocampus_reference <- ModuleExprScore(
  hippocampus_reference,
  n_genes = 25,
  method='Seurat'
)

hippocampus_reference <- ModuleConnectivity(hippocampus_reference, assay="RNA", slot="data")

# saving the query Seurat object after creating/projecting modules 
saveRDS(hippocampus_reference, file=paste0("/lustre/home/ms739/spatial_transcriptomics/thesis/module_preservation", 'HippRef_ModPreservation.rds'))

# to extract the module eigengenes for downstream analyses
projected_hMEs <- GetModules(hippocampus_reference)

# visualisation
# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  hippocampus_reference,
  wgcna_name = "hdWGCNA_projected",
  features='hMEs',
  order=TRUE
)

# stitch together with patchwork
pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/module_preservation/projectedModule_plots.pdf", height = 10, width = 8)
wrap_plots(plot_list, ncol=3)
dev.off()

# set expression matrix for reference dataset
hippocampus <- SetDatExpr(
  hippocampus,
)

# set expression matrix for query dataset:
hippocampus_reference <- SetDatExpr(
  hippocampus_reference,
)

# run module preservation analysis 
# this step takes overnight - submit via slurm Rscript in Linux 
hippocampus_reference <- ModulePreservation(
  hippocampus_reference,
  seurat_ref = hippocampus,
  name="ModPreservation",
  verbose=3,
  n_permutations=250
)

# getthe module preservation table
mod_pres <- GetModulePreservation(hippocampus_reference, "ModPreservation")$Z
obs <- GetModulePreservation(hippocampus_reference, "module_preservation")$obs

grep('summary', colnames(mod_pres))
grep('summary', colnames(obs))

plot_list <- PlotModulePreservation(
  hippocampus_reference,
  name="ModPreservation",
  statistics = "summary"
)

pdf(file = "/lustre/home/ms739/spatial_transcriptomics/thesis/module_preservation/summaryPlots.pdf", height = 4, width = 6)
wrap_plots(plot_list, ncol=2)
dev.off()
 
