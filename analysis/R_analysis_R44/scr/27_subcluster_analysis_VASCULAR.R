# AIM ---------------------------------------------------------------------
# run the subcluster analysis for the VAS cells from 

# LIBRARIES ---------------------------------------------------------------
library(scater)
library(Seurat)
library(tidyverse)
library(robustbase)
# library(SeuratData)
library(patchwork)
library(harmony)
library(ggrepel)
library(ComplexHeatmap)

# read in the data --------------------------------------------------------
# read in the cells from run 02
sobj <- readRDS("../../out/object/analysis_R44/26_sobj_integrated_cleanup_manualAnnotation.rds")

DimPlot(sobj,raster = T,group.by = "cell_id",label = T)

# wrangling ---------------------------------------------------------------
# subset only the cells of interest.
# martina suggested to keep IMM and LYM all togheter for this task
sobj_subset <- subset(sobj,subset = cell_id %in% c("VAS"))
DimPlot(sobj_subset,group.by = "cell_id",label=T)

# define the minimal metadata
meta_minimal <- sobj_subset@meta.data %>% select(project_id:demyelination_short, -seurat_clusters,- contains("RNA_snn_res"))

# notice it is critacal that all matrices have the same dimension
# I need to create a single object to add the cell cycle scoring and other metadata. I decided not to trimm further the dataset for genes content
sobj_total <- CreateSeuratObject(counts = sobj_subset@assays$RNA$counts,
                                 project = "VAS_subset",
                                 meta.data = meta_minimal,
                                 min.cells = 0, min.features = 0) %>%
  # normalize before the cell cycle assignament
  Seurat::NormalizeData(verbose = T)

# after creating the object I do not need the list and the merged matrix, free up some space
remove("meta_minimal","sobj","sobj_subset")
gc()

# add the cell cycle analysis
DefaultAssay(sobj_total) <- "RNA"
# s.genes <- cc.genes$s.genes
# g2m.genes <- cc.genes$g2m.genes
# sobj_total <- CellCycleScoring(sobj_total, s.features = s.genes, g2m.features = g2m.genes)
# sobj_total$percent.mt <- PercentageFeatureSet(sobj_total, pattern = "^MT-")
# sobj_total$percent.ribo <- PercentageFeatureSet(sobj_total, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
# sobj_total$percent.globin <- Seurat::PercentageFeatureSet(sobj_total,pattern = "^HB[^(P)]")

# check the scale matrix
sobj_total@assays$RNA$scale.data
# pull all the genes to scale
# all.genes <- rownames(sobj_total)

# rescale the data for regressing out the sources of variation do not scale all the genes. if needed I can scale them before the heatmap call. for speeding up the computation I will keep 
sobj_total <- sobj_total %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 4000) %>%
  # I can scale the missing features afterwards now focus on the highly variable one for speed purposes
  ScaleData(vars.to.regress = c("percent.mt","nCount_RNA","S.Score","G2M.Score"), verbose = T) %>% 
  # run this if you want to scale all the variables
  # ScaleData(vars.to.regress = c("percent.mt.harmony","nCount_RNA.harmony","S.Score.harmony","G2M.Score.harmony"), verbose = T,features = all.genes) %>% 
  RunPCA(npcs = 30, verbose = T) %>% 
  RunUMAP(reduction = "pca", dims = 1:30,return.model = TRUE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = seq(0.1, 1, by = 0.1)) %>%
  identity()

# check the status of dataset preintegration
DimPlot(sobj_total,group.by = "sample_id",raster = T)

# Run Harmony -------------------------------------------------------------
# The simplest way to run Harmony is to pass the Seurat object and specify which variable(s) to integrate out. RunHarmony returns a Seurat object, updated with the corrected Harmony coordinates. Let's set plot_convergence to TRUE, so we can make sure that the Harmony objective function gets better with each round.
sobj_total_h <- sobj_total %>%
  RunHarmony("sample_id", plot_convergence = TRUE)

# Downstream analysis -----------------------------------------------------
# Many downstream analyses are performed on low dimensional embeddings, not gene expression. To use the corrected Harmony embeddings rather than PCs, set reduction = 'harmony'. For example, let's perform the UMAP and Nearest Neighbor analyses using the Harmony embeddings.
sobj_total_h <- sobj_total_h %>%
  RunUMAP(reduction = "harmony", dims = 1:30,return.model = TRUE) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = seq(0.1, 1, by = 0.1)) %>%
  # FindClusters(resolution = 0.5) %>%
  identity()

# verify that all the relevant slots are filled
sobj_total_h@assays$RNA$counts[1:20,1:10]
sobj_total_h@assays$RNA$data[1:20,1:10]
sobj_total_h@assays$RNA$scale.data[1:20,1:10]

dim(sobj_total_h@assays$RNA$counts)
dim(sobj_total_h@assays$RNA$data)
dim(sobj_total_h@assays$RNA$scale.data)

DimPlot(sobj_total_h,raster = T,group.by = "sample_id",label = T)
DimPlot(sobj_total_h,raster = T,group.by = "RNA_snn_res.0.4",label = T)

# sobj_total_h <- readRDS("../../out/object/data.combined_harmonySkipIntegration_AllSoupX_01000_6000_15.rds")
saveRDS(sobj_total_h,"../../out/object/analysis_R44/27_VAS_subcluster_HarmonySample.rds")
