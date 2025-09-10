# AIM ---------------------------------------------------------------------
# rerun the integration analysis on the full dataset after removing the cells flagged by Martina

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

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat5 workflow
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 1000 * 1024^2)

# read in the data --------------------------------------------------------
# load the LUT
# LUT <- read_csv("../../data/LUT_samples.csv")

# read in the full dataset.
# here I am loading the object from the previous analysis
sobj <- readRDS("../../out/object/manualClean/20_sobj_integrated_cleanup.rds")

DimPlot(sobj,group.by = "RNA_snn_res.0.4",label = T)

sobj@meta.data %>%
  filter(RNA_snn_res.0.4 == "19") %>%
  summarise(n = n())

# 195718 - 1201 = 194517

# add the filter varibale to the full object
sobj$barcode <- colnames(sobj)

# read in the metadata containing of all the subclusters
meta_vector <- dir("../../out/object/subclusters/") %>% str_subset("clean") %>% str_subset("23_")
meta_all_keep <- lapply(meta_vector, function(x){
  readRDS(paste0("../../out/object/subclusters/",x))
}) %>%
  bind_rows()

dim(meta_all_keep)

sobj@meta.data

# keep only the cells in the filtered metadata
sobj_filter <- subset(sobj,subset = barcode %in% rownames(meta_all_keep))

# create the seurat object ------------------------------------------------
# variable to bring over the new object
meta_keep <- sobj_filter@meta.data[,c("project_id","original_sample_name","sample_number","sample_id","nbb","autopsy","cohort","sex","age","braak","amyloid","braaklb","pmd","ph","weight","csf","apoe","barcode","id","tissuecode","iduit","datumuit","uitvraag","recipient","dcodeprot","dcode","diagnosis","dcodewk","wcode","region","location","pathological_stage","demyelination","demyelination.score","specific","ocode","storage","cell_id_old")]

# new AssayV5
# do not filter for quality
sobj_total <- CreateSeuratObject(counts = sobj_filter@assays$RNA@counts,
                                 project = "cleanup_V5",
                                 meta.data = meta_keep,
                                 min.cells = 0, min.features = 0) %>%
  # NOTICE: If I normalize before the cell cycle assignament I obtain a different call for the cell cycle identity.
  Seurat::NormalizeData(verbose = T)

# after creating the object I do not need the list and the merged matrix, free up some space
remove(sobj,sobj_filter,meta_keep,meta_all_keep)
gc()

# add the cell cycle analysis
DefaultAssay(sobj_total) <- "RNA"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
sobj_total <- CellCycleScoring(sobj_total, s.features = s.genes, g2m.features = g2m.genes)
sobj_total$percent.mt <- PercentageFeatureSet(sobj_total, pattern = "^MT-")
sobj_total$percent.ribo <- PercentageFeatureSet(sobj_total, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
# add also the percentage of globin. in this dataset it is not meaningful as there is no blood
sobj_total$percent.globin <- Seurat::PercentageFeatureSet(sobj_total,pattern = "^Hb[^(p)]")

# standard processing -----------------------------------------------------
# check the scale matrix
sobj_total@assays$RNA$scale.data

# standard processing
sobj_total <- sobj_total %>%
  # skip the normalizatio that has been already performed at the beginning
  # Seurat::NormalizeData(verbose = T) %>%
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
DimPlot(sobj_total,group.by = "orig.ident",raster = T)

# Run Harmony -------------------------------------------------------------
# The simplest way to run Harmony is to pass the Seurat object and specify which variable(s) to integrate out. RunHarmony returns a Seurat object, updated with the corrected Harmony coordinates. Let's set plot_convergence to TRUE, so we can make sure that the Harmony objective function gets better with each round.
sobj_total_h <- sobj_total %>%
  RunHarmony("orig.ident", plot_convergence = TRUE)

# Harmony with two or more covariates
# Do the same with your Seurat object:
# seuratObject <- RunHarmony(seuratObject, c("dataset", "donor", "batch_id"))
# To directly access the new Harmony embeddings, use the Embeddings command.
# harmony_embeddings <- Embeddings(sobj_total_h, 'harmony')
# harmony_embeddings[1:5, 1:5]
# Let's make sure that the datasets are well integrated in the first 2 dimensions after Harmony.
# DimPlot(object = sobj_total_h, reduction = "harmony", pt.size = .1, group.by = "sample_id")

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

# DimPlot(sobj_total_h,group.by = "ID",raster = F)

# save the object
saveRDS(sobj_total_h,"../../out/object/analysis_R44/24_sobj_integrated_cleanup.rds")