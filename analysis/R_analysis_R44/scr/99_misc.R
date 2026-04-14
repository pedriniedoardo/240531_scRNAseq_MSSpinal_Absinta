# AIM ---------------------------------------------------------------------
# misc task onthe object

# libraries ---------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(SeuratWrappers)
library(limma)
library(ggrepel)
library(presto)

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
# options(Seurat.object.assay.version = "v5")

# read in the data --------------------------------------------------------
data.combined <- readRDS("../../out/object/analysis_R44/26_sobj_integrated_cleanup_manualAnnotation.rds")

# confirm the identity of the object
# Idents(data.combined)
DimPlot(data.combined,label = T,order = T,group.by = "cell_id")
# DimPlot(data.combined,label = T,order = T)

# save the metadata -------------------------------------------------------
# save the full metadata of the obkect
df_meta <- data.combined@meta.data %>%
  rownames_to_column("barcodes")

df_meta %>%
  write_tsv("../../out/table/analysis_R44/99_meta_sobj_integrated_cleanup_manualAnnotation.tsv")
