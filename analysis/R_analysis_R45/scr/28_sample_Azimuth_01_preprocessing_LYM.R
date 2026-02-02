# AIM ---------------------------------------------------------------------
# sample processing of Azimuth processing for reference labelling
# this part produce the file needed for the actual label transfer 
# this is for the LYM subset dataset

# libraries ---------------------------------------------------------------
library(Seurat)
library(Azimuth)
library(SeuratData)
library(patchwork)
library(tidyverse)
library(ComplexHeatmap)
library(cowplot)

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
# options(Seurat.object.assay.version = "v5")

# read in the dataset -----------------------------------------------------
# this is the query dataset
data.combined <- readRDS(file = "../../out/object/analysis_R44/27_LYM_subcluster_HarmonySample.rds")
DimPlot(data.combined)

# generate the diet object and save it
# save the coordiantes of the UMAP
data.combined@reductions$umap@cell.embeddings %>%
  data.frame() %>%
  rownames_to_column("barcode") %>%
  write_tsv("../../out/table/analysis_R45/28_LYM_subcluster_HarmonySample_coordUMAP.tsv")

DefaultAssay(data.combined) <- "RNA"
object_diet <- DietSeurat(object = data.combined, assays = "RNA")
saveRDS(object_diet,"../../out/object/analysis_R45/28_LYM_subcluster_HarmonySample_diet.rds")

# available datasets ------------------------------------------------------
# pick a reference
available_data <- AvailableData()

# install the picked reference dataset
# InstallData("pbmcref")

# if not working it is possible to install manually from the tar.gz file
# download.file(url = "http://seurat.nygenome.org/src/contrib/humancortexref.SeuratData_1.0.0.tar.gz",
#               destfile = "../data/azimuth/humancortexref.SeuratData_1.0.0.tar.gz")
# install.packages("../data/azimuth/humancortexref.SeuratData_1.0.0.tar.gz")

# load the reference dataset
reference <- LoadReference(path = "renv/library/linux-rocky-9.5/R-4.5/x86_64-conda-linux-gnu/pbmcref.SeuratData/azimuth/")
DimPlot(reference$plot,group.by = "celltype.l1")

# save the reference umap
df_point <- reference$map@reductions$refUMAP@cell.embeddings %>%
  data.frame() %>%
  rownames_to_column("barcode")
write_tsv(df_point,"../../out/table/analysis_R45/azimuth_pbmcref_CoordUMAP.tsv")

# save the meta of the ref
df_meta <- reference$map@meta.data %>%
  data.frame() %>%
  rownames_to_column("barcode")
write_tsv(df_meta,"../../out/table/analysis_R45/azimuth_pbmcref_Metadata.tsv")
