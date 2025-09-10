# AIM ---------------------------------------------------------------------
# run a statistical test on the proportion differences across clusters, between conditions

# libraries ---------------------------------------------------------------
library(Seurat)
library(SeuratData)
library(tidyverse)
library(cacoa)
library(conos)
library(sccore)
library(coda.base)
library(psych)

# build the cacoa object from seurat --------------------------------------
# load the seurat object
so <- readRDS("../../out/object/data.combined_harmonySkipIntegration_AllSoupX_00500_07000_05_AnnotationSCType.rds")
DimPlot(so,label = T,group.by = "RNA_snn_res.0.1",raster=T)

# wrangling ---------------------------------------------------------------
meta <- so@meta.data %>%
  rownames_to_column()

# run cacoa ---------------------------------------------------------------
# Cacoa currently supports inputs in several formats (see below). Most of them require the following metadata:
#
# sample.groups: vector with condition labels per sample named with sample ids
sample.groups <- meta$diagnosis
names(sample.groups) <- meta$orig.ident

# cell.groups: cell type annotation vector named by cell ids
cell.groups <- meta$RNA_snn_res.0.1
names(cell.groups) <- meta$rowname

# sample.per.cell: vector with sample labels per cell named with cell ids
sample.per.cell <- meta$orig.ident
names(sample.per.cell) <- meta$rowname

# ref.level: id of the condition, corresponding to the reference (i.e. control)
table(sample.groups)
ref.level <- "Non-demented control"

# target.level: id of the condition, corresponding to the target (i.e. case)
target.level <- "Multiple sclerosis"

# Additionally, embedding parameter containing a matrix or data.frame with a cell embedding can be provided. Rownames should match to the cell ids. It is used for visualization and some cluster-free analysis.
embedding <- so@reductions$umap@cell.embeddings %>%
  data.frame()

# Parameter graph.name is required for cluster-free analysis, and must contain a name of joint graph in Seurat object. For that, the Seurat object must have a joint graph estimated (see FindNeighbors). For visualization purposes, Seurat also must have cell embedding estimated or the embedding data frame must be provided in the embedding parameter.
names(so@graphs)
graph.name <- "RNA_nn"

# Seurat object so
cao <- Cacoa$new(so,
                 sample.groups=sample.groups,
                 cell.groups=cell.groups,
                 sample.per.cell=sample.per.cell,
                 ref.level=ref.level,
                 target.level=target.level,
                 embedding = embedding,
                 graph.name=graph.name)

# run cacoa ---------------------------------------------------------------

# Cluster-based changes ---------------------------------------------------
# The fastest way to visualize changes in the dataset is to show, what cell types changed their abundance and expression patterns.

# Estimate cluster-based changes
cao$estimateCellLoadings()
cao$estimateExpressionShiftMagnitudes()

# Plot compositional changes
cao$plotCellLoadings(show.pvals=FALSE)
# save plot in case all the sample are used
ggsave("../../out/plot/cacoa_diagnosis.pdf",width = 10,height = 6)

# The red line here shows statistical significance.
# Plot expression changes
cao$plotExpressionShiftMagnitudes()

# Here, y-axis shows magnitude of changes, cells show both expression and composition changes.
