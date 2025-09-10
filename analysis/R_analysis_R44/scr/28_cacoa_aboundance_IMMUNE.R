# AIM ---------------------------------------------------------------------
# run a statistical test on the proportion differences across clusters, between conditions
# in this case I am running the analysis on the IMMUNE subset

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(cacoa)
# library(conos)
# library(sccore)
library(coda.base)
library(psych)
library(quadprog)

# read in the dataset -----------------------------------------------------
# load the seurat object
# test <- readRDS("../../out/object/analysis_R44/26_sobj_integrated_cleanup_manualAnnotation.rds")
# DimPlot(test, group.by = "cell_id")
# scobj_subset <- subset(test,subset = cell_id %in% c("IMMUNE"))

obj_test <- readRDS("../../out/object/analysis_R44/27_IMMUNE_subcluster_HarmonySample.rds")
DimPlot(obj_test, group.by = "RNA_snn_res.0.4")

# sum(!colnames(obj_test) %in% colnames(scobj_subset))

# wraingling --------------------------------------------------------------
# pull the metadata from the object
meta_test <- obj_test@meta.data

# run cacao ---------------------------------------------------------------
# sample.groups: vector with condition labels per sample named with sample ids
sample.groups_test <- meta_test$diagnosis_short
names(sample.groups_test) <- meta_test$sample_id

# cell.groups: cell type annotation vector named by cell ids
cell.groups_test <- meta_test$RNA_snn_res.0.4
names(cell.groups_test) <- rownames(meta_test)

# sample.per.cell: vector with sample labels per cell named with cell ids
sample.per.cell_test <- meta_test$sample_id
names(sample.per.cell_test) <- rownames(meta_test)

# ref.level: id of the condition, corresponding to the reference (i.e. control)
table(sample.groups_test)
ref.level <- "CTRL"

# target.level: id of the condition, corresponding to the target (i.e. case)
target.level <- "MS"

# creat the cacao object
cao_test <- Cacoa$new(data.object = obj_test,
                      sample.groups=sample.groups_test,
                      cell.groups=cell.groups_test,
                      sample.per.cell=sample.per.cell_test,
                      ref.level=ref.level,
                      target.level=target.level)

# run cacao Cluster-based changes
# The fastest way to visualize changes in the dataset is to show, what cell types
# are changing as a function of disease.
# Cacoa_estimateCellLoadings()
# estimateCellLoadings()
cao_test$estimateCellLoadings()

# Plot compositional changes
cao_test$plotCellLoadings(show.pvals=FALSE)
