# AIM ---------------------------------------------------------------------
# martina asked to furter celanup the object to remove cluster 19 from res 0.4

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)

# input -------------------------------------------------------------------
# read in the object
sobj_total_h <- readRDS("../../out/object/manualClean/20_sobj_integrated_cleanup.rds")
DimPlot(sobj_total_h,group.by = "RNA_snn_res.0.4",label = T)

# cleanup -----------------------------------------------------------------
# remove cluster 19 of res 0.4
# flag the barcodes
sobj_total_h$br_keep <- !sobj_total_h@meta.data$RNA_snn_res.0.4 %in% c(19)
sobj_total_h_clean <- subset(sobj_total_h,subset = br_keep == T)

# check the cleanup version
DimPlot(sobj_total_h_clean,group.by = "RNA_snn_res.0.4",label = T)

# save the metadata of the cleanup
saveRDS(sobj_total_h_clean@meta.data,file = "../../out/object/subclusters/23_meta_sobj_clean.rds")
