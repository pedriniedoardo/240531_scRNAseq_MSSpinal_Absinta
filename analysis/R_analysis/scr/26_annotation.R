# AIM ---------------------------------------------------------------------
# define the final onnotation of the clusters following automatic and manual annotation informations

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)

# read in the data --------------------------------------------------------
scobj <- readRDS("../../out/object/manualClean/24_sobj_integrated_cleanup.rds")
Idents(scobj) <- "cell_id_old"
Idents(scobj) <- "RNA_snn_res.0.1"

# show the UMAP with the former anntoation
DimPlot(scobj,label = T)

# Add the cell_id based on martina's annotation
scobj$cell_id <- scobj@meta.data %>%
  mutate(cell_id = case_when(RNA_snn_res.0.1 %in% c(0)~"OLIGO",
                             RNA_snn_res.0.1 %in% c(4)~"OPC",
                             RNA_snn_res.0.1 %in% c(2,12)~"ASTRO",
                             RNA_snn_res.0.1 %in% c(1)~"IMMUNE",
                             RNA_snn_res.0.1 %in% c(9)~"B CELLS",
                             RNA_snn_res.0.1 %in% c(5)~"T CELLS",
                             RNA_snn_res.0.1 %in% c(8)~"VAS",
                             RNA_snn_res.0.1 %in% c(11)~"SCHWANN",
                             RNA_snn_res.0.1 %in% c(10)~"EPENDYMA",
                             RNA_snn_res.0.1 %in% c(3,6)~"STROMAL",
                             RNA_snn_res.0.1 %in% c(7)~"NEU")) %>%
  pull(cell_id)

# confirm the annotation
DimPlot(scobj,group.by = "cell_id",raster=T,label=T)

# save the objects --------------------------------------------------------
# save the full metadata for the annotated object
# meta <- scobj@meta.data %>%
#   rownames_to_column()
# write_tsv(meta,"../../out/table/manualClean/meta.data_harmonySkipIntegration_AllSoupX_00500_07000_05_AnnotationSCType_manualAnnotation.tsv")

# save the object
saveRDS(scobj,"../../out/object/manualClean/26_sobj_integrated_cleanup_manualAnnotation.rds")
