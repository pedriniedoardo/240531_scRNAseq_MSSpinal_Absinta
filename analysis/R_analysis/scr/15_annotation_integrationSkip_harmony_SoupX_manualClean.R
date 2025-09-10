# AIM ---------------------------------------------------------------------
# define the final onnotation of the clusters following automatic and manual annotation informations

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)

# read in the data --------------------------------------------------------
scobj <- readRDS("../../out/object/manualClean/data.combined_harmonySkipIntegration_AllSoupX_00500_07000_05_AnnotationSCType.rds")

# Add the cell_id based on martina's annotation
scobj$cell_id <- scobj@meta.data %>%
  mutate(cell_id = case_when(RNA_snn_res.0.2 %in% c(0,3,11)~"OLIGO",
                             RNA_snn_res.0.2 %in% c(4)~"OPC",
                             RNA_snn_res.0.2 %in% c(1)~"ASTRO",
                             RNA_snn_res.0.2 %in% c(2)~"IMMUNE",
                             RNA_snn_res.0.2 %in% c(12)~"B CELLS",
                             RNA_snn_res.0.2 %in% c(5)~"T CELLS",
                             RNA_snn_res.0.2 %in% c(10)~"VAS",
                             RNA_snn_res.0.2 %in% c(15)~"SCHWANN",
                             RNA_snn_res.0.2 %in% c(13)~"EPENDYMA",
                             RNA_snn_res.0.2 %in% c(6,7,8)~"STROMAL",
                             RNA_snn_res.0.2 %in% c(9,14)~"NEU")) %>%
  pull(cell_id)

# confirm the annotation
DimPlot(scobj,group.by = "cell_id",raster=T,label=T)

# save the objects --------------------------------------------------------
# save the full metadata for the annotated object
meta <- scobj@meta.data %>%
  rownames_to_column()
write_tsv(meta,"../../out/table/manualClean/meta.data_harmonySkipIntegration_AllSoupX_00500_07000_05_AnnotationSCType_manualAnnotation.tsv")

# save the object
saveRDS(scobj,"../../out/object/manualClean/data.combined_harmonySkipIntegration_AllSoupX_00500_07000_05_AnnotationSCType_manualAnnotation.rds")
