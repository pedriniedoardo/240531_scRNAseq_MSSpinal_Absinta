# AIM ---------------------------------------------------------------------
# apply the fourth-round cleanup/annotation calls (harvested in 02_filter_annotate_subclusters.R) to the full integrated reference object.
# The fourth-round LUT is a complete, self-contained consensus recomputed for all 8 populations directly from the raw 27_*_subcluster_HarmonySample.rds objects (not incremental on top of round 3, it re-flags the same round-3 clusters for removal too), so we start from the pristine pre-cleanup object rather than round 3's already-cleaned output - this keeps the "before cleanup" baseline and the "Cells removed" count honest across the full consensus.
# Adds cell_id_subcluster (short label), cell_id_subcluster2 (short label + cluster_id), and cell_id_subcluster3 (consensus long-form label) metadata columns, defaulting to the coarse cell_id where no subcluster call is available.
# Only a minimal number of cells are removed, so we do not rerun PCA/Harmony: the existing "harmony" reduction is reused and only the UMAP embedding is recomputed on the filtered object.

# LIBRARIES ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(patchwork)

# PARAMETERS ------------------------------------------------------------------
options(stringsAsFactors = FALSE)
options(Seurat.object.assay.version = "v5")

# read in the data --------------------------------------------------------
sobj <- readRDS("../../out/object/analysis_R44/001_sobj_integrated_cleanup_manualAnnotation_updateDemyelination.rds")

# list_out is a list with one element per subcluster, each a list(remove = <df>, annotate = <df>)
list_out <- readRDS("../../out/object/analysis_R45_pixi/02_barcodes_annotation_v2.rds")

# extract the tables
barcodes_rm <- lapply(list_out,function(x){
  x[["remove"]]
}) %>%
  bind_rows()

meta_annotate <- lapply(list_out,function(x){
  x[["annotate"]]
}) %>%
  bind_rows()

DimPlot(sobj, group.by = "cell_id", label = TRUE, raster = TRUE) + ggtitle("before cleanup")

# wrangling -----------------------------------------------------------------
# default the refined annotation to the coarse cell_id, then override with the subcluster calls
meta_new <- sobj@meta.data %>%
  select(barcode, cell_id) %>%
  left_join(meta_annotate %>% select(barcode, cell_id_subcluster,cell_id_subcluster2,cell_id_subcluster3), by = "barcode") %>%
  mutate(cell_id_subcluster = case_when(
    !is.na(cell_id_subcluster) ~ cell_id_subcluster,
    TRUE ~ cell_id
  )) %>%
  mutate(cell_id_subcluster2 = case_when(
    !is.na(cell_id_subcluster2) ~ cell_id_subcluster2,
    TRUE ~ cell_id
  )) %>%
  mutate(cell_id_subcluster3 = case_when(
    !is.na(cell_id_subcluster3) ~ cell_id_subcluster3,
    TRUE ~ cell_id
  )) %>%
  select(barcode, cell_id_subcluster,cell_id_subcluster2,cell_id_subcluster3)

sobj <- AddMetaData(sobj, meta_new)

# confirm the annotation update
table(sobj$cell_id_subcluster,sobj$cell_id, useNA = "ifany")
table(sobj$cell_id_subcluster2,sobj$cell_id, useNA = "ifany")
table(sobj$cell_id_subcluster3,sobj$cell_id, useNA = "ifany")

# remove the flagged contaminant barcodes ------------------------------------
message("Cells before cleanup: ", ncol(sobj))
sobj_clean <- subset(sobj, cells = setdiff(colnames(sobj), barcodes_rm$barcode))
message("Cells removed: ", ncol(sobj) - ncol(sobj_clean))
message("Cells after cleanup: ", ncol(sobj_clean))

# rerun only the UMAP, reuse the existing harmony reduction -----------------
sobj_clean <- sobj_clean %>%
  RunUMAP(reduction = "harmony", dims = 1:30, return.model = TRUE)

table(sobj_clean$cell_id_subcluster,sobj_clean$cell_id, useNA = "ifany")
table(sobj_clean$cell_id_subcluster2,sobj_clean$cell_id, useNA = "ifany")
table(sobj_clean$cell_id_subcluster3,sobj_clean$cell_id, useNA = "ifany")

# save the new reference ------------------------------------------------------
saveRDS(sobj_clean, "../../out/object/analysis_R45_pixi/03_sobj_integrated_cleanup_manualAnnotation_subclusterFiltered_v2.rds")

# diagnostic plots -------------------------------------------------------------
p01 <- DimPlot(sobj, group.by = "cell_id", label = TRUE, raster = TRUE) + ggtitle("before cleanup")
p02 <- DimPlot(sobj_clean, group.by = "cell_id", label = TRUE, raster = TRUE) + ggtitle("after cleanup - rerun UMAP")
ggsave(filename = "../../out/plot/analysis_R45_pixi/03_UMAP_fullObject_beforeAfter_cleanup_v2.pdf", plot = p01 | p02, height = 6, width = 14)
