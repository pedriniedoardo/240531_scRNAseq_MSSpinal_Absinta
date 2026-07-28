# AIM ---------------------------------------------------------------------
# Aletta reviewed the 27_*_subcluster_HarmonySample.R / _plot.R outputs and flagged, per subpopulation, (1) contaminant clusters to remove and (2) clusters that need a refined annotation label.
# The calls is in ../../data/LUT_subcluster_filter_annotate.csv (subcluster_id, col_meta, cluster_id, action, new_label).
# This script does NOT touch the subcluster objects themselves: it only harvests barcodes + new metadata so they can be merged back onto the full integrated object (001_sobj_integrated_cleanup_manualAnnotation_updateDemyelination.rds) in 29_apply_filter_annotate_fullObject.R, where the object is subset and only the UMAP is recomputed (removal is minimal, so PCA/Harmony are not rerun).

# LIBRARIES ---------------------------------------------------------------
library(Seurat)
library(tidyverse)

# PARAMETERS ------------------------------------------------------------------
options(stringsAsFactors = FALSE)
options(Seurat.object.assay.version = "v5")

# path per subcluster object, keyed on the same subcluster_id used in the LUT
rds_lut <- c(
  ASTRO = "../../out/object/analysis_R44/27_ASTRO_subcluster_HarmonySample.rds",
  IMMUNE = "../../out/object/analysis_R44/27_IMMUNE_subcluster_HarmonySample.rds",
  LYM = "../../out/object/analysis_R44/27_LYM_subcluster_HarmonySample.rds",
  NEU = "../../out/object/analysis_R44/27_NEU_subcluster_HarmonySample.rds",
  OLIGO = "../../out/object/analysis_R44/27_OLIGO_subcluster_HarmonySample.rds",
  OPC = "../../out/object/analysis_R44/27_OPC_subcluster_HarmonySample.rds",
  STROMAL = "../../out/object/analysis_R44/27_STROMAL_subcluster_HarmonySample.rds",
  VAS = "../../out/object/analysis_R44/27_VAS_subcluster_HarmonySample.rds")

# read in the colleague's calls -------------------------------------------------
# load the annotation summary produced by aletta
# lut <- read_csv("../../data/LUT_subcluster_filter_annotate.csv") %>%
#   mutate(cluster_id = as.character(cluster_id))
lut <- read_csv("../../data/260813_spinal_subcluster_annotation_aletta.csv") %>%
  mutate(cluster_id = as.character(cluster_id)) %>%
  mutate(new_label = paste0(subcluster_id,"|",anno_aletta_short))

# extraction loop -------------------------------------------------------------
# pop_name <- "STROMAL"
list_out <- lapply(names(rds_lut), function(pop_name) {
  message("\n========== ", pop_name, " ==========")

  lut_pop <- lut %>% filter(subcluster_id == pop_name)

  # no Aletta review yet for this subcluster: keep every barcode, single flat label
  if (nrow(lut_pop) == 0) {
    message("  No Aletta annotation for this subcluster - keeping all barcodes, single label: ", pop_name)
    sobj <- readRDS(rds_lut[[pop_name]])

    df_remove <- data.frame(subcluster_id = character(0), barcode = character(0))
    df_annotate <- data.frame(barcode = colnames(sobj), subcluster_id = pop_name, cell_id_subcluster = pop_name)

    rm(sobj)
    gc()

    return(list(remove = df_remove, annotate = df_annotate))
  }

  col_meta <- unique(lut_pop$col_meta)

  output_plot <- paste0("../../out/plot/analysis_R44/29_", pop_name, "_subcluster_filterAnnotate.pdf")

  sobj <- readRDS(rds_lut[[pop_name]])
  Idents(sobj) <- col_meta

  # clusters to remove -------------------------------------------------------
  bad_clusters <- lut_pop %>% filter(action == "remove") %>% pull(cluster_id)
  if (length(bad_clusters) == 0) {
    message("  No clusters flagged for removal")
    barcodes_rm <- character(0)
    df_remove <- data.frame(subcluster_id = character(0), barcode = barcodes_rm)
  } else {
    barcodes_rm <- WhichCells(sobj, idents = bad_clusters)
    message("  Cells flagged for removal: ", length(barcodes_rm), " / ", ncol(sobj))
    df_remove <- data.frame(subcluster_id = pop_name, barcode = barcodes_rm)
  }

  # clusters to (re)annotate --------------------------------------------------
  lut_annotate <- lut_pop %>% filter(action == "none")
  if (nrow(lut_annotate) == 0) {
    message("  No clusters flagged for re-annotation")
    df_annotate <- tibble(barcode = character(0), subcluster_id = character(0), cell_id_subcluster = character(0))
  } else {
    df_annotate <- pmap_dfr(lut_annotate, function(cluster_id, new_label, ...) {
      data.frame(
        barcode = WhichCells(sobj, idents = cluster_id),
        subcluster_id = pop_name,
        cell_id_subcluster = new_label
      )
    })
    message("  Cells re-annotated: ", nrow(df_annotate))
  }
  
  # shortlit the table
  test_meta <- lut_annotate %>% dplyr::select(cluster_id,anno_aletta_short,action,note,concern,new_label)
  
  # add the new annotation
  meta_fix <- sobj@meta.data %>%
    data.frame() %>%
    rownames_to_column("barcode_id") %>%
    left_join(test_meta,by = setNames("cluster_id", col_meta[1])) %>%
    column_to_rownames("barcode_id")
  
  # update the meta
  sobj@meta.data <- meta_fix
  
  # diagnostic plot -----------------------------------------------------------
  sobj$flag_remove <- colnames(sobj) %in% barcodes_rm
  p01 <- DimPlot(sobj, group.by = "new_label", label = TRUE) + ggtitle("aletta new annotation")
  p02 <- DimPlot(sobj, group.by = "flag_remove", cols = c("TRUE" = "red", "FALSE" = "grey80")) + ggtitle("flagged for removal")
  ggsave(filename = output_plot, plot = p01 | p02, height = 6, width = 12)

  rm(sobj)
  gc()

  tot <- list(remove = df_remove, annotate = df_annotate)
  return(tot)
})

# combined tables --------------------------------------------------------------
# save the full list of barcode to remove and to reannotate
saveRDS(list_out,"../../out/object/analysis_R44/29_barcodes_annotation.rds")
