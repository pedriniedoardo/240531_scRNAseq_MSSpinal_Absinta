# AIM ---------------------------------------------------------------------
# explore etherogeneity of the plasmablast population in the LYM subcluster

# renv integration --------------------------------------------------------
# to load the packages
# source(".Rprofile")

# libraries -----------------------------------------------------------------
# renv::install("phipsonlab/speckle")
# renv::install("statmod")
library(speckle)
library(limma)
library(statmod)
library(cowplot)
library(ggrepel)
library(finalfit)
library(Seurat)
library(tidyverse)
library(presto)
library(pals)
library(patchwork)
library(harmony)
library(scCustomize)
library(Nebulosa)
library(viridis)

# define the inputs ----------------------------------------------------------
# subset dataset
sobj_input_id <- "../../out/object/analysis_R44/27_LYM_subcluster_HarmonySample.rds"
message("input single cell object: ", sobj_input_id)

sobjFull_input_id <- "../../out/object/analysis_R44/29_sobj_integrated_cleanup_manualAnnotation_subclusterFiltered.rds"
message("input single cell full object: ", sobjFull_input_id)

# read in the object ----------------------------------------------------------
sobj <- readRDS(sobj_input_id)
sobj_full <- readRDS(sobjFull_input_id)

DefaultAssay(sobj) <- "RNA"
DefaultAssay(sobj_full) <- "RNA"

# confirm the identity of the dataset
DimPlot(sobj, label = T, raster = F,repel = T, group.by = "cell_id")
DimPlot(sobj_full, label = T, raster = T,repel = T, group.by = "cell_id")

# wrangling ---------------------------------------------------------------
# add the annotation to the subcluster
sobj <- AddMetaData(sobj,metadata = sobj_full@meta.data[,c("cell_id_subcluster",
                                                           "cell_id_subcluster2")])
DimPlot(sobj,group.by = "cell_id_subcluster",label = T, repel = T) + NoLegend()

# subset only the plasmablast and rerun the processing
sobj_subcluster <- subset(sobj,subset = cell_id_subcluster %in% c("LYM|Plasma"))
DimPlot(sobj_subcluster,group.by = "cell_id_subcluster",label = T, repel = T) + NoLegend()

# rescale the data for regressing out the sources of variation do not scale all the genes. if needed I can scale them before the heatmap call. for speeding up the computation I will keep 
# drop parent clustering columns — keep only sample-level metadata
cols_drop <- grep("RNA_snn_res|^seurat_clusters$", colnames(sobj_subcluster@meta.data), value = TRUE)
meta_keep <- sobj_subcluster@meta.data[, setdiff(colnames(sobj_subcluster@meta.data), cols_drop), drop = FALSE]

ct <- "Plasma"
# build fresh object from raw counts only
sobj_subcluster_reprocess <- CreateSeuratObject(
  counts = sobj_subcluster@assays$RNA$counts,
  project = paste0(ct, "_subcluster"),
  meta.data = meta_keep,
  min.cells = 0,
  min.features = 0
) %>%
  NormalizeData(verbose = FALSE)

gc()

# sanity check: orig.ident must retain real per-sample values here, not the single CreateSeuratObject() `project` label, or RunHarmony() below will have nothing to correct for
table(sobj_subcluster_reprocess$sample_id)
table(sobj_subcluster_reprocess$orig.ident)

# subcluster QC / integration parameters -----------------------------------
# Tune the parameters.
# This is a small, rare population (~2-3k cells across ~23 samples), so the generic defaults below (tuned for the much larger parent population) drive Harmony to over-fragment the data into noise-sized clusters.
min_cells_per_sample <- 20  # diagnostic-only threshold (no package default, no filtering applied) — flags samples too sparse to reliably inform Harmony's correction; all cells/donors are still kept, imbalance is instead handled via harmony_nclust/harmony_theta below
n_features <- 1500          # FindVariableFeatures(nfeatures=); Seurat default 2000, script previously used 3000 (more features than cells at this N invites PCA to fit noise)
n_pcs <- 15                 # RunPCA(npcs=) and dims= used downstream; Seurat default npcs = 50, script previously used 30 (PCs beyond ~15 at this N are mostly technical/sample variance)
harmony_theta <- 1          # RunHarmony(theta=); harmony default 2 (higher = more aggressive forcing of every batch into every cluster)
harmony_nclust <- 10        # RunHarmony(nclust=); harmony default min(round(N/30), 100) (~80 here) — scales Harmony's internal soft-clustering to the few real states expected in a rare, fairly homogeneous population
cluster_res <- seq(0.1, 0.5, by = 0.1)  # Seurat FindClusters default resolution = 0.8; script previously swept up to 1.0, which reliably over-clusters a dataset this size

# diagnose per-sample imbalance before Harmony -------------------------------
sample_counts <- sort(table(sobj_subcluster_reprocess$sample_id))
print(sample_counts)
n_below <- sum(sample_counts < min_cells_per_sample)

# NOTE: keeping all cells / all donors here on purpose (no sample_id filtering) - the imbalance flagged above is instead handled downstream via harmony_nclust/harmony_theta

# pre-Harmony pipeline ----------------------------------------------------
set.seed(42)
sobj_subcluster_reprocess <- sobj_subcluster_reprocess %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = n_features, verbose = FALSE) %>%
  ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA", "S.Score", "G2M.Score"),
            verbose = FALSE) %>%
  RunPCA(npcs = n_pcs, verbose = FALSE) %>%
  # optional pre-Harmony EDA: RunHarmony() below only consumes the "pca" reduction, so this UMAP/neighbors/clustering block is not required for integration — it's kept only as an optional "before integration" preview (e.g. to eyeball batch effects) and can be commented out if not needed
  RunUMAP(reduction = "pca", dims = 1:n_pcs, return.model = TRUE, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:n_pcs, verbose = FALSE) %>%
  FindClusters(resolution = cluster_res, verbose = FALSE)

sobj_subcluster_reprocess

DimPlot(sobj_subcluster_reprocess,group.by = "RNA_snn_res.0.1") + DimPlot(sobj_subcluster_reprocess,group.by = "sample_id")
ggsave("../../out/plot/analysis_R45_pixi/01_sample_UMAP.pdf",width = 12,height = 5)

# Harmony integration -----------------------------------------------------
set.seed(42)
sobj_subcluster_reprocess_h <- sobj_subcluster_reprocess %>%
  RunHarmony("sample_id", theta = harmony_theta, nclust = harmony_nclust, plot_convergence = FALSE) %>%
  RunUMAP(reduction = "harmony", dims = 1:n_pcs, return.model = TRUE, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = 1:n_pcs, verbose = FALSE) %>%
  FindClusters(resolution = cluster_res, verbose = FALSE)

sobj_subcluster_reprocess_h

DimPlot(sobj_subcluster_reprocess_h,group.by = "RNA_snn_res.0.1") + DimPlot(sobj_subcluster_reprocess_h,group.by = "sample_id")
ggsave("../../out/plot/analysis_R45_pixi/01_sample_UMAP_harmony.pdf",width = 12,height = 5)

# verification: confirm clusters are not just single-sample artifacts -------
# a cluster that's dominated (e.g. >70-80%) by one sample_id is still batch-driven noise, not biology
table(sobj_subcluster_reprocess_h$RNA_snn_res.0.1, sobj_subcluster_reprocess_h$sample_id)
round(prop.table(table(sobj_subcluster_reprocess_h$RNA_snn_res.0.1, sobj_subcluster_reprocess_h$sample_id), margin = 1), 2)

# set primary resolution as default active identity
Idents(sobj_subcluster_reprocess_h) <- "RNA_snn_res.0.1"

# DGE: cluster 1 vs cluster 0 (resolution 0.1, Harmony-corrected) -----------
markers_1vs0 <- FindMarkers(
  sobj_subcluster_reprocess_h,
  ident.1 = "1",
  ident.2 = "0",
  logfc.threshold = 0,
  min.pct = 0.1,
  verbose = FALSE
) %>%
  rownames_to_column("gene") %>%
  arrange(p_val_adj)

write_tsv(markers_1vs0, "../../out/table/analysis_R45_pixi/01_DGE_cluster1_vs_cluster0_res0.1_harmony.tsv")

# volcano plot ----------------------------------------------------------
volcano_logfc_cut <- 0.5   # |avg_log2FC| threshold for calling a gene up/down
volcano_padj_cut <- 0.05   # p_val_adj threshold for significance
n_label <- 15              # top genes per direction (by p_val_adj) to label on the volcano

df_volcano <- markers_1vs0 %>%
  mutate(
    direction = case_when(
      p_val_adj < volcano_padj_cut & avg_log2FC >  volcano_logfc_cut ~ "up in cluster 1",
      p_val_adj < volcano_padj_cut & avg_log2FC < -volcano_logfc_cut ~ "up in cluster 0",
      TRUE ~ "ns"
    )
  )

genes_label <- df_volcano %>%
  filter(direction != "ns") %>%
  group_by(direction) %>%
  slice_min(p_val_adj, n = n_label) %>%
  ungroup()

df_volcano %>%
  ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj), col = direction)) +
  geom_point(alpha = 0.5) +
  geom_text_repel(data = genes_label, aes(label = gene), size = 3, max.overlaps = Inf) +
  geom_vline(xintercept = c(-volcano_logfc_cut, volcano_logfc_cut), linetype = "dashed") +
  geom_hline(yintercept = -log10(volcano_padj_cut), linetype = "dashed") +
  scale_color_manual(values = c("up in cluster 1" = "firebrick", "up in cluster 0" = "steelblue", "ns" = "grey70")) +
  theme_bw() +
  labs(title = "Cluster 1 vs Cluster 0 (RNA_snn_res.0.1, Harmony)", col = NULL)
ggsave("../../out/plot/analysis_R45_pixi/01_volcano_cluster1_vs_cluster0_res0.1_harmony.pdf", width = 7, height = 6)

# is there in any cluster our favourite set of markers?
df_volcano %>%
  filter(gene %in% c("IGKC","IGHG1"))

# FeaturePlot(sobj_subcluster_reprocess_h,
#             features = c("IGKC","IGHG1"),
#             reduction = "umap") &
#   scale_color_viridis_c(option = "inferno")

pal <- viridis(n = 10,option = "B",direction = 1)
FeaturePlot_scCustom(seurat_object = sobj_subcluster_reprocess_h,
                     features = c("IGKC","IGHG1"),
                     order = T,
                     na_cutoff = 0,
                     colors_use = pal)

Plot_Density_Custom(seurat_object = sobj_subcluster_reprocess_h,
                    features = c("IGKC","IGHG1"))

# Stacked_VlnPlot(seurat_object = sobj_subcluster_reprocess_h,
#                 features = c("IGKC","IGHG1"),
#                 x_lab_rotate = TRUE)

VlnPlot(sobj_subcluster_reprocess_h,
        features = c("IGKC","IGHG1"))


# quick pathway analysis (EnrichR) on the significant DGE genes -------------
library(enrichR)

dbs_db <- c("KEGG_2021_Human", "MSigDB_Hallmark_2020", "Reactome_Pathways_2024",
            # cell-type marker libraries: Azimuth matches the HumanPBMCRef reference already
            # used elsewhere in this project; CellMarker_2024 is more granular/curated;
            # PanglaoDB_Augmented_2021 is broader tissue coverage but coarser
            "Azimuth_Cell_Types_2021", "CellMarker_2024", "PanglaoDB_Augmented_2021")

run_enrichr_safe <- function(genes, dbs) {
  if (length(genes) < 5) {
    message("fewer than 5 genes, skipping EnrichR for this direction")
    return(tibble())
  }
  out_enrich <- enrichr(genes, dbs)
  keep <- vapply(out_enrich, nrow, integer(1)) > 0
  out_enrich[keep] %>% bind_rows(.id = "annotation")
}

genes_up_c1 <- df_volcano %>% filter(direction == "up in cluster 1") %>% pull(gene)
genes_up_c0 <- df_volcano %>% filter(direction == "up in cluster 0") %>% pull(gene)

enrich_up_c1 <- run_enrichr_safe(genes_up_c1, dbs_db)
enrich_up_c0 <- run_enrichr_safe(genes_up_c0, dbs_db)

write_tsv(enrich_up_c1, "../../out/table/analysis_R45_pixi/01_enrichR_cluster1up_res0.1_harmony.tsv")
write_tsv(enrich_up_c0, "../../out/table/analysis_R45_pixi/01_enrichR_cluster0up_res0.1_harmony.tsv")

plot_enrich <- function(df, title) {
  df %>%
    group_by(annotation) %>%
    arrange(P.value) %>%
    dplyr::slice(1:10) %>%
    mutate(Term = str_sub(Term, 1, 40)) %>%
    mutate(Term = fct_reorder(Term, Combined.Score)) %>%
    ggplot(aes(y = Term, x = Combined.Score, size = Odds.Ratio, col = Adjusted.P.value)) +
    geom_point() +
    facet_wrap(~annotation, scales = "free", ncol = 1) +
    theme_bw() +
    theme(strip.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA)) +
    ggtitle(title)
}

plot_enrich(enrich_up_c1, "Up in cluster 1") + plot_enrich(enrich_up_c0, "Up in cluster 0")
ggsave("../../out/plot/analysis_R45_pixi/01_enrichR_cluster1_vs_cluster0_res0.1_harmony.pdf", width = 16, height = 10)

# EDA ---------------------------------------------------------------------
# pull more makers from the terms of interest from Azimuth
enrich_up_c0 %>%
  filter(str_detect(annotation,pattern = "Azimuth")) %>%
  filter(str_detect(Term,pattern = "Plasma"))

panel_01 <- c("CD79A","DERL3","MZB1","JCHAIN","TXNDC5")

# check the stats
df_volcano %>%
  filter(gene %in% panel_01)

# Create Plots
pal <- viridis(n = 10,option = "B")
FeaturePlot_scCustom(seurat_object = sobj_subcluster_reprocess_h,
                     features = panel_01,
                     order = T,
                     na_cutoff = 0,
                     colors_use = pal,num_columns = 3)

Plot_Density_Custom(seurat_object = sobj_subcluster_reprocess_h, features = panel_01)

Stacked_VlnPlot(seurat_object = sobj_subcluster_reprocess_h,
                features = panel_01,
                x_lab_rotate = TRUE)


# pull more makers from the terms of interest from MsigDB
enrich_up_c0 %>%
  filter(str_detect(annotation,pattern = "MSigDB")) %>%
  filter(str_detect(Term,pattern = "Interferon"))

panel_02 <- c("BST2", "CD74", "MTHFD2", "BANK1", "IRF4", "TXNIP", "CD38", "XAF1", "NCOA7", "PARP9")

# check the stats
df_volcano %>%
  filter(gene %in% panel_02)

# Create Plots
# pal <- viridis(n = 10,option = "B")
FeaturePlot_scCustom(seurat_object = sobj_subcluster_reprocess_h,
                     features = panel_02,
                     order = T,
                     na_cutoff = 0,
                     colors_use = pal,num_columns = 4)

Plot_Density_Custom(seurat_object = sobj_subcluster_reprocess_h, features = panel_02)

Stacked_VlnPlot(seurat_object = sobj_subcluster_reprocess_h,
                features = panel_02,
                x_lab_rotate = TRUE)


