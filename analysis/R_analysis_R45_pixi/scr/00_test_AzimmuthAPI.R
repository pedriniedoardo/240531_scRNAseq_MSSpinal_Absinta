# AIM ---------------------------------------------------------------------
# test azimuthAPI for the evaluation of the automatic annotation of the dataset

# REF ---------------------------------------------------------------------
# Pan-human Azimuth (AzimuthAPI) — intro vignette as a runnable script
# Source: https://satijalab.org/pan_human_azimuth/r/articles/intro_azimuthapi_vignette

# This vignette demonstrates how to install and use AzimuthAPI, an interface for cell type annotation using the Pan-human Azimuth neural network.
# The API takes a Seurat object and returns hierarchical cell type predictions, confidence scores, and low-dimensional embeddings useful for visualization.

# This is the test to run the tool locally. before preceeding make sure to:
# - create the conda env (I decided to create a local project folder, but it can be also a general conda env). The yaml is located in: "../env//env_panhumanpy.yml"
# - download the model weights ('https://zenodo.org/records/20401417/files/panhumanpy_inference_model_v1.keras?download=1') and link them into the expected cache folder (~/.cache/panhumanpy/v1/inference_model/inference_model.keras)
# - activate the python env inside the R script (Sys.setenv(RETICULATE_PYTHON = "./.conda-pan-human-azimuth/bin/python"))

# libraries ---------------------------------------------------------------
library(AzimuthAPI)
library(Seurat)
library(SeuratData)
library(reticulate)
library(tidyverse)
library(presto)
library(patchwork)
library(ComplexHeatmap)
library(viridis)
library(grid)

# custom functions --------------------------------------------------------
# define the jaccard score function
jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  return (intersection/union)
}

# build the dataset for the correlation plot crossing the two annotations per cell, and compute the jaccard score per pair (factored into a helper since we need this twice: azimuth clusters vs. labels, and PCA clusters vs. labels)
compute_jaccard_matrix <- function(meta, id_ref_col, id_query_col) {
  df_crossing <- crossing(id_ref = unique(meta[[id_ref_col]]),
                          id_query = unique(meta[[id_query_col]]))
  
  df_jaccard_score <- pmap(list(id_ref = df_crossing$id_ref,
                                id_query = df_crossing$id_query), function(id_ref, id_query){
                                  
                                  a <- meta %>%
                                    filter(.data[[id_ref_col]] == id_ref) %>% pull(cell_barcode)
                                  
                                  b <- meta %>%
                                    filter(.data[[id_query_col]] == id_query) %>% pull(cell_barcode)
                                  
                                  jaccard_score <- jaccard(a, b)
                                  
                                  df <- data.frame("id_ref" = id_ref,
                                                   "id_query" = id_query,
                                                   "jaccard_score" = jaccard_score)
                                  return(df)
                                }) %>%
    bind_rows()
  
  # shape it as a matrix
  df_jaccard_score %>%
    pivot_wider(names_from = id_ref, values_from = jaccard_score) %>%
    column_to_rownames("id_query")
}

# build a Jaccard-score Heatmap paired with the side-by-side DimPlots it compares (ref_col vs query_col, plus an optional third "majority vote" panel), so the same figure can be produced for multiple annotation/cluster comparisons.
# ref_col and query_col double as both the DimPlot group.by and the heatmap column/row titles.
# query_reduction defaults to ref_reduction, but pass both explicitly when the two annotations live natively in different reductions (e.g. an original annotation on "umap" vs an azimuth-derived one on "azimuth_umap").
# majority_col, when given, is a metadata column relabelling every cell of a query_col cluster with that cluster's dominant ref_col annotation -- plotted as a third UMAP so mixed vs. clean clusters are easy to spot against the raw per-cell ref/query views.
plot_jaccard_comparison <- function(mat_jaccard,
                                     sobj,
                                     ref_col,
                                     query_col,
                                     ref_reduction,
                                     query_reduction = ref_reduction,
                                     majority_col = NULL,
                                     majority_reduction = ref_reduction,
                                     label_size = 2.5) {
  ht <- Heatmap(mat_jaccard,
                name = "Jaccard score",
                col = viridis::viridis(option = "turbo", n = 20),
                row_names_side = "right",
                row_names_gp = gpar(fontsize = 8),
                column_names_side = "bottom",
                column_names_gp = gpar(fontsize = 8),
                row_dend_reorder = FALSE,
                column_dend_reorder = FALSE,
                row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                column_title = ref_col,
                row_title = query_col,
                show_column_names = TRUE,
                show_row_names = TRUE)

  p_ref <- DimPlot(sobj, group.by = ref_col, label.size = label_size, label = TRUE,
                    reduction = ref_reduction, repel = TRUE) + NoLegend()

  p_query <- DimPlot(sobj, group.by = query_col, label.size = label_size, label = TRUE,
                      reduction = query_reduction, repel = TRUE) + NoLegend()

  plot_list <- list(p_ref, p_query)
  width_list <- c(1, 1)

  if (!is.null(majority_col)) {
    p_majority <- DimPlot(sobj, group.by = majority_col, label.size = label_size, label = TRUE,
                          reduction = majority_reduction, repel = TRUE) + NoLegend()
    plot_list <- c(plot_list, list(p_majority))
    width_list <- c(width_list, 1)
  }

  # grid.grabExpr(draw(x)) turns the ComplexHeatmap output into a grob patchwork can lay out alongside the ggplots
  wrap_plots(plot_list, nrow = 1) + grid.grabExpr(draw(ht)) +
    plot_layout(widths = c(width_list, 3))
}

# parameters --------------------------------------------------------------
# threshold of confidence for the annotation
thr_azimuth <- 0.5

# set seurat compatible with seurat5 workflow
# options(Seurat.object.assay.version = "v5")

# future parameters
options(future.globals.maxSize = 1000 * 1024^2 * 2)

# 2. Point reticulate at the local panhumanpy Python environment -------------
# Pin RETICULATE_PYTHON directly to the env's interpreter (must be set before reticulate initializes any Python session).
# This is more robust than use_condaenv(name) on this machine, since that relies on `conda env list` to discover environments, which does not reliably see micromamba envs.
# Adjust the path if your pan-human-azimuth env lives elsewhere (check with: micromamba env list).
Sys.setenv(RETICULATE_PYTHON = "./.conda-pan-human-azimuth/bin/python")

# Confirm reticulate actually loaded THIS interpreter and that it can see panhumanpy inside it. 
py_config()
# might expose token
# Sys.getenv()
Sys.getenv(c("RETICULATE_PYTHON", "CONDA_PREFIX", "CONDA_DEFAULT_ENV"))

cat("Active Python:", py_config()$python, "\n")
cat("panhumanpy version:", as.character(import("panhumanpy")$`__version__`), "\n")

# read in the data --------------------------------------------------------
# the input dataset is already normalized
sobj <- readRDS("../../out/object/analysis_R44/29_sobj_integrated_cleanup_manualAnnotation_subclusterFiltered.rds")

# explore the dataset
sobj@meta.data %>%
  head()

DimPlot(sobj,group.by = "cell_id",repel = T,label = T)

# ---- Annotate a dataset with Pan-human Azimuth ------------------------------
# Users have two options for running Pan-human Azimuth on their dataset: cloud-based or local.
# The cloud-based option is the easiest to use, and requires no additional setup.
# The local option requires setting up a Python environment with panhumanpy (Python package for the Pan-human Azimuth neural network) and its dependencies.
# Annotation results, regardless of the method used, are stored in the Seurat object as cell-level metadata, and the embeddings generated by the underlying neural network are stored in a new Seurat reduction.

# Option 2: run `ANNotate` (requires panhumanpy in a Python env reticulate can see)
# run the minimal outpur that should replicate the cloud/remote version of the analysis
t_init_minimal <- Sys.time()
sobj_ann <- ANNotate(sobj, output_mode = "minimal")
t_end_minimal <- Sys.time()
# time
t_end_minimal - t_init_minimal

saveRDS(sobj_ann, "../../out/object/analysis_R45_pixi/00_sobj_integrated_cleanup_manualAnnotation_subclusterFiltered_ann.rds")
# sobj_ann <- readRDS("../../out/object/analysis_R45_pixi/00_sobj_integrated_cleanup_manualAnnotation_subclusterFiltered_ann.rds")

# ---- Inspect the annotation output -----------------------------------------
# ANNotate() (process_obj = TRUE by default) adds: metadata columns final_level_labels,
# final_level_confidence, level_1_labels...level_N_labels (detailed mode), azimuth_label;
# and reductions "azimuth_embed" (raw NN embedding) and "umapANNazimuth_umap" (panhumanpy's
# own internal UMAP of that embedding, not used below - we build our own via RunUMAP instead).
sobj_ann@meta.data %>% head()

# check the distribution of the confidcence score
sobj_ann@meta.data %>%
  ggplot(aes(x=final_level_confidence)) + geom_histogram()

# We see that in addition to cell type annotations organized hierarchically, Pan-human Azimuth provides confidence scores for each annotation; lower confidence scores indicate uncertainty in hierarchy assignments. We can set a confidence threshold to retain only high-confidence annotations for downstream tasks.

# make the annotation
sobj_ann$azimuth_confidence <- !is.na(sobj_ann$final_level_confidence) & sobj_ann$final_level_confidence > thr_azimuth

# check what cell would be lost if performing the filtering
DimPlot(sobj_ann,group.by = "cell_id",split.by = "azimuth_confidence")
# plot the confidence score on the original dim reduction
FeaturePlot(sobj_ann,features = "final_level_confidence")

# sobj_ann_hc <- subset(sobj_ann,
#                       !is.na(final_level_confidence) & final_level_confidence > thr_azimuth)

# check the levels after filtering
sobj_ann@meta.data %>%
  group_by(final_level_labels) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

# saveRDS(sobj_azimuth_hc, "../../out/object/analysis_R45_pixi/00_sobj_azimuth_detailed_hc.rds")

# ---- Generate a UMAP from the azimuth_embed reduction and visualize --------
# The azimuth_embed reduction stores the Pan-human Azimuth embeddings. This is useful in downstream analyses. For example, we can visualize the cell type predictions returned by Pan-human Azimuth in embedding space:
sobj_ann <- RunUMAP(sobj_ann,
                    dims = 1:ncol(Embeddings(sobj_ann, "azimuth_embed")),
                    reduction = "azimuth_embed",
                    reduction.name = "azimuth_umap",
                    reduction.key = "UMAPazimuth_") %>%
  # NOTE: FindNeighbors()'s default graph.name is derived from DefaultAssay(), which is still "RNA" here -- without an explicit graph.name this would silently overwrite the "pca_nn"/"pca_snn" graphs from the PCA FindNeighbors() call above.
  FindNeighbors(dims = 1:ncol(Embeddings(sobj_ann, "azimuth_embed")),
                reduction = "azimuth_embed",
                graph.name = c("azimuth_nn", "azimuth_snn")) %>%
  FindClusters(graph.name = "azimuth_snn",
               resolution = seq(0.2, 1, by = 0.2))

# Save the final object
saveRDS(sobj_ann, "../../out/object/analysis_R45_pixi/00_sobj_ann_AzimuthAPI.rds")
# sobj_ann <- readRDS("../../out/object/analysis_R45_pixi/00_sobj_ann_AzimuthAPI.rds")

# plot the range of resolutions from the new embedding
# try to plot both azimuth annotation and the original one in the new embedding
p_azimuth_umap <- DimPlot(sobj_ann,
                          group.by = "final_level_labels",
                          label.size = 2.5,
                          label = TRUE,
                          reduction = "azimuth_umap",
                          repel = TRUE) + NoLegend()
p_orig_umap <- DimPlot(sobj_ann,
                       group.by = "final_level_labels",
                       label.size = 2.5,
                       label = TRUE,
                       reduction = "umap",
                       repel = TRUE) + NoLegend()

p_azimuth_umap | p_orig_umap

# ggsave(plot = p_azimuth_umap,
#        filename = "../../out/plot/analysis_R45_pixi/00_UMAP_azimuth_embed_finalLevelLabels.pdf",
#        height = 6, width = 7)

# How many levels are available for each annotation
sobj_ann@meta.data %>%
  group_by(azimuth_broad) %>%
  summarise(n = n())

sobj_ann@meta.data %>%
  group_by(azimuth_medium) %>%
  summarise(n = n())

sobj_ann@meta.data %>%
  group_by(azimuth_fine) %>%
  summarise(n = n())

sobj_ann@meta.data %>%
  group_by(final_level_labels) %>%
  summarise(n = n())

# try to simplify the annotation to show only robusts annotation
# use the helper prep labels to simplyfy the annotaiton if there are less than 200 cells label them as low-freq-labels
sobj_ann <- PrepLabel(sobj_ann,
                       label_id = "final_level_labels",
                       newid = "final_level_labels_200",
                       cutid = "low-freq-labels",
                       cutoff = 200)

# confirm the labels
sobj_ann@meta.data %>%
  group_by(final_level_labels_200) %>%
  summarise(n = n(),.groups = "drop") %>%
  arrange(desc(n)) %>%
  print(n = 50)

# ---- Jaccard similarity between cluster assignments and cell type labels ----
# use the Jaccard score to measure the cross-cluster similarity per cell (how similar are the clusters from the query compared to the annotation derived from the reference)

# note: sobj already has a metadata column literally named "barcode", so use a
# different name here for the rownames-to-column step to avoid a collision.
meta_hc <- sobj_ann@meta.data %>% rownames_to_column("cell_barcode")

# Heatmap A: azimuth_embed-based clusters vs. Pan-human Azimuth cell type labels

# compare the automatic annotation with simplified markers with the one suggested by Aletta
mat_jaccard_01 <- compute_jaccard_matrix(meta_hc, "cell_id_subcluster2","final_level_labels_200")

p_jaccard_azimuth_01 <- plot_jaccard_comparison(mat_jaccard = mat_jaccard_01,
                                                sobj = sobj_ann,
                                                ref_col = "cell_id_subcluster2",
                                                query_col = "final_level_labels_200",
                                                ref_reduction = "umap")

# doing this on the full dataset is not really usefuls
p_jaccard_azimuth_01

# ---- Jaccard: panhuman annotation vs cluster assignment -- per-cell-type subclusters ----
# every 27_<POP>_subcluster_HarmonySample.rds object predates this panhuman annotation run, so transfer final_level_labels onto it from sobj_ann's metadata by barcode, then compare against a cluster_id column on its own "umap" reduction. Each population is processed as its own block (not a shared function) so any population-specific tweaks stay easy to make later.
# STROMAL/NEU have no manual subcluster annotation yet, so cluster_id is the raw RNA_snn_res.0.4 clustering. ASTRO/IMMUNE/LYM/OLIGO/OPC/VAS already carry Aletta's manual subcluster annotation (data/260813_spinal_subcluster_annotation_aletta.csv, itself keyed on RNA_snn_res.0.4), so cluster_id is cell_id_subcluster2 for those.

sobj_stromal <- readRDS("../../out/object/analysis_R44/27_STROMAL_subcluster_HarmonySample.rds")
# add the metadata from the full objetc after panhuman labelling
sobj_stromal <- AddMetaData(sobj_stromal,metadata = sobj_ann@meta.data[,c("cell_id_subcluster",
                                                                          "cell_id_subcluster2",
                                                                          "full_hierarchical_labels",
                                                                          "final_level_labels",
                                                                          "final_level_confidence",
                                                                          "full_consistent_hierarchy",
                                                                          "azimuth_broad",
                                                                          "azimuth_medium",
                                                                          "azimuth_fine",
                                                                          "azimuth_label",
                                                                          "azimuth_confidence")])

# simplify the annotation also in the subcluster
sobj_stromal <- PrepLabel(sobj_stromal,
                          label_id = "final_level_labels",
                          newid = "final_level_labels_50",
                          cutid = "z-low-freq-labels",
                          cutoff = 50)

# define a cluster_id for the grouping
cluster_id <- "RNA_snn_res.0.4"

meta_stromal <- sobj_stromal@meta.data %>% rownames_to_column("cell_barcode")
mat_jaccard_stromal <- compute_jaccard_matrix(meta_stromal,
                                              id_query_col = "final_level_labels_50",
                                              id_ref_col = cluster_id)

# pull the most frequent annotation per cluster_id based on panhuman, and relabel every cell of that cluster with it
majority_lookup_stromal <- meta_stromal %>%
  group_by(final_level_labels_50,.data[[cluster_id]]) %>%
  summarise(n = n(),.groups = "drop") %>%
  arrange(.data[[cluster_id]],desc(n)) %>%
  group_by(.data[[cluster_id]]) %>%
  slice_max(n = 1,order_by = n,with_ties = F) %>%
  ungroup() %>%
  # select(-n) %>%
  rename(majority_label = final_level_labels_50)

meta_stromal2 <- meta_stromal %>%
  left_join(majority_lookup_stromal, by = cluster_id)

sobj_stromal$majority_label <- meta_stromal2$majority_label

p_jaccard_stromal <- plot_jaccard_comparison(mat_jaccard = mat_jaccard_stromal,
                                              sobj = sobj_stromal,
                                              ref_col = "final_level_labels_50",
                                              query_col = cluster_id,
                                              ref_reduction = "umap",
                                              majority_col = "majority_label")

ggsave(plot = p_jaccard_stromal,
       filename = paste0("../../out/plot/analysis_R45_pixi/00_jaccard_STROMAL_finalLevelLabels50_vs_",cluster_id,".pdf"),
       height = 6, width = 20)

sobj_neu <- readRDS("../../out/object/analysis_R44/27_NEU_subcluster_HarmonySample.rds")
# add the metadata from the full objetc after panhuman labelling
sobj_neu <- AddMetaData(sobj_neu,metadata = sobj_ann@meta.data[,c("cell_id_subcluster",
                                                                    "cell_id_subcluster2",
                                                                    "full_hierarchical_labels",
                                                                    "final_level_labels",
                                                                    "final_level_confidence",
                                                                    "full_consistent_hierarchy",
                                                                    "azimuth_broad",
                                                                    "azimuth_medium",
                                                                    "azimuth_fine",
                                                                    "azimuth_label",
                                                                    "azimuth_confidence")])

# simplify the annotation also in the subcluster
sobj_neu <- PrepLabel(sobj_neu,
                      label_id = "final_level_labels",
                      newid = "final_level_labels_50",
                      cutid = "z-low-freq-labels",
                      cutoff = 50)

# define a cluster_id for the grouping
cluster_id <- "RNA_snn_res.0.4"

meta_neu <- sobj_neu@meta.data %>% rownames_to_column("cell_barcode")
mat_jaccard_neu <- compute_jaccard_matrix(meta_neu,
                                          id_query_col = "final_level_labels_50",
                                          id_ref_col = cluster_id)

# pull the most frequent annotation per cluster_id based on panhuman, and relabel every cell of that cluster with it
majority_lookup_neu <- meta_neu %>%
  group_by(final_level_labels_50,.data[[cluster_id]]) %>%
  summarise(n = n(),.groups = "drop") %>%
  arrange(.data[[cluster_id]],desc(n)) %>%
  group_by(.data[[cluster_id]]) %>%
  slice_max(n = 1,order_by = n,with_ties = F) %>%
  ungroup() %>%
  rename(majority_label = final_level_labels_50)

meta_neu2 <- meta_neu %>%
  left_join(majority_lookup_neu, by = cluster_id)

sobj_neu$majority_label <- meta_neu2$majority_label

p_jaccard_neu <- plot_jaccard_comparison(mat_jaccard = mat_jaccard_neu,
                                          sobj = sobj_neu,
                                          ref_col = "final_level_labels_50",
                                          query_col = cluster_id,
                                          ref_reduction = "umap",
                                          majority_col = "majority_label")

ggsave(plot = p_jaccard_neu,
       filename = paste0("../../out/plot/analysis_R45_pixi/00_jaccard_NEU_finalLevelLabels50_vs_",cluster_id,".pdf"),
       height = 6, width = 20)

# -------------------------------------------------------------------------
# here the one that aletta has provided an annotation

sobj_astro <- readRDS("../../out/object/analysis_R44/27_ASTRO_subcluster_HarmonySample.rds")
# add the metadata from the full objetc after panhuman labelling
sobj_astro <- AddMetaData(sobj_astro,metadata = sobj_ann@meta.data[,c("cell_id_subcluster",
                                                                      "cell_id_subcluster2",
                                                                        "full_hierarchical_labels",
                                                                        "final_level_labels",
                                                                        "final_level_confidence",
                                                                        "full_consistent_hierarchy",
                                                                        "azimuth_broad",
                                                                        "azimuth_medium",
                                                                        "azimuth_fine",
                                                                        "azimuth_label",
                                                                        "azimuth_confidence")])

# simplify the annotation also in the subcluster
sobj_astro <- PrepLabel(sobj_astro,
                        label_id = "final_level_labels",
                        newid = "final_level_labels_50",
                        cutid = "z-low-freq-labels",
                        cutoff = 50)

# define a cluster_id for the grouping -- already-annotated populations use Aletta's manual subcluster call
cluster_id <- "cell_id_subcluster2"

meta_astro <- sobj_astro@meta.data %>% rownames_to_column("cell_barcode")

# some barcodes have no panhuman call (NA final_level_labels_50) or no cluster assignment -- drop them for the jaccard/majority-vote computation rather than let NA pollute the crossing
meta_astro_hc <- meta_astro %>%
  filter(!is.na(.data[[cluster_id]]), !is.na(final_level_labels_50))

mat_jaccard_astro <- compute_jaccard_matrix(meta_astro_hc,
                                            id_query_col = "final_level_labels_50",
                                            id_ref_col = cluster_id)

# pull the most frequent annotation per cluster_id based on panhuman, and relabel every cell of that cluster with it
majority_lookup_astro <- meta_astro_hc %>%
  group_by(final_level_labels_50,.data[[cluster_id]]) %>%
  summarise(n = n(),.groups = "drop") %>%
  arrange(.data[[cluster_id]],desc(n)) %>%
  group_by(.data[[cluster_id]]) %>%
  slice_max(n = 1,order_by = n,with_ties = F) %>%
  ungroup() %>%
  rename(majority_label = final_level_labels_50)

meta_astro2 <- meta_astro %>%
  left_join(majority_lookup_astro, by = cluster_id)

sobj_astro$majority_label <- meta_astro2$majority_label

p_jaccard_astro <- plot_jaccard_comparison(mat_jaccard = mat_jaccard_astro,
                                            sobj = sobj_astro,
                                            ref_col = "final_level_labels_50",
                                            query_col = cluster_id,
                                            ref_reduction = "umap",
                                            majority_col = "majority_label")

ggsave(plot = p_jaccard_astro,
       filename = paste0("../../out/plot/analysis_R45_pixi/00_jaccard_ASTRO_finalLevelLabels50_vs_",cluster_id,".pdf"),
       height = 6, width = 20)

sobj_immune <- readRDS("../../out/object/analysis_R44/27_IMMUNE_subcluster_HarmonySample.rds")
# add the metadata from the full objetc after panhuman labelling
sobj_immune <- AddMetaData(sobj_immune,metadata = sobj_ann@meta.data[,c("cell_id_subcluster",
                                                                          "cell_id_subcluster2",
                                                                          "full_hierarchical_labels",
                                                                          "final_level_labels",
                                                                          "final_level_confidence",
                                                                          "full_consistent_hierarchy",
                                                                          "azimuth_broad",
                                                                          "azimuth_medium",
                                                                          "azimuth_fine",
                                                                          "azimuth_label",
                                                                          "azimuth_confidence")])

# simplify the annotation also in the subcluster
sobj_immune <- PrepLabel(sobj_immune,
                         label_id = "final_level_labels",
                         newid = "final_level_labels_50",
                         cutid = "z-low-freq-labels",
                         cutoff = 50)

# define a cluster_id for the grouping -- already-annotated populations use Aletta's manual subcluster call
cluster_id <- "cell_id_subcluster2"

meta_immune <- sobj_immune@meta.data %>% rownames_to_column("cell_barcode")

# some barcodes have no panhuman call (NA final_level_labels_50) or no cluster assignment -- drop them for the jaccard/majority-vote computation rather than let NA pollute the crossing
meta_immune_hc <- meta_immune %>%
  filter(!is.na(.data[[cluster_id]]), !is.na(final_level_labels_50))

mat_jaccard_immune <- compute_jaccard_matrix(meta_immune_hc,
                                             id_query_col = "final_level_labels_50",
                                             id_ref_col = cluster_id)

# pull the most frequent annotation per cluster_id based on panhuman, and relabel every cell of that cluster with it
majority_lookup_immune <- meta_immune_hc %>%
  group_by(final_level_labels_50,.data[[cluster_id]]) %>%
  summarise(n = n(),.groups = "drop") %>%
  arrange(.data[[cluster_id]],desc(n)) %>%
  group_by(.data[[cluster_id]]) %>%
  slice_max(n = 1,order_by = n,with_ties = F) %>%
  ungroup() %>%
  rename(majority_label = final_level_labels_50)

meta_immune2 <- meta_immune %>%
  left_join(majority_lookup_immune, by = cluster_id)

sobj_immune$majority_label <- meta_immune2$majority_label

p_jaccard_immune <- plot_jaccard_comparison(mat_jaccard = mat_jaccard_immune,
                                             sobj = sobj_immune,
                                             ref_col = "final_level_labels_50",
                                             query_col = cluster_id,
                                             ref_reduction = "umap",
                                             majority_col = "majority_label")

ggsave(plot = p_jaccard_immune,
       filename = paste0("../../out/plot/analysis_R45_pixi/00_jaccard_IMMUNE_finalLevelLabels50_vs_",cluster_id,".pdf"),
       height = 6, width = 20)

sobj_lym <- readRDS("../../out/object/analysis_R44/27_LYM_subcluster_HarmonySample.rds")
# add the metadata from the full objetc after panhuman labelling
sobj_lym <- AddMetaData(sobj_lym,metadata = sobj_ann@meta.data[,c("cell_id_subcluster",
                                                                    "cell_id_subcluster2",
                                                                    "full_hierarchical_labels",
                                                                    "final_level_labels",
                                                                    "final_level_confidence",
                                                                    "full_consistent_hierarchy",
                                                                    "azimuth_broad",
                                                                    "azimuth_medium",
                                                                    "azimuth_fine",
                                                                    "azimuth_label",
                                                                    "azimuth_confidence")])

# simplify the annotation also in the subcluster
sobj_lym <- PrepLabel(sobj_lym,
                      label_id = "final_level_labels",
                      newid = "final_level_labels_50",
                      cutid = "z-low-freq-labels",
                      cutoff = 50)

# define a cluster_id for the grouping -- already-annotated populations use Aletta's manual subcluster call
cluster_id <- "cell_id_subcluster2"

meta_lym <- sobj_lym@meta.data %>% rownames_to_column("cell_barcode")

# some barcodes have no panhuman call (NA final_level_labels_50) or no cluster assignment -- drop them for the jaccard/majority-vote computation rather than let NA pollute the crossing
meta_lym_hc <- meta_lym %>%
  filter(!is.na(.data[[cluster_id]]), !is.na(final_level_labels_50))

mat_jaccard_lym <- compute_jaccard_matrix(meta_lym_hc,
                                          id_query_col = "final_level_labels_50",
                                          id_ref_col = cluster_id)

# pull the most frequent annotation per cluster_id based on panhuman, and relabel every cell of that cluster with it
majority_lookup_lym <- meta_lym_hc %>%
  group_by(final_level_labels_50,.data[[cluster_id]]) %>%
  summarise(n = n(),.groups = "drop") %>%
  arrange(.data[[cluster_id]],desc(n)) %>%
  group_by(.data[[cluster_id]]) %>%
  slice_max(n = 1,order_by = n,with_ties = F) %>%
  ungroup() %>%
  rename(majority_label = final_level_labels_50)

meta_lym2 <- meta_lym %>%
  left_join(majority_lookup_lym, by = cluster_id)

sobj_lym$majority_label <- meta_lym2$majority_label

p_jaccard_lym <- plot_jaccard_comparison(mat_jaccard = mat_jaccard_lym,
                                          sobj = sobj_lym,
                                          ref_col = "final_level_labels_50",
                                          query_col = cluster_id,
                                          ref_reduction = "umap",
                                          majority_col = "majority_label")

ggsave(plot = p_jaccard_lym,
       filename = paste0("../../out/plot/analysis_R45_pixi/00_jaccard_LYM_finalLevelLabels50_vs_",cluster_id,".pdf"),
       height = 6, width = 20)

sobj_oligo <- readRDS("../../out/object/analysis_R44/27_OLIGO_subcluster_HarmonySample.rds")
# add the metadata from the full objetc after panhuman labelling
sobj_oligo <- AddMetaData(sobj_oligo,metadata = sobj_ann@meta.data[,c("cell_id_subcluster",
                                                                        "cell_id_subcluster2",
                                                                        "full_hierarchical_labels",
                                                                        "final_level_labels",
                                                                        "final_level_confidence",
                                                                        "full_consistent_hierarchy",
                                                                        "azimuth_broad",
                                                                        "azimuth_medium",
                                                                        "azimuth_fine",
                                                                        "azimuth_label",
                                                                        "azimuth_confidence")])

# simplify the annotation also in the subcluster
sobj_oligo <- PrepLabel(sobj_oligo,
                        label_id = "final_level_labels",
                        newid = "final_level_labels_50",
                        cutid = "z-low-freq-labels",
                        cutoff = 50)

# define a cluster_id for the grouping -- already-annotated populations use Aletta's manual subcluster call
cluster_id <- "cell_id_subcluster2"

meta_oligo <- sobj_oligo@meta.data %>% rownames_to_column("cell_barcode")

# some barcodes have no panhuman call (NA final_level_labels_50) or no cluster assignment -- drop them for the jaccard/majority-vote computation rather than let NA pollute the crossing
meta_oligo_hc <- meta_oligo %>%
  filter(!is.na(.data[[cluster_id]]), !is.na(final_level_labels_50))

mat_jaccard_oligo <- compute_jaccard_matrix(meta_oligo_hc,
                                            id_query_col = "final_level_labels_50",
                                            id_ref_col = cluster_id)

# pull the most frequent annotation per cluster_id based on panhuman, and relabel every cell of that cluster with it
majority_lookup_oligo <- meta_oligo_hc %>%
  group_by(final_level_labels_50,.data[[cluster_id]]) %>%
  summarise(n = n(),.groups = "drop") %>%
  arrange(.data[[cluster_id]],desc(n)) %>%
  group_by(.data[[cluster_id]]) %>%
  slice_max(n = 1,order_by = n,with_ties = F) %>%
  ungroup() %>%
  rename(majority_label = final_level_labels_50)

meta_oligo2 <- meta_oligo %>%
  left_join(majority_lookup_oligo, by = cluster_id)

sobj_oligo$majority_label <- meta_oligo2$majority_label

p_jaccard_oligo <- plot_jaccard_comparison(mat_jaccard = mat_jaccard_oligo,
                                            sobj = sobj_oligo,
                                            ref_col = "final_level_labels_50",
                                            query_col = cluster_id,
                                            ref_reduction = "umap",
                                            majority_col = "majority_label")

ggsave(plot = p_jaccard_oligo,
       filename = paste0("../../out/plot/analysis_R45_pixi/00_jaccard_OLIGO_finalLevelLabels50_vs_",cluster_id,".pdf"),
       height = 6, width = 20)

sobj_opc <- readRDS("../../out/object/analysis_R44/27_OPC_subcluster_HarmonySample.rds")
# add the metadata from the full objetc after panhuman labelling
sobj_opc <- AddMetaData(sobj_opc,metadata = sobj_ann@meta.data[,c("cell_id_subcluster",
                                                                    "cell_id_subcluster2",
                                                                    "full_hierarchical_labels",
                                                                    "final_level_labels",
                                                                    "final_level_confidence",
                                                                    "full_consistent_hierarchy",
                                                                    "azimuth_broad",
                                                                    "azimuth_medium",
                                                                    "azimuth_fine",
                                                                    "azimuth_label",
                                                                    "azimuth_confidence")])

# simplify the annotation also in the subcluster
sobj_opc <- PrepLabel(sobj_opc,
                      label_id = "final_level_labels",
                      newid = "final_level_labels_50",
                      cutid = "z-low-freq-labels",
                      cutoff = 50)

# define a cluster_id for the grouping -- already-annotated populations use Aletta's manual subcluster call
cluster_id <- "cell_id_subcluster2"

meta_opc <- sobj_opc@meta.data %>% rownames_to_column("cell_barcode")

# some barcodes have no panhuman call (NA final_level_labels_50) or no cluster assignment -- drop them for the jaccard/majority-vote computation rather than let NA pollute the crossing
meta_opc_hc <- meta_opc %>%
  filter(!is.na(.data[[cluster_id]]), !is.na(final_level_labels_50))

mat_jaccard_opc <- compute_jaccard_matrix(meta_opc_hc,
                                          id_query_col = "final_level_labels_50",
                                          id_ref_col = cluster_id)

# pull the most frequent annotation per cluster_id based on panhuman, and relabel every cell of that cluster with it
majority_lookup_opc <- meta_opc_hc %>%
  group_by(final_level_labels_50,.data[[cluster_id]]) %>%
  summarise(n = n(),.groups = "drop") %>%
  arrange(.data[[cluster_id]],desc(n)) %>%
  group_by(.data[[cluster_id]]) %>%
  slice_max(n = 1,order_by = n,with_ties = F) %>%
  ungroup() %>%
  rename(majority_label = final_level_labels_50)

meta_opc2 <- meta_opc %>%
  left_join(majority_lookup_opc, by = cluster_id)

sobj_opc$majority_label <- meta_opc2$majority_label

p_jaccard_opc <- plot_jaccard_comparison(mat_jaccard = mat_jaccard_opc,
                                          sobj = sobj_opc,
                                          ref_col = "final_level_labels_50",
                                          query_col = cluster_id,
                                          ref_reduction = "umap",
                                          majority_col = "majority_label")

ggsave(plot = p_jaccard_opc,
       filename = paste0("../../out/plot/analysis_R45_pixi/00_jaccard_OPC_finalLevelLabels50_vs_",cluster_id,".pdf"),
       height = 6, width = 20)

sobj_vas <- readRDS("../../out/object/analysis_R44/27_VAS_subcluster_HarmonySample.rds")
# add the metadata from the full objetc after panhuman labelling
sobj_vas <- AddMetaData(sobj_vas,metadata = sobj_ann@meta.data[,c("cell_id_subcluster",
                                                                    "cell_id_subcluster2",
                                                                    "full_hierarchical_labels",
                                                                    "final_level_labels",
                                                                    "final_level_confidence",
                                                                    "full_consistent_hierarchy",
                                                                    "azimuth_broad",
                                                                    "azimuth_medium",
                                                                    "azimuth_fine",
                                                                    "azimuth_label",
                                                                    "azimuth_confidence")])

# simplify the annotation also in the subcluster
sobj_vas <- PrepLabel(sobj_vas,
                      label_id = "final_level_labels",
                      newid = "final_level_labels_50",
                      cutid = "z-low-freq-labels",
                      cutoff = 50)

# define a cluster_id for the grouping -- already-annotated populations use Aletta's manual subcluster call
cluster_id <- "cell_id_subcluster2"

meta_vas <- sobj_vas@meta.data %>% rownames_to_column("cell_barcode")

# some barcodes have no panhuman call (NA final_level_labels_50) or no cluster assignment -- drop them for the jaccard/majority-vote computation rather than let NA pollute the crossing
meta_vas_hc <- meta_vas %>%
  filter(!is.na(.data[[cluster_id]]), !is.na(final_level_labels_50))

mat_jaccard_vas <- compute_jaccard_matrix(meta_vas_hc,
                                          id_query_col = "final_level_labels_50",
                                          id_ref_col = cluster_id)

# pull the most frequent annotation per cluster_id based on panhuman, and relabel every cell of that cluster with it
majority_lookup_vas <- meta_vas_hc %>%
  group_by(final_level_labels_50,.data[[cluster_id]]) %>%
  summarise(n = n(),.groups = "drop") %>%
  arrange(.data[[cluster_id]],desc(n)) %>%
  group_by(.data[[cluster_id]]) %>%
  slice_max(n = 1,order_by = n,with_ties = F) %>%
  ungroup() %>%
  rename(majority_label = final_level_labels_50)

meta_vas2 <- meta_vas %>%
  left_join(majority_lookup_vas, by = cluster_id)

sobj_vas$majority_label <- meta_vas2$majority_label

p_jaccard_vas <- plot_jaccard_comparison(mat_jaccard = mat_jaccard_vas,
                                          sobj = sobj_vas,
                                          ref_col = "final_level_labels_50",
                                          query_col = cluster_id,
                                          ref_reduction = "umap",
                                          majority_col = "majority_label")

ggsave(plot = p_jaccard_vas,
       filename = paste0("../../out/plot/analysis_R45_pixi/00_jaccard_VAS_finalLevelLabels50_vs_",cluster_id,".pdf"),
       height = 6, width = 20)

