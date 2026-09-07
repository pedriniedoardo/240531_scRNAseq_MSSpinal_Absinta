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

# build a Jaccard-score Heatmap paired with the two side-by-side DimPlots it compares (ref_col vs query_col), so the same figure can be produced for multiple annotation/cluster comparisons.
# ref_col and query_col double as both the DimPlot group.by and the heatmap column/row titles.
# query_reduction defaults to ref_reduction, but pass both explicitly when the two annotations live natively in different reductions (e.g. an original annotation on "umap" vs an azimuth-derived one on "azimuth_umap").
plot_jaccard_comparison <- function(mat_jaccard,
                                     sobj,
                                     ref_col,
                                     query_col,
                                     ref_reduction,
                                     query_reduction = ref_reduction,
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

  # grid.grabExpr(draw(x)) turns the ComplexHeatmap output into a grob patchwork can lay out alongside the ggplots
  (p_ref + p_query + grid.grabExpr(draw(ht))) +
    plot_layout(widths = c(1, 1, 3))
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

# ---- Jaccard similarity between cluster assignments and cell type labels ----
# use the Jaccard score to measure the cross-cluster similarity per cell (how similar are the clusters from the query compared to the annotation derived from the reference)

# note: sobj already has a metadata column literally named "barcode", so use a
# different name here for the rownames-to-column step to avoid a collision.
meta_hc <- sobj_ann@meta.data %>% rownames_to_column("cell_barcode")

# Heatmap A: azimuth_embed-based clusters vs. Pan-human Azimuth cell type labels

# compare the broad cell type we provided with the finer annotation produced by azimuth
mat_jaccard_01 <- compute_jaccard_matrix(meta_hc, "cell_type", "final_level_labels")
# compare the broard annotation with the congruent annotation priduced by azimuth
mat_jaccard_02 <- compute_jaccard_matrix(meta_hc, "cell_type", "azimuth_fine")
# compare the fine annotation produced by azimuth with the clustering calculated on the azimuth embedding
mat_jaccard_03 <- compute_jaccard_matrix(meta_hc, "final_level_labels", "azimuth_snn_res.0.4")
# compare the fine annotation produced by azimuth with the finer annotatoin priduced by aletta
mat_jaccard_04 <- compute_jaccard_matrix(meta_hc, "final_level_labels", "cell_id_subcluster")
mat_jaccard_05 <- compute_jaccard_matrix(meta_hc, "final_level_labels", "cell_id_subcluster2")


p_jaccard_azimuth <- plot_jaccard_comparison(mat_jaccard = mat_jaccard_03,
                                              sobj = sobj_ann,
                                              ref_col = "final_level_labels",
                                              query_col = "azimuth_snn_res.0.4",
                                              ref_reduction = "azimuth_umap")

p_jaccard_azimuth

ggsave(plot = p_jaccard_azimuth,
       filename = "../../out/plot/analysis_R45_pixi/00_jaccard_azimuthClusters_vs_finalLevelLabels.pdf",
       height = 6, width = 16)

# other comparisons (mat_jaccard_01/02/04/05) can be plotted the same way,
# e.g.:
# plot_jaccard_comparison(mat_jaccard_01, sobj_ann, "cell_type", "final_level_labels",
#                          ref_reduction = "umap", query_reduction = "azimuth_umap")

# ---- Save the final benchmarked object --------------------------------------
saveRDS(sobj_ann, "../../out/object/analysis_R45_pixi/00_sobj_ann_benchmarked.rds")
