# AIM ---------------------------------------------------------------------
# Differential abundance analysis MS vs CTRL (diagnosis_short) with miloR (https://github.com/MarioniLab/miloR) on the final consensus object generated in 03_apply_filter_annotate_fullObject.R.
# 1. build the Milo object from the Seurat object, re-using the SNN graph computed by Seurat on the harmony embedding (dims 1:30)
# 2. define the neighbourhoods, count the cells per neighbourhood and sample (sample_id)
# 3. test the differential abundance with the design ~ diagnosis_short (CTRL is the reference, positive logFC = enriched in MS)
# 4. annotate the neighbourhoods with cell_id and cell_id_subcluster3 and plot the DA per annotation level
# readapted from scRNA_SOD_3rd/analysis/R_analysis_R45_pixi/scr/general/29_miloR_SOD_vs_CT_differentialAbundance.R
# NOTICE:
# - the graph is built on the harmony embedding (not on the PCA), so the harmony embedding is the one used by makeNhoods for the refinement of the index cells
# - the SNN graph of Seurat stores the Jaccard weights. Following https://github.com/MarioniLab/miloR/issues/351#issuecomment-2462720639 the graph is converted to a binary dgCMatrix before buildFromAdjacency, to avoid the RAM spike on big objects (miloR only uses the graph topology)
# - the object is NOT converted with as.SingleCellExperiment: the SCE is built by hand with only the RNA logcounts and the reductions, to keep the RAM low
# - the design is ~ diagnosis_short only (no covariate, e.g. location), as the comparison is run on the sample as unit of replication

# setup -------------------------------------------------------------------
# miloR is not in the renv library yet, install it once
# BiocManager::install("miloR")
# renv::snapshot()

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(patchwork)
library(miloR)
library(SingleCellExperiment)
# plotUMAP is from scater, not miloR
library(scater)
library(ggbeeswarm)

# parameters --------------------------------------------------------------
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 2000 * 1024^2)

# comparison: the first level is the reference, positive logFC = enriched in the second level
var_contrast <- "diagnosis_short"
levels_contrast <- c("CTRL", "MS")

# the sample identifier (unit of replication of the DA test)
id_sample <- "sample_id"

# annotation levels to test
list_level <- c("cell_id", "cell_id_subcluster3")

# name of the Seurat graph to re-use. The default graph name of FindNeighbors on the RNA assay is "RNA_snn". If it is not in the object, the graph is recomputed on the harmony embedding with the k_graph and d_graph below
name_graph <- "RNA_snn"

# k of the Seurat graph (FindNeighbors default k.param = 20, change it if the graph was computed with a different value) and number of harmony dimensions used to compute the graph
k_graph <- 20
d_graph <- 30

# proportion of cells sampled as index of the neighbourhoods
prop_nhood <- 0.1

# significance threshold on the spatial FDR
fdr_thr <- 0.05

# minimum fraction of cells with the same label to assign the label to a neighbourhood
frac_thr <- 0.7

folder_table <- "../../out/table/analysis_R45_pixi/08_miloR/"
folder_plot <- "../../out/plot/analysis_R45_pixi/08_miloR/"
folder_object <- "../../out/object/analysis_R45_pixi/08_miloR/"
dir.create(folder_table, recursive = TRUE, showWarnings = FALSE)
dir.create(folder_plot, recursive = TRUE, showWarnings = FALSE)
dir.create(folder_object, recursive = TRUE, showWarnings = FALSE)

contrast_label <- paste0(levels_contrast[2], "_vs_", levels_contrast[1])

# read in the data --------------------------------------------------------
sobj <- readRDS("../../out/object/analysis_R45_pixi/03_sobj_integrated_cleanup_manualAnnotation_subclusterFiltered_v2.rds")

# make sure the correct default assay is loaded, should be RNA
DefaultAssay(sobj) <- "RNA"

# join the layers in case they are split per sample
if(length(Layers(sobj, search = "data")) > 1){
  sobj <- JoinLayers(sobj)
}

# confirm the graph and the reductions are in the object
Graphs(sobj)
Reductions(sobj)
stopifnot(all(c("harmony", "umap") %in% Reductions(sobj)))
stopifnot(all(c(var_contrast, id_sample, list_level) %in% colnames(sobj@meta.data)))

# the design of the test: samples and condition. every sample has to belong to a single level of the contrast
sobj@meta.data %>%
  group_by(.data[[var_contrast]], .data[[id_sample]]) %>%
  summarise(n = n(), .groups = "drop") %>%
  as.data.frame() %>%
  print()

stopifnot(all(sobj@meta.data[[var_contrast]] %in% levels_contrast))
stopifnot(all(rowSums(table(sobj@meta.data[[id_sample]], sobj@meta.data[[var_contrast]]) > 0) == 1))

# the missing labels would break annotateNhoods, make them an explicit label
for(level in list_level){
  sobj@meta.data[[level]] <- as.character(sobj@meta.data[[level]])
  sobj@meta.data[[level]][is.na(sobj@meta.data[[level]])] <- "NA"
}

# graph -------------------------------------------------------------------
# SNN graph computed by Seurat on the harmony embedding, it is the KNN graph on which the neighbourhoods are defined
if(!name_graph %in% Graphs(sobj)){
  message(name_graph, " not in the object: recompute the graph on harmony, dims 1:", d_graph, ", k ", k_graph)
  sobj <- FindNeighbors(sobj, reduction = "harmony", dims = 1:d_graph, k.param = k_graph, graph.name = c("milo_nn", "milo_snn"), verbose = FALSE)
  name_graph <- "milo_snn"
}
graph_obj <- sobj@graphs[[name_graph]]

# the graph has to follow the order of the cells in the object, otherwise the neighbourhoods would be assigned to the wrong cells
stopifnot(identical(rownames(graph_obj), colnames(sobj)))

# size of the graph as stored by Seurat (Graph class, with the Jaccard weights)
print(object.size(graph_obj), units = "Mb")

# convert it to a plain sparse matrix and check there are no NAs in the weights
graph_obj_fix <- as(graph_obj, "dgCMatrix")
sum(is.na(graph_obj_fix@x))

# binarize it: every edge gets weight 1. miloR only needs the topology of the graph (who is neighbour of whom)
graph_obj_fix@x <- rep(1, length(graph_obj_fix@x))
print(object.size(graph_obj_fix), units = "Mb")

# convert the adjacency matrix to the igraph object used by miloR, k is the k used to build the Seurat graph
graph_milo <- buildFromAdjacency(graph_obj_fix, k = k_graph, is.binary = TRUE)
print(object.size(graph_milo), units = "Mb")

# milo object -------------------------------------------------------------
# lightweight SCE built by hand instead of as.SingleCellExperiment (which would copy the whole Seurat object):
# - logcounts: RNA normalized data, not used by the DA test but needed if the neighbourhood markers (findNhoodMarkers) are computed later
# - colData: all the metadata of the Seurat object
# - reducedDims: names in upper case, as in as.SingleCellExperiment. HARMONY is used by makeNhoods, UMAP by the plots
list_reduced <- list(HARMONY = Embeddings(sobj, "harmony"),
                     UMAP = Embeddings(sobj, "umap"))
if("pca" %in% Reductions(sobj)){
  list_reduced$PCA <- Embeddings(sobj, "pca")
}

sce <- SingleCellExperiment(assays = list(logcounts = GetAssayData(sobj, assay = "RNA", layer = "data")),
                            colData = sobj@meta.data,
                            reducedDims = list_reduced)

# Milo extends the SingleCellExperiment class with the slots for the graph and the neighbourhoods
milo <- Milo(sce)

# add the graph computed in Seurat to the object, instead of recomputing it with buildGraph
miloR::graph(milo) <- miloR::graph(graph_milo)

# free the RAM: the Seurat object and the intermediate graphs are not needed anymore
rm(sobj, sce, graph_obj, graph_obj_fix, graph_milo)
gc()

# confirm the embedding and the metadata are correctly loaded
p_umap_check <- plotUMAP(milo, colour_by = list_level[1])
colData(milo)

# neighbourhoods ----------------------------------------------------------
# sample the index cells of the neighbourhoods: a random proportion (prop) of the cells is refined (refined = TRUE, recommended) by moving each index to the cell closest to the centre of its neighbours, in the embedding given by reduced_dims and d. Each neighbourhood is the index cell plus its neighbours in the graph.
# the refinement uses the same embedding used to build the graph (harmony, 30 dims). The random sampling is seeded to be reproducible
set.seed(123)
milo <- makeNhoods(milo, prop = prop_nhood, k = k_graph, d = d_graph, refined = TRUE, reduced_dims = "HARMONY")

# empirically the distribution of the size should peak between 50 and 100 cells. Otherwise rerun makeNhoods increasing k and/or prop (k has to be consistent with the graph, so it is better to start from the prop)
p_size <- plotNhoodSizeHist(milo)
p_size

ggsave(paste0(folder_plot, "08_miloR_nhoodSizeHist.pdf"), plot = p_size, width = 6, height = 4)

# count the cells per neighbourhood and sample ----------------------------
# the result is a matrix neighbourhoods x samples with the number of cells of each sample in each neighbourhood, this is the input of the DA test
milo <- countCells(milo, meta.data = data.frame(colData(milo)), samples = id_sample)
head(nhoodCounts(milo))

# differential abundance testing ------------------------------------------
# design: one row per sample, in the order of the columns of the count matrix. The first level of the contrast is the reference
design_milo <- data.frame(colData(milo))[, c(id_sample, var_contrast)] %>%
  distinct() %>%
  mutate(!!var_contrast := factor(.data[[var_contrast]], levels = levels_contrast))
rownames(design_milo) <- design_milo[[id_sample]]
design_milo <- design_milo[colnames(nhoodCounts(milo)), , drop = FALSE]
design_milo

# NB GLM (edgeR) on the neighbourhood counts, one test per neighbourhood. The p values are corrected with the spatial FDR, which weights each test by the connectivity of the neighbourhood to account for their overlap. "graph-overlap" computes the weights from the overlap of the neighbourhoods in the graph, so there is no need of calcNhoodDistance.
# the output has one row per neighbourhood: logFC (MS vs CTRL), PValue, FDR and SpatialFDR
da_results <- testNhoods(milo, design = as.formula(paste("~", var_contrast)), design.df = design_milo, fdr.weighting = "graph-overlap")

da_results %>%
  arrange(SpatialFDR) %>%
  head()

# diagnostic plots: uncorrected p values (the test should be balanced) and volcano
p_pval <- ggplot(da_results, aes(PValue)) +
  geom_histogram(bins = 50) +
  theme_bw()

p_volcano <- ggplot(da_results, aes(logFC, -log10(SpatialFDR))) +
  geom_point(shape = 1, alpha = 0.5) +
  geom_hline(yintercept = -log10(fdr_thr), linetype = "dashed", col = "gray") +
  theme_bw() +
  theme(strip.background = element_blank())

ggsave(paste0(folder_plot, "08_miloR_diagnostic_pvalue_volcano.pdf"), plot = p_pval | p_volcano, width = 12, height = 5)

# neighbourhood graph on the UMAP -----------------------------------------
# abstracted graph: each node is a neighbourhood (placed on the UMAP at the position of its index cell), each edge is the number of cells shared by two neighbourhoods
milo <- buildNhoodGraph(milo)

# UMAP of the single cells, colored by condition, to compare with the neighbourhood graph
p_umap <- plotUMAP(milo, colour_by = var_contrast)

# same UMAP split by condition (one panel per level). other_fields adds the column to the plot data so it can be used for the facets
p_umap_split <- plotUMAP(milo, colour_by = var_contrast, other_fields = var_contrast) +
  facet_wrap(as.formula(paste("~", var_contrast)))

# neighbourhood graph colored by the logFC. blue = enriched in the reference, red = enriched in the second level
scale_logfc <- scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, na.value = "grey50", name = paste0("logFC\n", levels_contrast[2], " vs ", levels_contrast[1]))

# only the significant neighbourhoods are colored (SpatialFDR < fdr_thr), the others are set to 0 (white) by plotNhoodGraphDA
p_nhood <- plotNhoodGraphDA(milo, da_results, layout = "UMAP", alpha = fdr_thr) +
  scale_logfc +
  ggtitle(paste0("Significant neighbourhoods only (SpatialFDR < ", fdr_thr, ")"))

# continuous logFC for all the neighbourhoods, not clipped by the significance: plotNhoodGraphDA sets the logFC to 0 when the SpatialFDR is above alpha, so the SpatialFDR is set to 0 in a copy of the table only for this plot (da_results is not modified)
p_nhood_cont <- plotNhoodGraphDA(milo, da_results %>% mutate(SpatialFDR = 0), layout = "UMAP", alpha = fdr_thr) +
  scale_logfc +
  ggtitle("All neighbourhoods (logFC, no significance threshold)")

# top row: single cells by condition (all together, split); bottom row: neighbourhood graph (significant only, continuous)
p_umap_all <- (p_umap | p_umap_split) / (p_nhood | p_nhood_cont)

ggsave(paste0(folder_plot, "08_miloR_nhoodGraphDA_UMAP.pdf"), plot = p_umap_all, width = 20, height = 14)

# annotation of the neighbourhoods ----------------------------------------
# for each level of the annotation:
# 1. annotateNhoods assigns to each neighbourhood the most abundant label among its cells, and saves the fraction of cells with that label (<level>_fraction)
# 2. neighbourhoods with a fraction below frac_thr are labelled "Mixed" (e.g. transitional states, or boundaries between cell types)
# 3. table of the DA results and count of the significant neighbourhoods (SpatialFDR < fdr_thr) per label, split in enriched in the second level (logFC > 0) and in the reference (logFC < 0)
# 4. beeswarm plot of the logFC of the neighbourhoods per label (labels ordered by the median logFC, significant neighbourhoods in red), next to the histogram of the label fraction used to choose frac_thr
# level <- "cell_id"
list_da_annotated <- lapply(list_level, function(level){
  da_level <- annotateNhoods(milo, da_results, coldata_col = level)

  # distribution of the fraction of cells with the label in each neighbourhood
  p_frac <- ggplot(da_level, aes(.data[[paste0(level, "_fraction")]])) +
    geom_histogram(bins = 50) +
    theme_bw() +
    labs(x = paste0(level, "_fraction"))

  # neighbourhoods with a mix of labels
  da_level <- da_level %>%
    mutate(!!level := ifelse(.data[[paste0(level, "_fraction")]] < frac_thr, "Mixed", .data[[level]]))

  da_level %>%
    write_tsv(paste0(folder_table, "08_miloR_DA_", contrast_label, "_", level, ".tsv"))

  # number of significant neighbourhoods per label
  da_level %>%
    group_by(.data[[level]]) %>%
    summarise(n_nhood = n(),
              n_up = sum(SpatialFDR < fdr_thr & logFC > 0),
              n_down = sum(SpatialFDR < fdr_thr & logFC < 0),
              .groups = "drop") %>%
    write_tsv(paste0(folder_table, "08_miloR_DA_summary_", contrast_label, "_", level, ".tsv"))

  # beeswarm of the logFC per label, labels ordered by the median logFC
  order_level <- da_level %>%
    group_by(.data[[level]]) %>%
    summarise(med = median(logFC), .groups = "drop") %>%
    arrange(med) %>%
    pull(.data[[level]])

  p_bee <- da_level %>%
    mutate(label = factor(.data[[level]], levels = order_level),
           color = ifelse(SpatialFDR < fdr_thr, "sig", "n.s.")) %>%
    arrange(color) %>%
    ggplot(aes(x = logFC, y = label, color = color)) +
    geom_quasirandom(alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", col = "gray") +
    scale_color_manual(values = c("n.s." = "gray", "sig" = "red")) +
    theme_bw() +
    labs(title = paste0("miloR DA ", levels_contrast[2], " vs ", levels_contrast[1], " - ", level),
         subtitle = paste0("SpatialFDR < ", fdr_thr, ", positive logFC = enriched in ", levels_contrast[2]),
         y = NULL)

  # the height follows the number of labels
  ggsave(paste0(folder_plot, "08_miloR_beeswarm_", level, ".pdf"), plot = p_bee | p_frac,
         width = 14, height = max(6, 0.25 * length(order_level) + 2), limitsize = FALSE)

  da_level
})
names(list_da_annotated) <- list_level

# save --------------------------------------------------------------------
# full DA table (one row per neighbourhood, unannotated) and the milo object to re-do the plots without recomputing the neighbourhoods
da_results %>%
  write_tsv(paste0(folder_table, "08_miloR_DA_", contrast_label, ".tsv"))

saveRDS(milo, paste0(folder_object, "08_miloR_object_", contrast_label, ".rds"))

message("\nAll done.")
