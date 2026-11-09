# AIM ---------------------------------------------------------------------
# custom plots on top of the miloR results of the design adjusted for location and sex (08c_miloR_DA_adjusted_fullObject.R: ~ location + sex + diagnosis_short). The milo object saved by 08_miloR_DA_fullObject.R and the annotated adjusted DA table are loaded, nothing is recomputed, the Seurat object is NOT needed.
# plots for one macro population (pop, default ASTRO), one panel per cell_id_subcluster3:
# 1. overview: all the cells of the population colored by subcluster, with the name of the subcluster on the cells
# 2. one panel per subcluster: the cells of the subcluster highlighted in red on the grey background of the whole population, with the name of the subcluster on the highlighted cells (and in the title, with the number of cells)
# 3. one panel per subcluster: the miloR neighbourhoods annotated with that subcluster, drawn at the position of their index cell and colored by logFC (positive = enriched in MS), the significant ones (SpatialFDR < fdr_thr) have a black outline
# NOTICE:
# - the milo object comes from 08_miloR_DA_fullObject.R (1.1 GB, load time about a minute); it holds the neighbourhoods and the counts, the DA results of the adjusted design are in the table of 08c (matched by Nhood)
# - nhoodIndex of a Milo object stores the POSITION of the index cell (column of the object), not the barcode
# - the DA table is the annotated one: cell_id_subcluster3 is the majority label of the neighbourhood, neighbourhoods below frac_thr were already relabelled "Mixed" in 08c_miloR_DA_adjusted_fullObject.R, so they are not in the subclusters of the population

# libraries ---------------------------------------------------------------
library(tidyverse)
library(patchwork)
library(miloR)
library(SingleCellExperiment)
library(cowplot)

# parameters --------------------------------------------------------------
# macro population to plot: the cells with cell_id == pop (all the subcluster3 labels of the cells are the ones of the population, e.g. ASTRO|...)
pop <- "ASTRO"
level_cell <- "cell_id"
level_sub <- "cell_id_subcluster3"

# number of columns of the panels
ncol_panel <- 3

# significance threshold on the spatial FDR (same as 08c_miloR_DA_adjusted_fullObject.R)
fdr_thr <- 0.05

contrast_label <- "MS_vs_CTRL"
levels_contrast <- c("CTRL", "MS")

# DA results and plots: design adjusted for location and sex (08c_miloR_DA_adjusted_fullObject.R)
folder_table <- "../../out/table/analysis_R45_pixi/08_miloR_adjusted/"
folder_plot <- "../../out/plot/analysis_R45_pixi/08_miloR_adjusted/"
# milo object (graph, neighbourhoods, counts) saved by 08_miloR_DA_fullObject.R. 08c does not save it again: the neighbourhoods are the same, only the test changes
folder_object <- "../../out/object/analysis_R45_pixi/08_miloR/"

# read in the data --------------------------------------------------------
milo <- readRDS(paste0(folder_object, "08_miloR_object_", contrast_label, ".rds"))

# annotated DA table of the adjusted design (one row per neighbourhood, from 08c_miloR_DA_adjusted_fullObject.R)
da_sub <- read_tsv(paste0(folder_table, "08c_miloR_DA_adjusted_", contrast_label, "_", level_sub, ".tsv"), show_col_types = FALSE)

# coordinates ------------------------------------------------------------
umap <- reducedDim(milo, "UMAP")

# one row per cell, only the columns needed for the plots
df_all <- tibble(barcode = colnames(milo),
                 UMAP1 = umap[, 1],
                 UMAP2 = umap[, 2],
                 cell_id = as.character(colData(milo)[[level_cell]]),
                 subcluster = as.character(colData(milo)[[level_sub]]))

# index cell of each neighbourhood: nhoodIndex is a list of positions in the object, the neighbourhood n is the n-th element and matches the column Nhood of the DA table
idx_nhood <- unlist(nhoodIndex(milo))
stopifnot(length(idx_nhood) == nrow(da_sub),
          all(sort(da_sub$Nhood) == seq_along(idx_nhood)))

df_nhood <- tibble(Nhood = seq_along(idx_nhood),
                   UMAP1 = umap[idx_nhood, 1],
                   UMAP2 = umap[idx_nhood, 2]) %>%
  left_join(da_sub, by = "Nhood") %>%
  mutate(sig = SpatialFDR < fdr_thr)

# free the RAM: the milo object is not needed anymore
rm(milo, umap)
gc()

# population --------------------------------------------------------------
df_pop <- df_all %>%
  filter(cell_id == pop)

sub_levels <- sort(unique(df_pop$subcluster))
sub_levels

# confirm the number of cells per subcluster
table(df_pop$subcluster)

# position of the label: median of the coordinates of the cells of each subcluster
df_lab <- df_pop %>%
  group_by(subcluster) %>%
  summarise(UMAP1 = median(UMAP1),
            UMAP2 = median(UMAP2),
            n = n(),
            .groups = "drop")

# 1. overview: all the subclusters, named on the cells --------------------
p_overview <- df_pop %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = subcluster)) +
  geom_point(size = 0.2) +
  ggrepel::geom_label_repel(data = df_lab, aes(label = subcluster), size = 3, fill = alpha("white", 0.7), label.size = 0, show.legend = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_bw() +
  labs(title = paste0(pop, ": ", level_sub))

# ggsave(paste0(folder_plot, "08d_miloR_", pop, "_overview_subcluster.pdf"), plot = p_overview, width = 10, height = 7)

# 2. one panel per subcluster: highlighted cells, named ---------------------
# sb <- sub_levels[1]
list_plot_hl <- map(sub_levels, function(sb){
  df_hl <- df_pop %>%
    filter(subcluster == sb)

  ggplot() +
    # background: all the cells of the population
    geom_point(data = df_pop, aes(x = UMAP1, y = UMAP2), color = "grey85", size = 0.3) +
    # selected cells
    geom_point(data = df_hl, aes(x = UMAP1, y = UMAP2), color = "red", size = 0.1) +
    # name of the subset on the selected cells
    # geom_label(data = df_lab %>% filter(subcluster == sb), aes(x = UMAP1, y = UMAP2, label = subcluster),
    #            size = 3, fill = alpha("white", 0.7), label.size = 0) +
    theme_cowplot() +
    labs(title = paste0(sb, " (n = ", nrow(df_hl), ")"))
})

p_hl <- wrap_plots(list_plot_hl, ncol = ncol_panel)

ggsave(paste0(folder_plot, "08d_miloR_", pop, "_highlight_subcluster.pdf"), plot = p_hl,
       width = 5 * ncol_panel, height = 4.5 * ceiling(length(sub_levels) / ncol_panel), limitsize = FALSE)

# 3. one panel per subcluster: miloR neighbourhoods colored by logFC --------
# neighbourhoods annotated with one of the subclusters of the population
df_nhood_pop <- df_nhood %>%
  filter(.data[[level_sub]] %in% sub_levels)

# same limits in all the panels, symmetric around 0, so the colors are comparable
# lim_logfc <- max(abs(df_nhood_pop$logFC))
lim_logfc <- 3

# sb <- sub_levels[1]
list_plot_nhood <- map(sub_levels, function(sb){
  df_nh <- df_nhood_pop %>%
    filter(.data[[level_sub]] == sb) %>%
    arrange(logFC)

  ggplot() +
    geom_point(data = df_pop, aes(x = UMAP1, y = UMAP2), color = "grey85", size = 0.3) +
    # significant neighbourhoods on top, with the black outline
    geom_point(data = df_nh %>% arrange(sig), aes(x = UMAP1, y = UMAP2, fill = logFC, color = sig), shape = 21, size = 2, stroke = 0.1,colour = "gray",alpha=0.8) +
    # geom_label(data = df_lab %>% filter(subcluster == sb), aes(x = UMAP1, y = UMAP2, label = subcluster),
    #            size = 3, fill = alpha("white", 0.7), label.size = 0) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-lim_logfc, lim_logfc),
                         name = paste0("logFC\n", levels_contrast[2], " vs ", levels_contrast[1]),oob = scales::squish) +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "grey60"), guide = "none") +
    theme_cowplot() +
    labs(title = paste0(sb, "\n (", nrow(df_nh), " nhoods, ", sum(df_nh$sig), " significant)"))
})

p_nhood <- wrap_plots(list_plot_nhood, ncol = ncol_panel) +
  plot_layout(guides = "collect") +
  plot_annotation(title = paste0("miloR neighbourhoods of the ", pop, " subclusters (", levels_contrast[2], " vs ", levels_contrast[1], ")"),
                  subtitle = paste0("each point is a neighbourhood at the position of its index cell"))

ggsave(paste0(folder_plot, "08d_miloR_", pop, "_nhoodLogFC_subcluster.pdf"), plot = p_nhood,
       width = 5 * ncol_panel + 1, height = 4.5 * ceiling(length(sub_levels) / ncol_panel) + 0.5, limitsize = FALSE)

# snippet: same highlight plot (panel 2) from a Seurat object -------------
# NOT run here (the script works on the milo object). Use it if the milo object is not available or on any other Seurat object.
# the final object is 03_sobj_integrated_cleanup_manualAnnotation_subclusterFiltered_v2.rds (7 GB, load it only if needed)
#
library(Seurat)
sobj <- readRDS("../../out/object/analysis_R45_pixi/03_sobj_integrated_cleanup_manualAnnotation_subclusterFiltered_v2.rds")

# subset to the population
test <- subset(sobj, cell_id == "ASTRO")

# barcodes of each subcluster, as a named list
split_cells <- split(colnames(test), test$cell_id_subcluster3)

# UMAP coordinates of the population, to place the name of the subset on the selected cells (median of the coordinates)
emb <- as.data.frame(Embeddings(test, "umap"))
colnames(emb) <- c("UMAP1", "UMAP2")

# one panel per subcluster: highlighted cells in red on the grey background of the whole population
list_plot <- imap(split_cells, function(x, nm){
  # position of the label: median of the selected cells
  df_lab <- data.frame(UMAP1 = median(emb[x, "UMAP1"]),
                       UMAP2 = median(emb[x, "UMAP2"]),
                       label = nm)

  DimPlot(test,
          reduction = "umap",
          cells.highlight = x,
          cols.highlight = "red",
          cols = "grey85",
          sizes.highlight = 0.1,
          pt.size = 0.3) +
    # name of the subset on the selected cells
    # geom_label(data = df_lab, aes(x = UMAP1, y = UMAP2, label = label), size = 3, fill = alpha("white", 0.7), label.size = 0, inherit.aes = FALSE) +
    ggtitle(paste0(nm, " (n = ", length(x), ")")) +
    NoLegend()
})

wrap_plots(list_plot, ncol = 3)

# NOTICE: split.by in DimPlot (as in the first line of the original snippet) shows only the cells of the subcluster in each panel, without the grey background of the rest of the population; cells.highlight keeps the background, which is easier to read

message("\nAll done.")
