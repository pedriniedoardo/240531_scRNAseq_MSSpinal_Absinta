# AIM ---------------------------------------------------------------------
# run propeller on the final consensus object generated in 03_apply_filter_annotate_fullObject.R.
# test the cell type proportion differences for diagnosis (MS vs CTRL), pathological stage and the updated demyelination classes (WM and GM).
# two levels of testing:
# - macro level: cell_id proportions over the whole object
# - subcluster level: cell_id_subcluster3 proportions computed within each macro population, defined as the prefix of the subcluster label (macro_pop: LYM, ASTRO, IMMUNE, ...; not cell_id, which is assigned before the subclustering) (share of that population per sample), so that the test is not driven by the abundance of the macro population itself
# readapted from R_analysis_R44/scr/26_propeller.R and R_analysis_R44/scr/002_update propeller_demyelination.R

# LIBRARIES ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(speckle)
library(limma)
library(statmod)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(patchwork)

# PARAMETERS ------------------------------------------------------------------
options(stringsAsFactors = FALSE)

# contrasts to test. the first level is the reference
list_contrast <- list(
  diagnosis = list(var = "diagnosis_short", levels = c("CTRL","MS"), ref = "CTRL"),
  pathology = list(var = "pathological_stage_short", levels = c("CTRL","IN","ACT"), ref = "CTRL"),
  demye_WM = list(var = "demye_WM_class", levels = c("no","low","high"), ref = "no"),
  demye_GM = list(var = "demye_GM_class", levels = c("no","low","high"), ref = "no"),
  demye_tot = list(var = "demye_tot_class", levels = c("no","low","high"), ref = "no"),
  location = list(var = "location", levels = c("cervical","thoracic","lumbar"), ref = "cervical")
)

# for the subcluster level, drop samples with too few cells of the given macro population, as their proportions are too noisy
min_cells_sample <- 10

# in the plots, flag (red point + sample label) the samples with fewer cells than this in the tested population
min_cells_flag <- 50

folder_table <- "../../out/table/analysis_R45_pixi/04_propeller/"
folder_plot <- "../../out/plot/analysis_R45_pixi/04_propeller/"
# dir.create(folder_table, recursive = TRUE, showWarnings = FALSE)
# dir.create(folder_plot, recursive = TRUE, showWarnings = FALSE)

# read in the data --------------------------------------------------------
sobj <- readRDS("../../out/object/analysis_R45_pixi/03_sobj_integrated_cleanup_manualAnnotation_subclusterFiltered_v2.rds")
sobj@meta.data

# confirm the dataset
DimPlot(sobj,group.by = "cell_id",label = T,repel = T)

# only the metadata is needed
meta_ref <- sobj@meta.data %>%
  mutate(across(all_of(map_chr(list_contrast, "var")), as.character)) %>%
  mutate(location = factor(location, levels = c("cervical","thoracic","lumbar"))) %>%
  # cell_id is assigned before the subclustering, cell_id_subcluster3 after it, so the two are not nested (e.g. LYM|IgG-secreting plasmablasts are in both B CELLS and T CELLS).
  # the macro population of a subcluster is the prefix of its label (LYM, ASTRO, IMMUNE, ...)
  mutate(macro_pop = str_remove(cell_id_subcluster3, "\\|.*$"))

rm(sobj)

# confirm that each subcluster maps to a single macro_pop
table(meta_ref$cell_id_subcluster3, meta_ref$macro_pop)
gc()

# confirm the covariates per sample
meta_ref %>%
  group_by(sample_id, location, diagnosis_short, pathological_stage_short, demye_WM_class, demye_GM_class, demye_tot_class) %>%
  summarise(n = n(), .groups = "drop") %>%
  print(n = Inf)

# function for a single test ----------------------------------------------
# run propeller for one covariate on one set of cells, save the tables and the plots
# - meta: cells to test (already restricted to the population of interest)
# - cluster_col: column with the labels whose proportions are tested
# - tag: used in the file names and in the plot titles
# - contrast_name: name of the covariate in list_contrast

# meta <- meta_ref
# cluster_col <- "cell_id"
# tag <- "cellid"
# contrast_name <- "diagnosis"
run_propeller <- function(meta, cluster_col, tag, contrast_name){
  ct <- list_contrast[[contrast_name]]

  meta_test <- meta %>%
    dplyr::select(ct$var,all_of(cluster_col),sample_id) %>%
    filter(!is.na(.data[[ct$var]])) %>%
    mutate(group = factor(.data[[ct$var]], levels = ct$levels) %>% droplevels(),
           cluster = .data[[cluster_col]])

  if(nlevels(meta_test$group) < 2 | n_distinct(meta_test$cluster) < 2){
    message("skip ", tag, " ", contrast_name, ": less than 2 groups or clusters")
    return(NULL)
  }

  # Run propeller testing for cell type proportion differences between the groups.
  # cluster is the cluster/celltype id
  # sample is the reference id of the biological replicates
  # group is the grouping id
  out <- propeller(clusters = meta_test$cluster,
                   sample = meta_test$sample_id,
                   group = meta_test$group) %>%
    rownames_to_column("cluster") %>%
    mutate(tag = tag,
           contrast = contrast_name,
           var = ct$var)

  out %>%
    write_tsv(paste0(folder_table, "04_propeller_", tag, "_", contrast_name, ".tsv"))

  # proportions per sample, include the missing clusters as 0
  df_summary <- meta_test %>%
    count(sample_id, group, cluster, name = "n") %>%
    complete(nesting(sample_id, group), cluster, fill = list(n = 0)) %>%
    group_by(sample_id) %>%
    mutate(tot = sum(n),
           prop = n/tot,
           # tot = cells of the sample in the tested population
           low_cells = tot < min_cells_flag) %>%
    ungroup()

  df_summary %>%
    write_tsv(paste0(folder_table, "04_propeller_prop_", tag, "_", contrast_name, ".tsv"))

  # scale the plot size to the number of clusters
  n_cluster <- n_distinct(df_summary$cluster)
  n_col <- ceiling(sqrt(n_cluster))
  n_row <- ceiling(n_cluster / n_col)

  # plot 01
  p01 <- df_summary %>%
    ggplot(aes(x = group, y = prop)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(aes(color = low_cells), position = position_jitter(width = 0.1, height = 0, seed = 1), shape = 1, alpha = 0.7) +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"), guide = "none") +
    facet_wrap(~cluster, scales = "free", ncol = n_col) +
    theme_bw() +
    theme(strip.background = element_blank(), axis.text.x = element_text(hjust = 1, angle = 45)) +
    labs(x = ct$var, title = paste(tag, contrast_name))
  ggsave(paste0(folder_plot, "04_propeller_plot01_", tag, "_", contrast_name, ".pdf"), plot = p01, width = 2.5 * n_col + 1, height = 2.5 * n_row + 1, limitsize = FALSE)

  # plot 02
  p02 <- df_summary %>%
    ggplot() +
    geom_boxplot(aes(x = cluster, y = prop, color = group), outlier.shape = NA) +
    geom_point(aes(x = cluster, y = prop, color = group), position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0, dodge.width = 0.8), alpha = 0.7) +
    theme_cowplot() +
    theme(axis.text.x = element_text(hjust = 1, angle = 90)) +
    scale_y_sqrt() +
    labs(color = ct$var, title = paste(tag, contrast_name))
  ggsave(paste0(folder_plot, "04_propeller_plot02_", tag, "_", contrast_name, ".pdf"), plot = p02, width = 0.5 * n_cluster + 4, height = 6, limitsize = FALSE)

  out
}

# macro level: cell_id over the whole object ------------------------------
# one test per covariate, all the cells are used
list_res_macro <- map(names(list_contrast), function(x){
  run_propeller(meta = meta_ref,
                cluster_col = "cell_id",
                tag = "cellid",
                contrast_name = x) 
})

# global level: cell_id_subcluster3 over the whole object -----------------
# proportions are relative to all the cells of the sample, so a subcluster can change only because its parent population changes
list_res_global <- map(names(list_contrast), function(x){
  run_propeller(meta = meta_ref,
                cluster_col = "cell_id_subcluster3",
                tag = "subcluster3_all",
                contrast_name = x)
})

# subcluster level: cell_id_subcluster3 within each macro population ------
# skip the populations that have a single label (e.g. not subclustered)
pop_test <- meta_ref %>%
  group_by(macro_pop) %>%
  summarise(n_sub = n_distinct(cell_id_subcluster3)) %>%
  filter(n_sub > 1) %>%
  pull(macro_pop)

pop_test

# cells of one population, dropping the samples with too few of them (also used for the heatmaps)
get_meta_pop <- function(pop){
  meta_ref %>%
    dplyr::select(cell_id_subcluster3, macro_pop, sample_id, all_of(map_chr(list_contrast, "var") %>% unname())) %>%
    filter(macro_pop == pop) %>%
    group_by(sample_id) %>%
    filter(n() >= min_cells_sample) %>%
    ungroup()
}

# tag made safe for the file names
make_tag <- function(pop){
  paste0("subcluster3_", str_replace_all(pop, "[^A-Za-z0-9]+", "_") %>% str_remove("_$"))
}

# for each population: keep its cells, drop the samples with too few of them, then one test per covariate
# pop <- "VAS"
list_res_sub <- map(pop_test, function(pop){
  meta_pop <- get_meta_pop(pop)
  tag <- make_tag(pop)

  map(names(list_contrast), function(x) {
    run_propeller(meta = meta_pop,
                  cluster_col = "cell_id_subcluster3",
                  tag = tag,
                  contrast_name = x)
    })
}) %>%
  unlist(recursive = FALSE)

# summary of all the tests ------------------------------------------------
# level: cell_id (whole object), subcluster3_global (whole object), subcluster3_withinPop (share within the macro population)
# FDR is the one from propeller (within each test). FDR_level is recomputed over all the tests of the same level and contrast,
# so for subcluster3_withinPop it is corrected across all the populations together
df_summary_all <- bind_rows(list_res_macro, list_res_global, list_res_sub) %>%
  mutate(level = case_when(tag == "cellid" ~ "cell_id",
                           tag == "subcluster3_all" ~ "subcluster3_global",
                           TRUE ~ "subcluster3_withinPop")) %>%
  group_by(level, contrast) %>%
  mutate(FDR_level = p.adjust(P.Value, method = "BH")) %>%
  ungroup() %>%
  relocate(level, tag, contrast, var, cluster)

df_summary_all %>%
  write_tsv(paste0(folder_table, "04_propeller_summary_allLevels.tsv"))

df_summary_all %>%
  filter(FDR < 0.05 | FDR_level < 0.05) %>%
  select(level, tag, contrast, cluster, P.Value, FDR, FDR_level) %>%
  print(n = Inf)

# heatmap of the proportions per sample -----------------------------------
# sample level metadata for the annotation
meta_sample <- meta_ref %>%
  distinct(sample_id, diagnosis_short, pathological_stage_short, location)

sample_prop <- meta_ref %>%
  count(sample_id, cell_id, name = "n") %>%
  complete(sample_id, cell_id, fill = list(n = 0)) %>%
  group_by(sample_id) %>%
  mutate(tot = sum(n),
         prop = n/tot) %>%
  ungroup() %>%
  group_by(cell_id) %>%
  # logit transform as in propeller (clamped to avoid infinite values for the zeros), then z-score per cell type
  mutate(logit_prop = qlogis(pmin(pmax(prop, 1e-4), 1 - 1e-4)),
         scale_prop = scale(prop) %>% .[,1],
         scale_logit_prop = scale(logit_prop) %>% .[,1]) %>%
  ungroup()

sample_prop_wide <- sample_prop %>%
  select(cell_id, sample_id, prop) %>%
  pivot_wider(names_from = sample_id, values_from = prop) %>%
  column_to_rownames("cell_id") %>%
  as.matrix()

sample_prop_wide2 <- sample_prop %>%
  select(cell_id, sample_id, scale_prop) %>%
  pivot_wider(names_from = sample_id, values_from = scale_prop) %>%
  column_to_rownames("cell_id") %>%
  as.matrix()

sample_prop_wide3 <- sample_prop %>%
  select(cell_id, sample_id, scale_logit_prop) %>%
  pivot_wider(names_from = sample_id, values_from = scale_logit_prop) %>%
  column_to_rownames("cell_id") %>%
  as.matrix()

# confirm the dimensions and the scaling
colSums(sample_prop_wide)
apply(sample_prop_wide2, 1, mean)
apply(sample_prop_wide2, 1, sd)
apply(sample_prop_wide3, 1, mean)
apply(sample_prop_wide3, 1, sd)

meta_sample_prop <- data.frame(sample_id = colnames(sample_prop_wide)) %>%
  left_join(meta_sample, by = "sample_id")

column_meta_sample_prop <- HeatmapAnnotation(diagnosis = meta_sample_prop$diagnosis_short,
                                             # pathology = meta_sample_prop$pathological_stage_short,
                                             # location = as.character(meta_sample_prop$location),
                                             col = list(
                                               diagnosis = c("CTRL" = "blue",
                                                             "MS" = "red")
                                               # pathology = c("CTRL" = "blue",
                                               #               "IN" = "orange",
                                               #               "ACT" = "red")
                                               # location = c("cervical" = "gray",
                                               #              "thoracic" = "gray30",
                                               #              "lumbar" = "black")
                                               ))

# ht1 <- Heatmap(sample_prop_wide, show_column_names = TRUE, show_row_dend = FALSE,
#                name = "prop_cell_type",
#                column_title = "sample",
#                col = viridis::viridis(option = "turbo", n = 10),
#                top_annotation = column_meta_sample_prop, show_row_names = TRUE)
# 
# pdf(paste0(folder_plot, "04_heatmap_cellType_prop.pdf"), width = 12, height = 6)
# draw(ht1, heatmap_legend_side = "left", annotation_legend_side = "left", padding = unit(c(30, 2, 2, 2), "mm"))
# dev.off()

# ht2 <- Heatmap(sample_prop_wide2[c("B CELLS"),,drop = F], show_column_names = TRUE, show_row_dend = FALSE,
#                name = "scaled_prop",
#                column_title = "sample",
#                col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
#                top_annotation = column_meta_sample_prop, show_row_names = TRUE)

ht2 <- Heatmap(sample_prop_wide2, show_column_names = TRUE, show_row_dend = FALSE,
               name = "scaled_prop",
               column_title = "sample",
               col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
               top_annotation = column_meta_sample_prop, show_row_names = TRUE)

# ht2 <- Heatmap(sample_prop_wide[c("B CELLS"),,drop = F], show_column_names = TRUE, show_row_dend = FALSE,
#                name = "scaled_prop",
#                column_title = "sample",
#                # col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
#                top_annotation = column_meta_sample_prop, show_row_names = TRUE)

pdf(paste0(folder_plot, "04_heatmap_cellType_prop_scaled.pdf"), width = 8, height = 6)
draw(ht2, heatmap_legend_side = "left", annotation_legend_side = "left", padding = unit(c(30, 2, 2, 2), "mm"))
dev.off()

# same as ht2 but on the logit scale (the one tested by propeller), so a few outlier samples do not dominate the color scale
ht3 <- Heatmap(sample_prop_wide3, show_column_names = TRUE, show_row_dend = FALSE,
               name = "scaled_logit_prop",
               column_title = "sample",
               col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
               top_annotation = column_meta_sample_prop, show_row_names = TRUE)

pdf(paste0(folder_plot, "04_heatmap_cellType_prop_scaled_logit.pdf"), width = 8, height = 6)
draw(ht3, heatmap_legend_side = "left", annotation_legend_side = "left", padding = unit(c(30, 2, 2, 2), "mm"))
dev.off()

# only for B cells
ht32 <- Heatmap(sample_prop_wide3["B CELLS",,drop = F], show_column_names = TRUE, show_row_dend = FALSE,
                name = "scaled_logit_prop",
                column_title = "sample",
                col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                cluster_columns = T,,
                top_annotation = column_meta_sample_prop, show_row_names = TRUE)

pdf(paste0(folder_plot, "04_heatmap_cellType_prop_scaled_logit_02.pdf"), width = 8, height = 5)
draw(ht32, heatmap_legend_side = "left", annotation_legend_side = "left", padding = unit(c(30, 2, 2, 2), "mm"))
dev.off()

# same, with the samples split by diagnosis and not clustered
ht4 <- Heatmap(sample_prop_wide3, show_column_names = TRUE, show_row_dend = FALSE,
               name = "scaled_logit_prop",
               column_title = "sample",
               col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
               cluster_columns = FALSE,
               column_split = factor(meta_sample_prop$diagnosis_short, levels = c("CTRL","MS")),
               top_annotation = column_meta_sample_prop, show_row_names = TRUE)

pdf(paste0(folder_plot, "04_heatmap_cellType_prop_scaled_logit_splitDiagnosis.pdf"), width = 8, height = 6)
draw(ht4, heatmap_legend_side = "left", annotation_legend_side = "left", padding = unit(c(30, 2, 2, 2), "mm"))
dev.off()

# boxplots of the significant comparisons ---------------------------------
# one pdf for each level/tag and contrast with at least one significant cluster (propeller FDR or FDR_level < fdr_boxplot)
# the proportions per sample are the ones saved by run_propeller, so they come from the same cells as the test
# strip label reports the propeller FDR and the FDR over the level
fdr_boxplot <- 0.05

df_sig <- df_summary_all %>%
  filter(FDR < fdr_boxplot | FDR_level < fdr_boxplot) %>%
  arrange(P.Value)

# tag <- "subcluster3_all"
# contrast <- "diagnosis"
list_plot <- pmap(distinct(.data = df_sig, tag, contrast), function(tag, contrast){
  tg <- tag
  cn <- contrast
  ct <- list_contrast[[cn]]

  df_lab <- df_sig %>%
    filter(tag == tg, contrast == cn) %>%
    mutate(label = paste0(cluster, "\nFDR=", signif(FDR, 2), " FDR_level=", signif(FDR_level, 2)))

  df_plot <- read_tsv(paste0(folder_table, "04_propeller_prop_", tg, "_", cn, ".tsv"), show_col_types = FALSE) %>%
    inner_join(df_lab %>% select(cluster, label), by = "cluster") %>%
    mutate(group = factor(group, levels = ct$levels),
           label = factor(label, levels = df_lab$label))

  n_cluster <- nrow(df_lab)
  n_col <- min(n_cluster, 4)
  n_row <- ceiling(n_cluster / n_col)

  p <- df_plot %>%
    ggplot(aes(x = group, y = prop)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(aes(color = low_cells), position = position_jitter(width = 0.1, height = 0, seed = 1), shape = 1, alpha = 0.7) +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"), guide = "none") +
    facet_wrap(~label, scales = "free", ncol = n_col) +
    theme_bw() +
    theme(strip.background = element_blank(), axis.text.x = element_text(hjust = 1, angle = 45)) +
    labs(x = ct$var, title = paste(tg, cn))
  # ggsave(paste0(folder_plot, "04_boxplot_sig_", tg, "_", cn, ".pdf"), plot = p, width = 3 * n_col + 1, height = 3 * n_row + 1, limitsize = FALSE)
  return(p)
})

p01 <- wrap_plots(list_plot)
ggsave(plot = p01,filename = paste0(folder_plot, "04_boxplot_significant_prop.pdf"), width = 15, height = 12)

# barplot of the cells per sample in each macro population ----------------
# one pdf per population in pop_test: cells of each sample, stacked by cell_id_subcluster3
# sample labels in red = removed from the subcluster tests (fewer than min_cells_sample cells), dashed line = min_cells_sample
# samples with no cells of the population are shown with an empty bar (they are removed as well)
all_sample <- distinct(meta_ref, sample_id)

# pop <- "B CELLS"
walk(pop_test, function(pop){
  df_count <- meta_ref %>%
    filter(macro_pop == pop) %>%
    count(sample_id, cell_id_subcluster3, name = "n")

  df_tot <- df_count %>%
    group_by(sample_id) %>%
    summarise(tot = sum(n), .groups = "drop") %>%
    right_join(all_sample, by = "sample_id") %>%
    mutate(tot = replace_na(tot, 0),
           filtered = tot < min_cells_sample) %>%
    arrange(tot)

  df_plot <- df_count %>%
    mutate(sample_id = factor(sample_id, levels = df_tot$sample_id))

  p <- ggplot() +
    geom_col(data = df_plot, aes(x = sample_id, y = n, fill = cell_id_subcluster3)) +
    # keep the samples without cells on the axis
    geom_blank(data = df_tot %>% mutate(sample_id = factor(sample_id, levels = df_tot$sample_id)), aes(x = sample_id, y = tot)) +
    geom_hline(yintercept = min_cells_sample, linetype = "dashed") +
    scale_x_discrete(drop = FALSE) +
    theme_cowplot() +
    theme(axis.text.x = element_text(hjust = 1, angle = 90, color = ifelse(df_tot$filtered, "red", "black"))) +
    labs(x = "sample", y = "n cells", fill = "cell_id_subcluster3",
         title = paste0(pop, ": cells per sample"),
         subtitle = paste0("red = removed (< ", min_cells_sample, " cells), ", sum(df_tot$filtered), "/", nrow(df_tot), " samples"))

  ggsave(paste0(folder_plot, "04_barplot_ncells_", make_tag(pop), ".pdf"), plot = p, width = 0.3 * nrow(df_tot) + 5, height = 6, limitsize = FALSE)
})
