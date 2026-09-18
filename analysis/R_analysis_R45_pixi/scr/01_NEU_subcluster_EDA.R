# AIM ---------------------------------------------------------------------
# script to run some specific EDA over the NEU subclsuter:
# proportion analysis
# produce a table of markers per cluster at all the resolutions

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

# define the inputs ----------------------------------------------------------
#
sobj_input_id <- "../../out/object/analysis_R44/27_NEU_subcluster_HarmonySample.rds"
message("input single cell object: ", sobj_input_id)

# this is already a subset of a specific cell type

# column with the resolution of interest (used for the propeller test)
col_subcluster <- "RNA_snn_res.0.1"

# column for the grouping condition
col_group <- "diagnosis_short"

# column for the sample id
col_donor <- "sample_id"

# output locations, following the project convention (out/<kind>/analysis_R45_pixi/)
dir_out_plot <- "../../out/plot/analysis_R45_pixi/"
dir_out_tab <- "../../out/table/analysis_R45_pixi/"
dir.create(dir_out_plot, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_out_tab, recursive = TRUE, showWarnings = FALSE)

# read in the object ----------------------------------------------------------
scobj <- readRDS(sobj_input_id)
# confirm the identity of the dataset
DimPlot(scobj, label = T, raster = F,repel = T, group.by = col_subcluster)

# wrangling ---------------------------------------------------------------
# work with fixed column names (group/donor/subcluster) for the rest of the script, so that only the three col_* variables above need to change
meta_test <- scobj@meta.data %>%
  transmute(group = .data[[col_group]],
            donor = .data[[col_donor]],
            subcluster = .data[[col_subcluster]]) %>%
  # it would be fine even just donor in this case as it is unique anyway. this implementation allows to handle cases where there is the same sample name across groups, and keep them as separated entities
  mutate(sample_prop = paste0(donor, group))

# run the test ------------------------------------------------------------
table(meta_test$donor, meta_test$group)

out_neu_subcluster <- propeller(
  clusters = meta_test$subcluster,
  sample = meta_test$sample_prop,
  group = meta_test$group
)

out_neu_subcluster %>%
  rownames_to_column("subcluster") %>%
  write_tsv(file.path(dir_out_tab, "01_propeller_NEU_subcluster.tsv"))

# plotting ------------------------------------------------------------------
# named palette, so each subcluster keeps its colour regardless of the legend order
pal_sub <- setNames(pals::glasbey(n_distinct(meta_test$subcluster)), levels(factor(meta_test$subcluster)))
scales::show_col(pal_sub)

# default plot
speckle::plotCellTypeProps(
  x = scobj,
  clusters = meta_test$subcluster,
  sample = meta_test$group
) + scale_fill_manual(values = pal_sub) + theme_minimal() + theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1, angle = 45))
ggsave(file.path(dir_out_plot, "01_plot_propeller_NEU_subcluster.pdf"), height = 5, width = 5)

# custom plot
# complete() adds the missing sample x subcluster combinations with n = 0, otherwise the zeros are dropped from the boxplots
df_summary_neu <- meta_test %>%
  count(sample_prop, donor, group, subcluster, name = "n") %>%
  complete(nesting(sample_prop, donor, group), subcluster, fill = list(n = 0)) %>%
  group_by(sample_prop) %>%
  mutate(tot = sum(n),
         prop = n / tot) %>%
  ungroup()

df_summary_neu %>%
  write_tsv(file.path(dir_out_tab, "01_df_summary_NEU_subcluster.tsv"))

# plot 01: one panel per NEU subcluster, boxplot + jittered points (one point per sample)
df_summary_neu %>%
  ggplot(aes(x = group, y = prop)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1), shape = 1, alpha = 0.7) +
  facet_wrap(~subcluster, scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank(), axis.text.x = element_text(hjust = 1, angle = 45))
ggsave(file.path(dir_out_plot, "01_plot_propeller_NEU_subcluster_replicates.pdf"), width = 10, height = 10)

# plot 02: all NEU subclusters side by side, colored/dodged by group
df_summary_neu %>%
  ggplot() +
  geom_boxplot(aes(x = subcluster, y = prop, color = group), outlier.shape = NA) +
  geom_point(aes(x = subcluster, y = prop, color = group), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), alpha = 0.7) +
  theme_cowplot() +
  theme(axis.text.x = element_text(hjust = 1, angle = 90)) +
  scale_y_sqrt()
ggsave(file.path(dir_out_plot, "01_plot_propeller_NEU_subcluster_byGroup.pdf"), width = 8, height = 5)

# markers per cluster at all the resolutions --------------------------------
# Seurat v5 runs the wilcoxon test through presto automatically when presto is installed, so FindAllMarkers is already the fast implementation (SeuratWrappers::RunPrestoAll is not needed)
DefaultAssay(scobj) <- "RNA"

# all the clustering resolutions stored in the object
id_resolution <- str_subset(colnames(scobj@meta.data), pattern = "^RNA_snn_res\\.") %>%
  sort()
message("resolutions found: ", paste(id_resolution, collapse = ", "))

# plot panels dim reductions
list_plot <- map(id_resolution,function(x){
  plot <- DimPlot(scobj,reduction = "umap",group.by = x,label = T,repel = T)+ NoLegend() 
  return(plot)
})

plot_panel <-  wrap_plots(list_plot)
ggsave(plot = plot_panel,filename = file.path(dir_out_plot, "01_plot_panel_NEU_subcluster.pdf"), width = 25, height = 15)

# one FindAllMarkers per resolution, the result is a list named by resolution the Idents assignment is local to the function, so scobj is not modified
list_markers <- id_resolution %>%
  set_names() %>%
  map(function(id_res) {
    message("FindAllMarkers at ", id_res)
    Idents(scobj) <- id_res

    # find markers for every cluster compared to all remaining cells, report only the positive ones
    FindAllMarkers(scobj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
      arrange(cluster, p_val_adj, desc(avg_log2FC))
  })

# save the tables, all the resolutions stacked in one table (the cluster ids overlap across resolutions, the resolution column tells them apart)
df_markers_all <- list_markers %>%
  bind_rows(.id = "resolution")

df_markers_all %>%
  write_tsv(file.path(dir_out_tab, "01_FindAllMarkers_allResolutions_NEU_subcluster.tsv"))

# top 100 per cluster, within each resolution
df_markers_all %>%
  group_by(resolution, cluster) %>%
  slice_head(n = 100) %>%
  ungroup() %>%
  write_tsv(file.path(dir_out_tab, "01_FindAllMarkers_allResolutions_NEU_subcluster_top100.tsv"))

# top 100 per cluster, within each resolution, after removing MT, ribosomal and globin genes
df_markers_all %>%
  filter(str_detect(gene, pattern = "^MT-", negate = T),
         str_detect(gene, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", negate = T),
         str_detect(gene, pattern = "^HB[^(P)]", negate = T)) %>%
  group_by(resolution, cluster) %>%
  slice_head(n = 100) %>%
  ungroup() %>%
  write_tsv(file.path(dir_out_tab, "01_FindAllMarkers_allResolutions_NEU_subcluster_top100_noRIBOandMT.tsv"))
