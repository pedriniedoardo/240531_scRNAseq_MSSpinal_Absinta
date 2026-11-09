# AIM ---------------------------------------------------------------------
# run the single cell DE analysis (FindMarkers, accelerated by presto) on the final consensus object generated in 03_apply_filter_annotate_fullObject.R.
# comparisons: diagnosis (MS vs CTRL) and pathological stage (IN/ACT vs CTRL).
# run on all the levels of the macro annotation (cell_id) and of the consensus subcluster annotation (cell_id_subcluster3).
# for each comparison and group the script also:
# 1. summarises the number of DEGs (padj < padj_thr and |avg_log2FC| > lfc_thr, split in up and down) and computes the "perturbation" score = number of DEGs / number of cells in the comparison. barplot of the score per group and UMAP colored by the score.
#    NOTICE: the number of DEGs depends on the statistical power (number of cells), dividing by the number of cells is a rough normalization, interpret the score as a relative ranking between groups.
# scDist on the same comparisons, and its comparison with the perturbation score, is run in a separate job (05b_scDist_fullObject.R) because of the memory usage. Run this script first: 05b reads 05_DE_sc_summary.tsv.
# readapted from R_analysis_R44/scr/30_DE_analysis_sc.R and from scRNA_SOD_3rd/analysis/R_analysis_R45_pixi/scr/general/24_DEG_SOD_vs_CT_perturbation_scDist_annotation.R
# NOTICE: the DE is run on the RNA assay (log normalized data layer).

# LIBRARIES ---------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(SeuratWrappers)
library(limma)
library(ggrepel)
library(presto)
library(patchwork)
library(cowplot)
library(scDist)

# PARAMETERS ------------------------------------------------------------------
options(stringsAsFactors = FALSE)
options(Seurat.object.assay.version = "v5")

# contrasts to test. the first level is the reference
list_contrast <- list(
  diagnosis = list(var = "diagnosis_short", levels = c("CTRL","MS"), ref = "CTRL"),
  pathology = list(var = "pathological_stage_short", levels = c("CTRL","IN","ACT"), ref = "CTRL")
)

# annotation levels to test
list_level <- c("cell_id", "cell_id_subcluster3")

# minimum number of cells per group to run the comparison
min_cells <- 10
# thresholds used to call a gene as DEG in the summary and in the volcano (the full tables are saved without any threshold)
padj_thr <- 0.05
lfc_thr <- 1

folder_table <- "../../out/table/analysis_R45_pixi/05_DE_sc/"
folder_plot <- "../../out/plot/analysis_R45_pixi/05_DE_sc/"

# read in the data --------------------------------------------------------
sobj <- readRDS("../../out/object/analysis_R45_pixi/03_sobj_integrated_cleanup_manualAnnotation_subclusterFiltered_v2.rds")

# make sure the correct default assay is loaded, should be RNA
DefaultAssay(sobj) <- "RNA"
class(sobj@assays$RNA)

# check that in the RNA slot the data object is indeed loaded with normalized values
LayerData(sobj, assay = "RNA", layer = "data")[1:10, 1:10]

# metadata and UMAP coordinates used for the score plots
meta_all <- sobj@meta.data
df_umap <- Embeddings(sobj, "umap") %>%
  as.data.frame() %>%
  set_names(c("UMAP_1", "UMAP_2")) %>%
  rownames_to_column("barcode")

# confirm the groups
lapply(list_contrast, function(ct){
  table(sobj@meta.data[[ct$var]], sobj$cell_id, useNA = "ifany")
})

# run DE ------------------------------------------------------------------
# avg_log2FC: log fold-change of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the first group (the non reference level)
df_run <- crossing(level = list_level, contrast = names(list_contrast))


# level <- "cell_id"
# contrast <- "diagnosis"
list_out <- pmap(df_run, function(level, contrast){
  ct <- list_contrast[[contrast]]

  # define the grouping for the comparison
  sobj$condition.annotation <- paste(sobj@meta.data[[ct$var]], sobj@meta.data[[level]], sep = "_")
  Idents(sobj) <- "condition.annotation"
  n_ident <- table(sobj$condition.annotation)

  clusters_id <- as.character(sort(unique(sobj@meta.data[[level]])))
  lvl_test <- setdiff(ct$levels, ct$ref)
  
  # lvl <- "MS"
  lapply(lvl_test, function(lvl){
    # x <- "ASTRO"
    df_res <- lapply(clusters_id, function(x){
      id_1 <- paste0(lvl, "_", x)
      id_2 <- paste0(ct$ref, "_", x)
      n_1 <- ifelse(id_1 %in% names(n_ident), n_ident[[id_1]], 0)
      n_2 <- ifelse(id_2 %in% names(n_ident), n_ident[[id_2]], 0)

      if(n_1 < min_cells | n_2 < min_cells){
        message("skip ", id_1, " vs ", id_2, ": ", n_1, " vs ", n_2, " cells")
        return(NULL)
      }

      message("run ", id_1, " vs ", id_2)
      tryCatch({
        # FindMarkers uses presto automatically when the package is installed (Seurat v5)
        FindMarkers(sobj, assay = "RNA", ident.1 = id_1, ident.2 = id_2, min.pct = 0.05, logfc.threshold = 0, verbose = FALSE) %>%
          rownames_to_column("gene") %>%
          mutate(id_1 = id_1,
                 id_2 = id_2,
                 n_1 = n_1,
                 n_2 = n_2,
                 cluster = x)
      }, error = function(e){
        message("failed ", id_1, " vs ", id_2, ": ", conditionMessage(e))
        NULL
      })
    }) %>%
      bind_rows()

    if(nrow(df_res) == 0){
      return(NULL)
    }

    df_res <- df_res %>%
      mutate(level = level,
             contrast = contrast,
             comparison = paste0(lvl, "_vs_", ct$ref))

    df_res %>%
      write_tsv(paste0(folder_table, "05_DE_sc_", level, "_", contrast, "_", lvl, "_vs_", ct$ref, ".tsv"))

    return(df_res)
  }) %>%
    bind_rows()
})

df_res_all <- bind_rows(list_out)

# plot the results --------------------------------------------------------
# one pvalue distribution and one volcano per level/contrast/comparison
df_res_all %>%
  split(f = list(.$level, .$contrast, .$comparison), drop = TRUE) %>%
  walk(function(df){
    n_cluster <- n_distinct(df$cluster)
    n_col <- ceiling(sqrt(n_cluster))
    n_row <- ceiling(n_cluster / n_col)
    file_id <- paste(unique(df$level), unique(df$contrast), unique(df$comparison), sep = "_")

    # show the distribution of the pvalues
    p01 <- df %>%
      ggplot(aes(x = p_val)) +
      geom_histogram(breaks = 0:20/20) +
      theme_bw() +
      facet_wrap(~cluster, ncol = n_col) +
      theme(strip.background = element_blank()) +
      ggtitle(file_id)
    ggsave(paste0(folder_plot, "05_dist_p_value_", file_id, ".pdf"), plot = p01, width = 3 * n_col + 1, height = 2.5 * n_row + 1, limitsize = FALSE)

    # render all of them as a volcano plot
    test_significant <- df %>%
      filter(abs(avg_log2FC) > lfc_thr & p_val_adj < padj_thr)

    # label only the top genes per cluster to keep the plot readable
    test_label <- test_significant %>%
      group_by(cluster) %>%
      arrange(p_val_adj) %>%
      dplyr::slice(1:10) %>%
      ungroup()

    p02 <- df %>%
      ggplot(aes(x = avg_log2FC, y = -log(p_val_adj))) +
      geom_point(alpha = 0.01) +
      geom_point(data = test_significant, col = "red", alpha = 0.5) +
      geom_vline(xintercept = c(-lfc_thr, lfc_thr), col = "gray", linetype = "dashed") +
      geom_hline(yintercept = (-log(padj_thr)), col = "gray", linetype = "dashed") +
      geom_text_repel(data = test_label, aes(label = gene), size = 2, segment.alpha = 0.4) +
      facet_wrap(~cluster, ncol = n_col) +
      theme_bw() +
      theme(strip.background = element_blank()) +
      ggtitle(file_id)
    ggsave(paste0(folder_plot, "05_volcano_", file_id, ".pdf"), plot = p02, width = 4 * n_col + 1, height = 3.5 * n_row + 1, limitsize = FALSE)
  })

# summary of the number of DEGs and perturbation score --------------------
# only the tested groups are in df_res_all. n_1 + n_2 = number of cells in the comparison
df_summary <- df_res_all %>%
  group_by(level, contrast, comparison, cluster, n_1, n_2) %>%
  summarise(n_up = sum(p_val_adj < padj_thr & avg_log2FC > lfc_thr, na.rm = TRUE),
            n_down = sum(p_val_adj < padj_thr & avg_log2FC < -lfc_thr, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(n_cells = n_1 + n_2,
         n_deg = n_up + n_down,
         score = n_deg / n_cells,
         score_up = n_up / n_cells,
         score_down = n_down / n_cells)

df_summary %>%
  write_tsv(paste0(folder_table, "05_DE_sc_summary.tsv"))

df_summary %>%
  print(n = Inf)

# one barplot of the score (stacked by direction) and one UMAP colored by the score per level/contrast/comparison
df_summary %>%
  split(f = list(.$level, .$contrast, .$comparison), drop = TRUE) %>%
  walk(function(df){
    level <- unique(df$level)
    file_id <- paste(level, unique(df$contrast), unique(df$comparison), sep = "_")

    p_bar <- df %>%
      mutate(cluster = fct_reorder(cluster, score)) %>%
      pivot_longer(cols = c(score_up, score_down), names_to = "direction", values_to = "score_dir") %>%
      mutate(direction = str_remove(direction, "score_")) %>%
      ggplot(aes(x = cluster, y = score_dir, fill = direction)) +
      geom_col() +
      scale_fill_manual(values = c(up = "firebrick", down = "steelblue")) +
      coord_flip() +
      theme_bw() +
      theme(strip.background = element_blank()) +
      labs(title = paste0("perturbation score ", file_id),
           subtitle = paste0("n DEGs (padj < ", padj_thr, ", |log2FC| > ", lfc_thr, ") / n cells"),
           x = NULL, y = "perturbation score")

    # the groups not tested are in grey
    p_umap <- df_umap %>%
      left_join(meta_all %>% select(barcode, cluster = all_of(level)), by = "barcode") %>%
      mutate(cluster = as.character(cluster)) %>%
      left_join(df %>% select(cluster, score), by = "cluster") %>%
      # plot the cells with a score on top of the not tested ones
      arrange(!is.na(score)) %>%
      ggplot(aes(x = UMAP_1, y = UMAP_2, color = score)) +
      geom_point(size = 0.2, shape = 16) +
      scale_color_viridis_c(na.value = "grey85", option = "inferno", direction = -1) +
      theme_cowplot() +
      labs(title = paste0("perturbation score ", file_id), color = "score")

    ggsave(paste0(folder_plot, "05_perturbation_score_", file_id, ".pdf"),
           plot = p_bar | p_umap, width = 20, height = max(8, 0.25 * nrow(df) + 2), limitsize = FALSE)
  })

message("\nAll done.")
