# AIM ---------------------------------------------------------------------
# run scDist (https://github.com/phillipnicol/scDist) on the final consensus object generated in 03_apply_filter_annotate_fullObject.R, using the sample as random effect.
# split from 05_DE_sc_fullObject.R to run the memory demanding scDist in a separate job.
# same comparisons and levels of 05_DE_sc_fullObject.R: diagnosis (MS vs CTRL) and pathological stage (IN/ACT vs CTRL), on cell_id and cell_id_subcluster3.
# the comparison of the scDist distance with the perturbation score of 05_DE_sc_fullObject.R is in 05c_perturbationScore_vs_scDist_fullObject.R (run it after 05 and 05b, the three scripts can be otherwise run in parallel).
# readapted from scRNA_SOD_3rd/analysis/R_analysis_R45_pixi/scr/general/24_DEG_SOD_vs_CT_perturbation_scDist_annotation.R
# NOTICE: scDist is run on the RNA assay (log normalized data layer).
# two models are run for each comparison (list_model):
# base: cond + (1|sample). Files keep the original names (read by 05c).
# adj: cond + sex + location + (1|sample), same covariates of the pseudobulk DE (06_DE_pbulk_deseq2_fullObject.R). Files have the suffix _adj.
# NOTICE: scDist takes the distance from the FIRST fixed effect (row 2 of the coefficient table), the other fixed effects only adjust the model: cond must stay first.
# the sample random effect does not correct for an imbalance of sex/location between the conditions, hence the adjusted model.
# scDist skips (with a warning) the groups where sex or location have a single level, so the adj tables can have fewer groups than the base ones.
# NOTICE: scDist converts the matrix to a dense one and fits each group independently: it is run one group at a time to cap the RAM (same results as a single call)

# setup -------------------------------------------------------------------
# scDist is not on CRAN, install it once from GitHub
# remotes::install_github("phillipnicol/scDist")

# LIBRARIES ---------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(patchwork)
library(cowplot)
library(scDist)

# PARAMETERS ------------------------------------------------------------------
options(stringsAsFactors = FALSE)
options(Seurat.object.assay.version = "v5")

# contrasts to test. the first level is the reference. keep in sync with 05_DE_sc_fullObject.R
list_contrast <- list(
  diagnosis = list(var = "diagnosis_short", levels = c("CTRL","MS"), ref = "CTRL"),
  pathology = list(var = "pathological_stage_short", levels = c("CTRL","IN","ACT"), ref = "CTRL")
)

# annotation levels to test
list_level <- c("cell_id", "cell_id_subcluster3")

# minimum number of cells per group to run the comparison
min_cells <- 10
# minimum number of samples per level of the contrast to run scDist
min_samples <- 3

# the random effect of scDist
id_sample <- "sample_id"

# models to run: additional fixed effects after cond (cond is always the first, see the NOTICE above) and suffix of the output files
list_model <- list(
  base = list(covariates = character(0), tag = ""),
  adj = list(covariates = c("sex", "location"), tag = "_adj")
)

# scDist:::condition_v_celltype_distance runs one scDist per pair of groups on the full matrix: very slow and memory demanding with many groups (cell_id_subcluster3), off by default
run_diagnostic <- FALSE

folder_table <- "../../out/table/analysis_R45_pixi/05_DE_sc/"
folder_plot <- "../../out/plot/analysis_R45_pixi/05_DE_sc/"

# read in the data --------------------------------------------------------
sobj <- readRDS("../../out/object/analysis_R45_pixi/03_sobj_integrated_cleanup_manualAnnotation_subclusterFiltered_v2.rds")

# make sure the correct default assay is loaded, should be RNA
DefaultAssay(sobj) <- "RNA"
class(sobj@assays$RNA)

# join the layers in case they are split per sample
if(length(Layers(sobj, search = "data")) > 1){
  sobj <- JoinLayers(sobj)
}

# normalized matrix (genes x cells) used by scDist, metadata and UMAP coordinates used for the plots
mat_norm <- LayerData(sobj, assay = "RNA", layer = "data")
mat_norm[1:10, 1:10]

meta_all <- sobj@meta.data %>%
  mutate(barcode = rownames(.))
df_umap <- Embeddings(sobj, "umap") %>%
  as.data.frame() %>%
  set_names(c("UMAP_1", "UMAP_2")) %>%
  rownames_to_column("barcode")

# the full object is not needed anymore, free the RAM before scDist
rm(sobj)
invisible(gc())

# run scDist --------------------------------------------------------------
df_run <- crossing(level = list_level, contrast = names(list_contrast))

# each comparison is a two-level test, the cells of the other levels are dropped. The groups are kept only if they have enough cells and samples in both levels
list_scdist <- pmap(df_run, function(level, contrast){
  ct <- list_contrast[[contrast]]
  lvl_test <- setdiff(ct$levels, ct$ref)

  pmap(crossing(lvl = lvl_test, model = names(list_model)), function(lvl, model){
    fixed_effects <- c("cond", list_model[[model]]$covariates)
    file_id <- paste0(level, "_", contrast, "_", lvl, "_vs_", ct$ref, list_model[[model]]$tag)
    message("run scDist ", file_id, ": ", paste(fixed_effects, collapse = " + "))

    # cells of the two levels of the comparison
    meta_cmp <- meta_all %>%
      filter(.data[[ct$var]] %in% c(lvl, ct$ref)) %>%
      mutate(group = as.character(.data[[level]]),
             cond = factor(.data[[ct$var]], levels = c(ct$ref, lvl)),
             sample = .data[[id_sample]])

    df_group <- meta_cmp %>%
      group_by(group) %>%
      summarise(n_1 = sum(cond == lvl),
                n_2 = sum(cond == ct$ref),
                s_1 = n_distinct(sample[cond == lvl]),
                s_2 = n_distinct(sample[cond == ct$ref]),
                .groups = "drop")

    groups_ok <- df_group %>%
      filter(n_1 >= min_cells, n_2 >= min_cells, s_1 >= min_samples, s_2 >= min_samples) %>%
      pull(group)

    if(length(groups_ok) < 2){
      message("  skip ", file_id, ": fewer than 2 groups with enough cells/samples")
      return(NULL)
    }

    meta_cmp <- meta_cmp %>%
      filter(group %in% groups_ok)

    # the first level of the fixed effect is the reference. the rownames are the barcodes
    meta_scdist <- meta_cmp[, unique(c(fixed_effects, "sample", "group"))]
    rownames(meta_scdist) <- meta_cmp$barcode

    # scDist converts the whole sparse matrix to a dense one (as.matrix + t) and then fits each group independently: to cap the RAM it is run one group at a time, so only the cells of the group are densified. The results are the same as a single call
    out_list <- map(groups_ok, function(g){
      cells_g <- rownames(meta_scdist)[meta_scdist$group == g]
      set.seed(123)
      out_g <- tryCatch({
        scDist(normalized_counts = mat_norm[, cells_g],
               meta.data = meta_scdist[cells_g, ],
               fixed.effects = fixed_effects,
               random.effects = "sample",
               clusters = "group")
      }, error = function(e){
        message("  failed ", file_id, " ", g, ": ", conditionMessage(e))
        NULL
      })
      # release the dense copies before the next group
      gc(verbose = FALSE)
      out_g
    }) %>%
      compact()

    if(length(out_list) == 0){
      return(NULL)
    }

    # merge in a single scDist-like object (a plain list, used by DistPlot and distGenes)
    out_scdist <- list(results = do.call(rbind, unname(map(out_list, "results"))),
                       vals = do.call(c, unname(map(out_list, "vals"))),
                       design = out_list[[1]]$design,
                       gene.names = out_list[[1]]$gene.names)

    # groups skipped inside scDist (single sampled level) have no row
    if(nrow(out_scdist$results) == 0){
      message("  skip ", file_id, ": no group with results")
      return(NULL)
    }

    df_res <- out_scdist$results %>%
      as.data.frame() %>%
      rownames_to_column("cluster") %>%
      left_join(df_group %>% dplyr::rename(cluster = group), by = "cluster") %>%
      mutate(level = level,
             contrast = contrast,
             comparison = paste0(lvl, "_vs_", ct$ref),
             model = model,
             design = paste(deparse(out_scdist$design), collapse = ""))

    df_res %>%
      write_tsv(paste0(folder_table, "05_scDist_", file_id, ".tsv"))

    # distance per group, next to the UMAP colored by the distance
    p_dist <- DistPlot(out_scdist) +
      theme_cowplot() +
      ggtitle(file_id)

    p_umap <- df_umap %>%
      left_join(meta_all %>% select(barcode, cluster = all_of(level)), by = "barcode") %>%
      mutate(cluster = as.character(cluster)) %>%
      left_join(df_res %>% select(cluster, dist = `Dist.`), by = "cluster") %>%
      arrange(!is.na(dist)) %>%
      ggplot(aes(x = UMAP_1, y = UMAP_2, color = dist)) +
      geom_point(size = 0.2, shape = 16) +
      scale_color_viridis_c(na.value = "grey85", option = "inferno", direction = -1) +
      theme_cowplot() +
      labs(title = paste0("scDist ", file_id), color = "Dist.")

    ggsave(paste0(folder_plot, "05_scDist_", file_id, ".pdf"),
           plot = p_dist | p_umap, width = 20, height = max(8, 0.3 * nrow(df_res) + 2), limitsize = FALSE)

    # genes driving the distance, one plot per group
    list_genes <- map(rownames(out_scdist$results) %>% set_names(), function(cl){
      distGenes(out_scdist, cluster = cl, num_genes = 5) +
        ggtitle(cl)
    })
    ggsave(paste0(folder_plot, "05_scDist_distGenes_", file_id, ".pdf"),
           plot = wrap_plots(list_genes, ncol = 3),
           width = 20, height = 4 * ceiling(length(list_genes) / 3), limitsize = FALSE)

    # distance between the groups vs the distance due to the contrast (red). NOTICE: full matrix, very memory demanding
    if(run_diagnostic && model == "base"){
      set.seed(123)
      out_ctdist <- scDist:::condition_v_celltype_distance(normalized_counts = mat_norm[, rownames(meta_scdist)],
                                                           meta.data = meta_scdist,
                                                           fixed.effects = "cond",
                                                           random.effects = "sample",
                                                           clusters = "group")
      out_ctdist$D %>%
        as.data.frame() %>%
        rownames_to_column("cluster") %>%
        write_tsv(paste0(folder_table, "05_scDist_condition_v_celltype_distance_", file_id, ".tsv"))

      ggsave(paste0(folder_plot, "05_scDist_condition_v_celltype_distance_", file_id, ".pdf"),
             plot = out_ctdist$plot + theme_cowplot() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle(file_id),
             width = 10, height = 6, limitsize = FALSE)
    }

    df_res
  }) %>%
    bind_rows()
})

df_scdist_all <- bind_rows(list_scdist)

# one table per model: the base one keeps the original name (read by 05c), the adjusted one has the suffix _adj
iwalk(list_model, function(m, name_model){
  df_scdist_all %>%
    filter(model == name_model) %>%
    write_tsv(paste0(folder_table, "05_scDist_allComparisons", m$tag, ".tsv"))
})

message("\nAll done.")
