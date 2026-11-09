# AIM ---------------------------------------------------------------------
# pseudobulk DE analysis (DESeq2) on the final consensus object generated in 03_apply_filter_annotate_fullObject.R.
# run for each macro population (cell_id). the consensus subclusters (cell_id_subcluster3) are not analysed for the moment.
# the raw counts are aggregated per sample with Seurat::AggregateExpression
# comparisons: diagnosis (MS vs CTRL) and pathological stage (IN/ACT vs CTRL).
# the model accounts for the location and the sex: ~ location + sex + treat
# the DEGs are called on the ashr shrunk LFC (padj < th_padj and |log2FC| > th_lfc)
# besides the per unit diagnostic plots (PCA, p-value histogram, volcano, MA plot, heatmap), at the end of the script:
# 1. barplot of the number of DEGs per cell type, stacked by direction
# 2. upsetplot of the DEGs shared between the single cell (05_DE_sc_fullObject.R) and the pseudobulk analysis, per cell type
#    NOTICE: a gene is counted as shared only if it is called in the same direction by the two approaches (the element of the sets is gene|direction). The same thresholds (th_padj, th_lfc) are applied to the single cell tables
# 3. correlation of the log2FC and of the adjusted p-value between single cell and pseudobulk, on the genes tested by both approaches
# readapted also from scRNA_SOD_3rd/analysis/R_analysis_R45_pixi/scr/general/28_pseudobulk_SOD_vs_CT_DESeq2_anno_l0.R
# readapted from R_analysis_R44/scr/30_DE_analysis_pbulk_deseq2_IMMUNE_MsVsCtrl.R and 35_DE_analysis_pbulk_deseq2_IMMUNE_ActiveVsCtrl.R
# the original routine is a readaptation of different resources
# https://github.com/sib-swiss/single-cell-training/
# https://hbctraining.github.io/scRNA-seq/lessons/pseudobulk_DESeq2_scrnaseq.html

# LIBRARIES ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(GGally)
library(cowplot)
library(ComplexHeatmap)
library(scales)
library(circlize)
library(DESeq2)
library(RNAseqQC)
library(limma)
library(ashr)
library(magick)
library(Matrix)
library(ggrepel)

# PARAMETERS ------------------------------------------------------------------
options(stringsAsFactors = FALSE)
options(Seurat.object.assay.version = "v5")

# contrasts to test. the first level is the reference
list_contrast <- list(
  diagnosis = list(var = "diagnosis_short", levels = c("CTRL","MS"), ref = "CTRL"),
  pathology = list(var = "pathological_stage_short", levels = c("CTRL","IN","ACT"), ref = "CTRL")
)

# drop the pseudobulk samples built from too few cells
min_cells_sample <- 10
# minimum number of samples per level of the contrast to run the test
min_samples <- 3

# thresholds for the DEGs
th_padj <- 0.05
th_lfc <- 0.5

folder_table <- "../../out/table/analysis_R45_pixi/06_pbulk/"
folder_plot <- "../../out/plot/analysis_R45_pixi/06_pbulk/"
folder_object <- "../../out/object/analysis_R45_pixi/06_pbulk/"
# single cell DE tables of 05_DE_sc_fullObject.R
folder_table_sc <- "../../out/table/analysis_R45_pixi/05_DE_sc/"

# make the labels safe to be used in the file names
make_safe <- function(x){
  x %>%
    str_replace_all(pattern = "[^A-Za-z0-9]+", replacement = "_") %>%
    str_remove(pattern = "_$")
}

# read in the data --------------------------------------------------------
sobj <- readRDS("../../out/object/analysis_R45_pixi/03_sobj_integrated_cleanup_manualAnnotation_subclusterFiltered_v2.rds")
DefaultAssay(sobj) <- "RNA"

meta_all <- sobj@meta.data %>%
  mutate(across(all_of(map_chr(list_contrast, "var")), as.character))

# sample level metadata. sample_id maps 1:1 to the clinical covariates
LUT_sample <- meta_all %>%
  distinct(sample_id, sex, age, weight, csf, apoe, location,
           diagnosis_short, pathological_stage_short)

# units to test: the macro populations (cell_id) only.
# the consensus subclusters (cell_id_subcluster3) are not analysed for the moment
df_unit <- data.frame(unit_col = "cell_id", unit_value = sort(unique(as.character(meta_all$cell_id)))) %>%
  mutate(unit_id = paste0(unit_col, "_", make_safe(unit_value)))

df_unit

# functions -----------------------------------------------------------------
# AggregateExpression sanitizes the group names: prefix "g" if the name starts with a digit and "_" replaced by "-".
# replicated here to match the sample_id with the column names of the aggregated matrix
fix_sample_names <- function(x){
  case_when(str_detect(x, "^\\d") ~ paste0("g", x), TRUE ~ x) %>%
    str_replace_all("_", "-")
}

## pseudobulk counts and sample metadata for one unit (subset of the seurat object) and one contrast
build_pseudobulk <- function(sobj_unit, ct){
  # n cells and average qc per sample, drop the samples with too few cells
  meta_sample <- sobj_unit@meta.data %>%
    group_by(sample_id) %>%
    summarise(n_cells = n(),
              avg_percent.mt = mean(percent.mt),
              avg_percent.ribo = mean(percent.ribo)) %>%
    mutate(sample_id = as.character(sample_id)) %>%
    filter(n_cells >= min_cells_sample)

  if(nrow(meta_sample) == 0){
    return(NULL)
  }

  # sum of the raw counts per sample
  agg <- AggregateExpression(sobj_unit, assays = "RNA", group.by = "sample_id", return.seurat = FALSE)$RNA
  counts <- agg[, fix_sample_names(meta_sample$sample_id), drop = FALSE] %>%
    as.matrix() %>%
    round()
  colnames(counts) <- meta_sample$sample_id
  storage.mode(counts) <- "integer"

  # sample metadata
  col_data <- meta_sample %>%
    left_join(LUT_sample, by = "sample_id") %>%
    mutate(M_reads = colSums(counts) / 10^6,
           treat = factor(.data[[ct$var]], levels = ct$levels),
           location = factor(location, levels = c("cervical","thoracic","lumbar")),
           sex = factor(sex)) %>%
    as.data.frame()
  rownames(col_data) <- col_data$sample_id

  list(counts = counts, col_data = droplevels(col_data))
}

## DESeq2 model for one unit: fixed design ~ location + sex + treat, low expressed genes dropped with filterByExpr
fit_pseudobulk_dds <- function(counts, col_data){
  design <- model.matrix(~ location + sex + treat, data = col_data)
  colnames(design)[1] <- "intercept"

  dds <- DESeqDataSetFromMatrix(countData = counts, colData = col_data, design = design)
  dds <- dds[edgeR::filterByExpr(counts(dds), group = col_data$treat), ]
  dds <- DESeq(dds)

  list(dds = dds, design = design)
}

## sample level QC plots for one unit: total reads, sample clustering, PCA. returns the variance stabilized object
plot_pseudobulk_qc <- function(dds, col_data, ct, file_id){
  # total reads per sample
  p_reads <- colSums(counts(dds)) %>%
    data.frame(tot_reads = .) %>%
    rownames_to_column("sample") %>%
    ggplot(aes(x = sample, y = tot_reads)) + geom_col() +
    theme_bw() +
    theme(axis.text.x = element_text(hjust = 1, angle = 90)) +
    ggtitle(file_id)
  ggsave(paste0(folder_plot, "06_pBulk_totReads_", file_id, ".pdf"), plot = p_reads, width = 10, height = 6)

  # scale the data
  if(nrow(dds) > 1000){
    vsd <- vst(dds, blind = TRUE)
  } else {
    vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
  }

  # clustering of the samples. the seed controls the random annotation colors
  pdf(paste0(folder_plot, "06_heatmap_SampleCluster_filterExp_", file_id, ".pdf"), width = 12, height = 6)
  set.seed(21)
  hm <- plot_sample_clustering(vsd,
                               anno_vars = c("treat","location","sex","diagnosis_short","pathological_stage_short","age","M_reads"),
                               distance = "euclidean")
  draw(hm, heatmap_legend_side = "left", annotation_legend_side = "left", padding = unit(c(2, 2, 2, 30), "mm"))
  dev.off()

  # PCA on the 500 genes with the highest variance
  rv <- rowVars(assay(vsd), useNames = TRUE)
  select_var <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
  pca <- prcomp(t(assay(vsd)[select_var, ]))
  percent_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  n_pc <- min(10, ncol(pca$x))

  df_pca <- pca$x[, 1:n_pc, drop = FALSE] %>%
    data.frame() %>%
    rownames_to_column("sample_id") %>%
    left_join(col_data, by = "sample_id")

  p_pca1 <- df_pca %>%
    ggplot(aes(x = PC1, y = PC2, col = location)) +
    geom_point(size = 3, alpha = 0.6) +
    scale_x_continuous(expand = expansion(mult = 0.1)) +
    theme_bw() + xlab(paste0("PC1: ", percent_var[1], "% variance")) + ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    ggtitle(file_id)
  ggsave(paste0(folder_plot, "06_PCA_pseudobulk_filterExp_", file_id, ".pdf"), plot = p_pca1, width = 7, height = 6)

  p_pca2 <- df_pca %>%
    ggplot(aes(x = PC1, y = PC2, col = treat, label = sample_id)) +
    geom_point(size = 3, alpha = 0.6) +
    ggrepel::geom_text_repel(show.legend = FALSE) +
    scale_x_continuous(expand = expansion(mult = 0.1)) +
    theme_bw() + xlab(paste0("PC1: ", percent_var[1], "% variance")) + ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    labs(col = ct$var) +
    ggtitle(file_id)
  ggsave(paste0(folder_plot, "06_PCA_pseudobulk_filterExp2_", file_id, ".pdf"), plot = p_pca2, width = 7, height = 6)

  # more PC by condition
  p_panel <- df_pca %>%
    ggpairs(columns = paste0("PC", 1:n_pc), ggplot2::aes(colour = treat), upper = "blank") +
    theme_bw() + theme(strip.background = element_blank(), axis.text.x = element_text(hjust = 1, angle = 45))
  ggsave(paste0(folder_plot, "06_PCA_pseudobulk_panel_", file_id, ".pdf"), plot = p_panel, width = 2 * n_pc, height = 2 * n_pc)

  # PC scores by sample metadata
  df_pc_long <- df_pca %>%
    dplyr::select(sample_id, all_of(paste0("PC", 1:n_pc))) %>%
    pivot_longer(names_to = "var_2", values_to = "value_2", -sample_id) %>%
    mutate(var_2 = str_remove(var_2, "PC") %>% str_pad(side = "left", width = 2, pad = "0") %>% paste0("PC", .))

  p_cat <- df_pca %>%
    dplyr::select(sample_id, treat, sex, location, diagnosis_short, pathological_stage_short) %>%
    mutate(across(-sample_id, as.character)) %>%
    pivot_longer(names_to = "var_1", values_to = "value_1", -sample_id) %>%
    left_join(df_pc_long, by = "sample_id", relationship = "many-to-many") %>%
    mutate(comparison = paste0(var_1, "_vs_", var_2)) %>%
    ggplot(aes(x = value_1, y = value_2)) +
    facet_wrap(~comparison, scales = "free", ncol = n_pc) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitter(width = 0.2)) +
    theme_bw() + theme(strip.background = element_blank(), axis.text.x = element_text(hjust = 1, angle = 45))
  ggsave(paste0(folder_plot, "06_panel_metadata_PC_pbulk_cat_", file_id, ".pdf"), plot = p_cat, width = 3 * n_pc, height = 18, limitsize = FALSE)

  p_cont <- df_pca %>%
    dplyr::select(sample_id, age, weight, csf, avg_percent.mt, avg_percent.ribo, M_reads, n_cells) %>%
    mutate(across(-sample_id, ~ suppressWarnings(as.numeric(.x)))) %>%
    pivot_longer(names_to = "var_1", values_to = "value_1", -sample_id) %>%
    left_join(df_pc_long, by = "sample_id", relationship = "many-to-many") %>%
    mutate(comparison = paste0(var_1, "_vs_", var_2)) %>%
    ggplot(aes(x = value_1, y = value_2)) +
    facet_wrap(~comparison, scales = "free", ncol = n_pc) +
    geom_smooth(method = "lm") +
    geom_point() +
    theme_bw() + theme(strip.background = element_blank(), axis.text.x = element_text(hjust = 1, angle = 45))
  ggsave(paste0(folder_plot, "06_panel_metadata_PC_pbulk_cont_", file_id, ".pdf"), plot = p_cont, width = 3 * n_pc, height = 21, limitsize = FALSE)

  # location and sex removal, for visualisation only ----------------------
  # mimics on the plots what location + sex do in the DE design (DE uses the raw counts, never these values).
  # the design ~ treat is passed to protect the effect of the contrast, otherwise the treatment signal is partly removed when location or sex are unbalanced between the groups.
  # a covariate with a single level can not be removed and is skipped
  cd_vsd <- as.data.frame(colData(vsd))
  list_batch <- list(location = cd_vsd$location, sex = cd_vsd$sex) %>%
    keep(~ nlevels(.x) >= 2)

  if(length(list_batch) > 0){
    vsd_batch <- vsd
    assay(vsd_batch) <- limma::removeBatchEffect(assay(vsd),
                                                 batch = list_batch[[1]],
                                                 batch2 = if(length(list_batch) > 1) list_batch[[2]] else NULL,
                                                 design = model.matrix(~ treat, data = cd_vsd))

    # clustering of the samples after the removal
    pdf(paste0(folder_plot, "06_heatmap_SampleCluster_batchRemoved_", file_id, ".pdf"), width = 12, height = 6)
    set.seed(21)
    hm_batch <- plot_sample_clustering(vsd_batch, anno_vars = c("treat","location","sex"), distance = "euclidean")
    draw(hm_batch, heatmap_legend_side = "left", annotation_legend_side = "left", padding = unit(c(2, 2, 2, 30), "mm"))
    dev.off()

    # PCA before (top) and after (bottom) the removal, coloured by treat (left) and location (right)
    plot_pca <- function(v, col_var, title){
      df <- plotPCA(v, intgroup = c("treat","location","sex"), returnData = TRUE)
      pv <- round(100 * attr(df, "percentVar"))
      ggplot(df, aes(x = PC1, y = PC2, col = .data[[col_var]], shape = sex)) +
        geom_point(size = 3, alpha = 0.7) +
        xlab(paste0("PC1: ", pv[1], "% variance")) + ylab(paste0("PC2: ", pv[2], "% variance")) +
        theme_bw() + ggtitle(title)
    }

    p_batch <- plot_grid(plot_pca(vsd, "treat", paste0(file_id, " - vst")),
                         plot_pca(vsd, "location", paste0(file_id, " - vst")),
                         plot_pca(vsd_batch, "treat", paste0(file_id, " - ", paste(names(list_batch), collapse = " + "), " removed")),
                         plot_pca(vsd_batch, "location", paste0(file_id, " - ", paste(names(list_batch), collapse = " + "), " removed")),
                         ncol = 2)
    ggsave(paste0(folder_plot, "06_PCA_pseudobulk_batchRemoved_", file_id, ".pdf"), plot = p_batch, width = 14, height = 11)
  }

  vsd
}

## one contrast per tested level vs the reference. returns the unshrunk and the ashr shrunk tables, one row per gene and level
test_contrasts <- function(dds, design, lvl_test, ref){
  list_res <- map(lvl_test, function(lvl){
    res <- results(dds, contrast = makeContrasts(contrasts = paste0("treat", lvl), levels = design)[, 1], alpha = th_padj)
    list(res = res, res_shr = lfcShrink(dds, res = res, type = "ashr"))
  }) %>%
    set_names(paste0(lvl_test, "_vs_", ref))

  tidy_res <- function(x){
    x %>%
      data.frame() %>%
      rownames_to_column(var = "symbol") %>%
      arrange(pvalue)
  }

  list(df_res = map(list_res, ~ tidy_res(.x$res)) %>% bind_rows(.id = "condVsRef"),
       df_res_shr = map(list_res, ~ tidy_res(.x$res_shr)) %>% bind_rows(.id = "condVsRef"))
}

## diagnostic plots of the DE results for one unit: p-value histogram, volcano, DEG heatmap, MA plot (unshrunk and shrunk), dispersion
plot_de <- function(dds, vsd, df_res, df_res_shr, col_data, n_test, file_id){
  # histogram of the p values, excluding genes with very small counts
  p_hist <- df_res %>%
    dplyr::filter(baseMean > 1) %>%
    ggplot(aes(x = pvalue)) + geom_histogram(breaks = 0:20/20) +
    facet_wrap(~condVsRef) +
    theme_bw() +
    theme(strip.background = element_blank()) +
    ggtitle(file_id)
  ggsave(paste0(folder_plot, "06_histogram_pvalue_pseudobulk_filterExp_", file_id, ".pdf"), plot = p_hist, width = 3 * n_test + 1, height = 4)

  # volcano, on the unshrunk (res) and on the shrunk (shr) estimates
  walk2(list(df_res, df_res_shr), c("", "shr_"), function(df, tag){
    df_plot <- df %>%
      mutate(col = ifelse(((padj < th_padj) & abs(log2FoldChange) > th_lfc & !is.na(symbol)), yes = 1, no = 0)) %>%
      filter(!is.na(col))

    p_volcano <- df_plot %>%
      ggplot(aes(x = log2FoldChange, y = -log(padj))) +
      geom_point(data = df_plot[df_plot$col == 0, ], aes(col = factor(col)), alpha = 0.05) +
      geom_point(data = df_plot[df_plot$col == 1, ], aes(col = factor(col)), alpha = 0.5) +
      geom_vline(xintercept = c(-th_lfc, th_lfc), col = "red", linetype = "dashed") +
      geom_hline(yintercept = (-log(th_padj)), col = "red", linetype = "dashed") +
      scale_color_manual(values = c("0" = "black", "1" = "red")) +
      ggrepel::geom_text_repel(
        data = df_plot %>%
          group_by(condVsRef) %>%
          arrange(padj) %>%
          dplyr::slice(1:30) %>%
          dplyr::filter(abs(log2FoldChange) > th_lfc) %>%
          dplyr::filter(padj < th_padj),
        aes(label = symbol), segment.alpha = 0.4) +
      facet_wrap(~condVsRef) +
      theme_bw() +
      theme(strip.background = element_blank()) +
      theme(legend.position = "none") +
      ggtitle(file_id)
    ggsave(paste0(folder_plot, "06_vulcano_plot_pseudobulk_filterExp_", tag, file_id, ".pdf"), plot = p_volcano, width = 6 * n_test + 1, height = 7)
  })

  # heatmap of the DEGs and MA plot, on the unshrunk (res) and on the shrunk (shr) estimates
  walk2(list(df_res, df_res_shr), c("", "shr_"), function(df, tag){
    DEG_2 <- df %>%
      dplyr::filter(padj < th_padj, abs(log2FoldChange) > th_lfc) %>%
      pull(symbol) %>%
      unique()

    if(length(DEG_2) > 1){
      mat_shr <- assay(vsd)[rownames(vsd) %in% DEG_2, , drop = FALSE]
      mat2_shr <- (mat_shr - rowMeans(mat_shr)) / rowSds(mat_shr, useNames = TRUE)

      column_ha_shr <- HeatmapAnnotation(treat = as.character(col_data[colnames(mat2_shr), "treat"]),
                                         origin = as.character(col_data[colnames(mat2_shr), "location"]),
                                         col = list(origin = c("lumbar" = "gray",
                                                               "thoracic" = "gray30",
                                                               "cervical" = "black")))

      ht2_shr <- Heatmap(mat2_shr, show_column_names = TRUE, raster_by_magick = TRUE, show_row_dend = FALSE, use_raster = TRUE,
                         name = "exp",
                         column_title = file_id,
                         row_names_gp = gpar(fontsize = ifelse(length(DEG_2) > 100, 3, 8)),
                         top_annotation = column_ha_shr, show_row_names = TRUE)
      pdf(paste0(folder_plot, "06_heatmap_DEG_plot_pseudobulk_filterExp_", tag, file_id, ".pdf"), width = 10, height = 8)
      draw(ht2_shr, heatmap_legend_side = "left", annotation_legend_side = "left", padding = unit(c(30, 2, 2, 2), "mm"))
      dev.off()
    }

    p_ma <- df %>%
      filter(!is.na(padj)) %>%
      mutate(color = factor(ifelse(padj < th_padj, yes = 1, no = 0), levels = c(0, 1))) %>%
      arrange(color) %>%
      ggplot(aes(x = baseMean, y = log2FoldChange, col = color)) +
      geom_point(alpha = 0.2) +
      scale_x_log10() +
      scale_color_manual(values = c("0" = "gray", "1" = "red"), drop = FALSE) +
      facet_wrap(~condVsRef) +
      theme_bw() +
      theme(strip.background = element_blank(), legend.position = "none") +
      ggtitle(file_id)
    ggsave(paste0(folder_plot, "06_MA_plot_pseudobulk_filterExp_", tag, file_id, ".pdf"), plot = p_ma, width = 6 * n_test + 1, height = 6)
  })

  # dispersion
  pdf(paste0(folder_plot, "06_dispersion_plot_pseudobulk_filterExp_", file_id, ".pdf"), width = 5, height = 5)
  plotDispEsts(dds)
  dev.off()
}

# main loop: unit x contrast ----------------------------------------------
df_run <- df_unit %>%
  crossing(contrast = names(list_contrast))

list_summary <- pmap(df_run, function(unit_col, unit_value, unit_id, contrast){
  ct <- list_contrast[[contrast]]
  file_id <- paste0(unit_id, "_", contrast)
  message("== ", file_id)
  df_tag <- data.frame(unit_id = unit_id, unit_col = unit_col, cluster = unit_value, contrast = contrast)

  tryCatch({
    # cells of the unit with an available classification for the contrast
    cells <- rownames(meta_all)[as.character(meta_all[[unit_col]]) == unit_value & !is.na(meta_all[[ct$var]])]
    if(length(cells) == 0){
      return(data.frame(df_tag, status = "skip: no cells"))
    }

    pb <- build_pseudobulk(subset(sobj, cells = cells), ct)
    if(is.null(pb)){
      return(data.frame(df_tag, status = "skip: no samples"))
    }

    # levels with enough samples to be tested against the reference
    n_level <- table(pb$col_data$treat)
    n_level_txt <- paste(names(n_level), n_level, sep = ":", collapse = "|")
    lvl_test <- setdiff(names(n_level), ct$ref)
    lvl_test <- lvl_test[n_level[lvl_test] >= min_samples]

    if(!(ct$ref %in% names(n_level)) || n_level[[ct$ref]] < min_samples || length(lvl_test) == 0){
      message("skip ", file_id, ": ", n_level_txt)
      return(data.frame(df_tag, n_samples = n_level_txt, status = "skip: too few samples"))
    }

    # keep only the samples of the levels to be tested
    pb$col_data <- droplevels(pb$col_data[pb$col_data$treat %in% c(ct$ref, lvl_test), ])
    pb$counts <- pb$counts[, rownames(pb$col_data)]

    saveRDS(pb$counts, file = paste0(folder_object, "06_counts_pBulk_", file_id, ".rds"))
    saveRDS(pb$col_data, file = paste0(folder_object, "06_colData_pBulk_", file_id, ".rds"))

    # model, QC and DE
    fit <- fit_pseudobulk_dds(pb$counts, pb$col_data)
    saveRDS(fit$dds, paste0(folder_object, "06_ddsHTSeq_pseudobulk_filterExp_", file_id, ".rds"))

    vsd <- plot_pseudobulk_qc(fit$dds, pb$col_data, ct, file_id)

    res <- test_contrasts(fit$dds, fit$design, lvl_test, ct$ref)
    df_res <- res$df_res %>% mutate(unit_col = unit_col, cluster = unit_value, contrast = contrast)
    df_res_shr <- res$df_res_shr %>% mutate(unit_col = unit_col, cluster = unit_value, contrast = contrast)

    write_tsv(df_res, paste0(folder_table, "06_DE_pseudobulk_filterExp_", file_id, ".tsv"))
    write_tsv(df_res_shr, paste0(folder_table, "06_DE_pseudobulk_filterExp_shr_", file_id, ".tsv"))

    plot_de(fit$dds, vsd, df_res, df_res_shr, pb$col_data, length(lvl_test), file_id)

    # summary of the run. n_up/n_down: shrunk LFC (used by the barplot and the comparison with the single cell), *_unshr: unshrunk LFC
    df_n_shr <- df_res_shr %>%
      group_by(condVsRef) %>%
      summarise(n_genes = n(),
                n_up = sum(padj < th_padj & log2FoldChange > th_lfc, na.rm = TRUE),
                n_down = sum(padj < th_padj & log2FoldChange < -th_lfc, na.rm = TRUE))

    df_n_unshr <- df_res %>%
      group_by(condVsRef) %>%
      summarise(n_up_unshr = sum(padj < th_padj & log2FoldChange > th_lfc, na.rm = TRUE),
                n_down_unshr = sum(padj < th_padj & log2FoldChange < -th_lfc, na.rm = TRUE))

    left_join(df_n_shr, df_n_unshr, by = "condVsRef") %>%
      mutate(unit_id = unit_id, unit_col = unit_col, cluster = unit_value, contrast = contrast,
             n_samples = n_level_txt,
             design = paste(colnames(fit$design), collapse = "+"),
             status = "done")
  }, error = function(e){
    message("failed ", file_id, ": ", conditionMessage(e))
    # close any device left open by the failure
    while(dev.cur() > 1) dev.off()
    data.frame(df_tag, status = paste0("error: ", conditionMessage(e)))
  })
})

df_summary <- bind_rows(list_summary) %>%
  relocate(unit_id, unit_col, cluster, contrast, status)

df_summary %>%
  write_tsv(paste0(folder_table, "06_DE_pbulk_summary.tsv"))

df_summary %>%
  print(n = Inf)

# barplot of the number of DEGs ---------------------------------------------
# only the units successfully run. The tables are read back from the output folder
df_done <- df_summary %>%
  filter(status == "done") %>%
  distinct(unit_id, unit_col, cluster, contrast)

df_pbulk_shr <- pmap(df_done, function(unit_id, unit_col, cluster, contrast){
  read_tsv(paste0(folder_table, "06_DE_pseudobulk_filterExp_shr_", unit_id, "_", contrast, ".tsv"),
           col_types = cols(cluster = col_character(), .default = col_guess()))
}) %>%
  bind_rows()

df_summary %>%
  filter(status == "done") %>%
  split(f = list(.$unit_col, .$contrast, .$condVsRef), drop = TRUE) %>%
  walk(function(df){
    file_id <- paste(unique(df$unit_col), unique(df$contrast), unique(df$condVsRef), sep = "_")

    p_bar <- df %>%
      mutate(cluster = fct_reorder(cluster, n_up + n_down)) %>%
      pivot_longer(cols = c(n_up, n_down), names_to = "direction", values_to = "n") %>%
      mutate(direction = str_remove(direction, "n_")) %>%
      ggplot(aes(x = cluster, y = n, fill = direction)) +
      geom_col() +
      scale_fill_manual(values = c(up = "firebrick", down = "steelblue")) +
      coord_flip() +
      theme_bw() +
      theme(strip.background = element_blank()) +
      labs(title = paste0("pseudobulk DESeq2 ", file_id),
           subtitle = paste0("n DEGs (padj < ", th_padj, ", |shrunk log2FC| > ", th_lfc, ")"),
           x = NULL, y = "n DEGs")
    ggsave(paste0(folder_plot, "06_barplot_nDEG_pseudobulk_", file_id, ".pdf"),
           plot = p_bar, width = 8, height = max(5, 0.25 * nrow(df) + 2), limitsize = FALSE)
  })

# single cell vs pseudobulk -----------------------------------------------
# single cell tables of 05_DE_sc_fullObject.R (one per level/contrast/comparison, the summary table is not matched by the pattern)
files_sc <- list.files(folder_table_sc, pattern = "^05_DE_sc_cell_id.*_vs_.*\\.tsv$", full.names = TRUE)

if(length(files_sc) == 0){
  message("no single cell table found in ", folder_table_sc, ": skip the comparison with the pseudobulk")
} else {
  key <- c("unit_col", "contrast", "condVsRef", "cluster")

  df_sc <- map(files_sc, function(f){
    read_tsv(f, col_types = cols(cluster = col_character(), .default = col_guess()))
  }) %>%
    bind_rows() %>%
    dplyr::rename(unit_col = level, condVsRef = comparison)

  # groups tested by both approaches
  groups_both <- inner_join(distinct(df_sc, across(all_of(key))),
                            distinct(df_pbulk_shr, across(all_of(key))),
                            by = key)

  # DEGs with the same thresholds, element = gene|direction
  df_deg_both <- bind_rows(
    df_sc %>%
      filter(p_val_adj < th_padj, abs(avg_log2FC) > th_lfc) %>%
      mutate(direction = ifelse(avg_log2FC > 0, "up", "down"), method = "single_cell") %>%
      select(all_of(key), gene, direction, method),
    df_pbulk_shr %>%
      filter(padj < th_padj, abs(log2FoldChange) > th_lfc) %>%
      mutate(direction = ifelse(log2FoldChange > 0, "up", "down"), method = "pseudobulk") %>%
      select(all_of(key), gene = symbol, direction, method)
  ) %>%
    semi_join(groups_both, by = key) %>%
    mutate(element = paste(gene, direction, sep = "|"))

  df_deg_both %>%
    write_tsv(paste0(folder_table, "06_DEG_singleCell_and_pseudobulk.tsv"))

  # overlap table, one row per group
  df_overlap <- groups_both %>%
    rowwise() %>%
    mutate(n_single_cell = sum(df_deg_both$method == "single_cell" & df_deg_both$unit_col == unit_col & df_deg_both$contrast == contrast & df_deg_both$condVsRef == condVsRef & df_deg_both$cluster == cluster),
           n_pseudobulk = sum(df_deg_both$method == "pseudobulk" & df_deg_both$unit_col == unit_col & df_deg_both$contrast == contrast & df_deg_both$condVsRef == condVsRef & df_deg_both$cluster == cluster),
           n_shared = {
             sel <- df_deg_both$unit_col == unit_col & df_deg_both$contrast == contrast & df_deg_both$condVsRef == condVsRef & df_deg_both$cluster == cluster
             length(intersect(df_deg_both$element[sel & df_deg_both$method == "single_cell"], df_deg_both$element[sel & df_deg_both$method == "pseudobulk"]))
           },
           jaccard = n_shared / max(n_single_cell + n_pseudobulk - n_shared, 1)) %>%
    ungroup()

  df_overlap %>%
    write_tsv(paste0(folder_table, "06_overlap_DEG_singleCell_vs_pseudobulk.tsv"))

  # one pdf per level/contrast/comparison, one upsetplot per group (one page each)
  df_deg_both %>%
    split(f = list(.$unit_col, .$contrast, .$condVsRef), drop = TRUE) %>%
    walk(function(df){
      file_id <- paste(unique(df$unit_col), unique(df$contrast), unique(df$condVsRef), sep = "_")
      pdf(paste0(folder_plot, "06_upset_DEG_singleCell_vs_pseudobulk_", file_id, ".pdf"), width = 6, height = 4)
      walk(sort(unique(df$cluster)), function(cl){
        list_set <- list(single_cell = df %>% filter(cluster == cl, method == "single_cell") %>% pull(element),
                         pseudobulk = df %>% filter(cluster == cl, method == "pseudobulk") %>% pull(element))
        comb_mat <- make_comb_mat(list_set)
        # the number of genes is added on top of the intersection bars and at the end of the set bars
        draw(UpSet(comb_mat, comb_order = order(comb_size(comb_mat), decreasing = TRUE),
                   top_annotation = upset_top_annotation(comb_mat, add_numbers = TRUE),
                   right_annotation = upset_right_annotation(comb_mat, add_numbers = TRUE),
                   column_title = paste0(file_id, " - ", cl)))
      })
      dev.off()
    })

  # correlation single cell vs pseudobulk ------------------------------------
  # genes present in both full tables (no threshold). NOTICE: the single cell table only contains the genes returned by FindMarkers (min.pct filter), so the correlation is computed on the genes tested by both approaches
  # the pseudobulk LFC is the ashr shrunk one. The significance is -log10 of the adjusted p-value, floored to avoid Inf for p_val_adj == 0
  padj_floor <- 1e-300

  df_cor <- df_sc %>%
    semi_join(groups_both, by = key) %>%
    select(all_of(key), gene, avg_log2FC, p_val_adj) %>%
    inner_join(df_pbulk_shr %>% select(all_of(key), gene = symbol, log2FoldChange, padj), by = c(key, "gene")) %>%
    filter(!is.na(avg_log2FC), !is.na(log2FoldChange), !is.na(p_val_adj), !is.na(padj)) %>%
    mutate(neglog10_padj_sc = -log10(pmax(p_val_adj, padj_floor)),
           neglog10_padj_pbulk = -log10(pmax(padj, padj_floor)))

  df_cor_summary <- df_cor %>%
    group_by(across(all_of(key))) %>%
    filter(n() >= 3) %>%
    summarise(n_genes = n(),
              pearson_lfc = cor(avg_log2FC, log2FoldChange, method = "pearson"),
              spearman_lfc = cor(avg_log2FC, log2FoldChange, method = "spearman"),
              pearson_padj = cor(neglog10_padj_sc, neglog10_padj_pbulk, method = "pearson"),
              spearman_padj = cor(neglog10_padj_sc, neglog10_padj_pbulk, method = "spearman"),
              .groups = "drop")

  df_cor_summary %>%
    write_tsv(paste0(folder_table, "06_correlation_singleCell_vs_pseudobulk.tsv"))

  df_cor_summary %>%
    print(n = Inf)

  df_lab_cor <- df_cor_summary %>%
    mutate(lab_lfc = paste0("r = ", round(pearson_lfc, 2), "\nrho = ", round(spearman_lfc, 2), "\nn = ", n_genes),
           lab_padj = paste0("r = ", round(pearson_padj, 2), "\nrho = ", round(spearman_padj, 2), "\nn = ", n_genes))

  # one pair of scatter plots per level/contrast/comparison, faceted by group
  df_cor %>%
    split(f = list(.$unit_col, .$contrast, .$condVsRef), drop = TRUE) %>%
    walk(function(df){
      file_id <- paste(unique(df$unit_col), unique(df$contrast), unique(df$condVsRef), sep = "_")
      df_lab <- df_lab_cor %>%
        semi_join(df %>% distinct(across(all_of(key))), by = key)
      n_cluster <- n_distinct(df$cluster)
      n_col <- ceiling(sqrt(n_cluster))
      n_row <- ceiling(n_cluster / n_col)

      p_lfc <- df %>%
        ggplot(aes(x = avg_log2FC, y = log2FoldChange)) +
        geom_hline(yintercept = 0, color = "grey70") +
        geom_vline(xintercept = 0, color = "grey70") +
        geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
        geom_point(alpha = 0.2, size = 0.7) +
        geom_text(data = df_lab, aes(label = lab_lfc), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1, size = 3, inherit.aes = FALSE) +
        facet_wrap(~cluster, scales = "free", ncol = n_col) +
        theme_bw() +
        theme(strip.background = element_blank()) +
        labs(title = paste0("log2FC single cell vs pseudobulk - ", file_id),
             subtitle = "r = Pearson, rho = Spearman",
             x = "single cell avg_log2FC", y = "pseudobulk log2FoldChange (shrunk)")
      ggsave(paste0(folder_plot, "06_scatter_log2FC_singleCell_vs_pseudobulk_", file_id, ".pdf"),
             plot = p_lfc, width = 4 * n_col + 1, height = 3.5 * n_row + 1, limitsize = FALSE)

      p_padj <- df %>%
        ggplot(aes(x = neglog10_padj_sc, y = neglog10_padj_pbulk)) +
        geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
        geom_point(alpha = 0.2, size = 0.7) +
        geom_text(data = df_lab, aes(label = lab_padj), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1, size = 3, inherit.aes = FALSE) +
        facet_wrap(~cluster, scales = "free", ncol = n_col) +
        theme_bw() +
        theme(strip.background = element_blank()) +
        labs(title = paste0("adjusted p-value single cell vs pseudobulk - ", file_id),
             subtitle = paste0("r = Pearson, rho = Spearman, p-values floored at ", padj_floor),
             x = "single cell -log10(p_val_adj)", y = "pseudobulk -log10(padj)")
      ggsave(paste0(folder_plot, "06_scatter_padj_singleCell_vs_pseudobulk_", file_id, ".pdf"),
             plot = p_padj, width = 4 * n_col + 1, height = 3.5 * n_row + 1, limitsize = FALSE)
    })
}

message("\nAll done.")
