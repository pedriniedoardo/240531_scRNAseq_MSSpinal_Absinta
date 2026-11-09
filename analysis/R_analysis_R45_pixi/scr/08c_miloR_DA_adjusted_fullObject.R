# AIM ---------------------------------------------------------------------
# Differential abundance MS vs CTRL with miloR using the same design of the pseudobulk DE (06_DE_pbulk_deseq2_fullObject.R): ~ location + sex + diagnosis_short, and comparison with the unadjusted design ~ diagnosis_short of 08_miloR_DA_fullObject.R.
# 1. load the milo object saved by 08_miloR_DA_fullObject.R: the neighbourhoods and the counts per sample are already in the object, so ONLY the test is rerun (the neighbourhoods are the same, the results of the two designs are directly comparable neighbourhood by neighbourhood)
# 2. test with the adjusted design (testNhoods, spatial FDR graph-overlap), annotate the neighbourhoods with cell_id and cell_id_subcluster3 (same "Mixed" rule of 08), save the tables and the plots
# 3. compare with the unadjusted results of 08: logFC and SpatialFDR neighbourhood by neighbourhood, and the number of significant neighbourhoods per label
# NOTICE:
# - in testNhoods the tested coefficient is the LAST term of the model, so diagnosis_short has to be the last term of the formula (checked below)
# - the results of 08 are NOT overwritten: all the outputs of this script are in the 08_miloR_adjusted folders and have the prefix 08c
# - the milo object is not saved again (1.1 GB), the only new object is the adjusted DA table

# libraries ---------------------------------------------------------------
library(tidyverse)
library(patchwork)
library(miloR)
library(SingleCellExperiment)
library(ggbeeswarm)

# parameters --------------------------------------------------------------
# comparison: the first level is the reference, positive logFC = enriched in the second level
var_contrast <- "diagnosis_short"
levels_contrast <- c("CTRL", "MS")

# covariates of the model, as in the pseudobulk DE. The variable of interest goes LAST
var_covariate <- c("location", "sex")
levels_location <- c("cervical", "thoracic", "lumbar")

# the sample identifier (unit of replication of the DA test)
id_sample <- "sample_id"

# annotation levels
list_level <- c("cell_id", "cell_id_subcluster3")

# significance threshold on the spatial FDR and minimum fraction of cells with the same label, as in 08_miloR_DA_fullObject.R
fdr_thr <- 0.05
frac_thr <- 0.7

contrast_label <- paste0(levels_contrast[2], "_vs_", levels_contrast[1])
formula_adj <- as.formula(paste("~", paste(c(var_covariate, var_contrast), collapse = " + ")))
formula_adj

# results of 08 (unadjusted design)
folder_table_ref <- "../../out/table/analysis_R45_pixi/08_miloR/"
folder_object_ref <- "../../out/object/analysis_R45_pixi/08_miloR/"

folder_table <- "../../out/table/analysis_R45_pixi/08_miloR_adjusted/"
folder_plot <- "../../out/plot/analysis_R45_pixi/08_miloR_adjusted/"
dir.create(folder_table, recursive = TRUE, showWarnings = FALSE)
dir.create(folder_plot, recursive = TRUE, showWarnings = FALSE)

# read in the data --------------------------------------------------------
milo <- readRDS(paste0(folder_object_ref, "08_miloR_object_", contrast_label, ".rds"))

# the neighbourhoods and the counts are in the object
milo
dim(nhoodCounts(milo))
stopifnot(ncol(nhoods(milo)) > 0, nrow(nhoodCounts(milo)) == ncol(nhoods(milo)))

# design ------------------------------------------------------------------
# one row per sample, in the order of the columns of the count matrix. sample_id maps 1:1 to the covariates
design_milo <- data.frame(colData(milo))[, c(id_sample, var_covariate, var_contrast)] %>%
  distinct()

# every sample has to have a single value for each variable
stopifnot(nrow(design_milo) == n_distinct(design_milo[[id_sample]]))

design_milo <- design_milo %>%
  mutate(location = factor(location, levels = levels_location),
         sex = factor(sex),
         !!var_contrast := factor(.data[[var_contrast]], levels = levels_contrast))
rownames(design_milo) <- design_milo[[id_sample]]
design_milo <- design_milo[colnames(nhoodCounts(milo)), , drop = FALSE]

# samples per level of the covariates and of the contrast
table(design_milo[[var_contrast]], design_milo$location)
table(design_milo[[var_contrast]], design_milo$sex)

# the tested coefficient is the last column of the model matrix and has to be the contrast
mm <- model.matrix(formula_adj, data = design_milo)
colnames(mm)
stopifnot(tail(colnames(mm), 1) == paste0(var_contrast, levels_contrast[2]),
          qr(mm)$rank == ncol(mm))

# differential abundance testing ------------------------------------------
# same settings of 08: NB GLM (edgeR) per neighbourhood, spatial FDR with graph-overlap weights. Only the formula changes
da_adj <- testNhoods(milo, design = formula_adj, design.df = design_milo, fdr.weighting = "graph-overlap")

da_adj %>%
  arrange(SpatialFDR) %>%
  head()

da_adj %>%
  write_tsv(paste0(folder_table, "08c_miloR_DA_adjusted_", contrast_label, ".tsv"))

# diagnostic plots: uncorrected p values and volcano
p_pval <- ggplot(da_adj, aes(PValue)) +
  geom_histogram(bins = 50) +
  theme_bw() +
  labs(title = paste(deparse(formula_adj), collapse = ""))

p_volcano <- ggplot(da_adj, aes(logFC, -log10(SpatialFDR))) +
  geom_point(shape = 1, alpha = 0.5) +
  geom_hline(yintercept = -log10(fdr_thr), linetype = "dashed", col = "gray") +
  theme_bw() +
  theme(strip.background = element_blank())

ggsave(paste0(folder_plot, "08c_miloR_adjusted_diagnostic_pvalue_volcano.pdf"), plot = p_pval | p_volcano, width = 12, height = 5)

# comparison with the unadjusted design of 08 -----------------------------
da_ref <- read_tsv(paste0(folder_table_ref, "08_miloR_DA_", contrast_label, ".tsv"), show_col_types = FALSE)

# same neighbourhoods in the two runs (same object), matched by Nhood
stopifnot(nrow(da_ref) == nrow(da_adj), all(da_ref$Nhood %in% da_adj$Nhood))

df_cmp <- inner_join(da_ref %>% dplyr::select(Nhood, logFC, PValue, SpatialFDR),
                     da_adj %>% dplyr::select(Nhood, logFC, PValue, SpatialFDR),
                     by = "Nhood", suffix = c("_ref", "_adj"))

# agreement of the logFC and number of significant neighbourhoods in the two designs
cor(df_cmp$logFC_ref, df_cmp$logFC_adj)
cor(df_cmp$logFC_ref, df_cmp$logFC_adj, method = "spearman")
table(ref = df_cmp$SpatialFDR_ref < fdr_thr, adjusted = df_cmp$SpatialFDR_adj < fdr_thr)

df_cmp %>%
  mutate(sig = case_when(SpatialFDR_ref < fdr_thr & SpatialFDR_adj < fdr_thr ~ "both",
                         SpatialFDR_ref < fdr_thr ~ "unadjusted only",
                         SpatialFDR_adj < fdr_thr ~ "adjusted only",
                         TRUE ~ "none")) %>%
  dplyr::count(sig) %>%
  write_tsv(paste0(folder_table, "08c_miloR_comparison_nSignificant_", contrast_label, ".tsv"))

p_logfc <- df_cmp %>%
  mutate(sig = case_when(SpatialFDR_ref < fdr_thr & SpatialFDR_adj < fdr_thr ~ "both",
                         SpatialFDR_ref < fdr_thr ~ "unadjusted only",
                         SpatialFDR_adj < fdr_thr ~ "adjusted only",
                         TRUE ~ "none")) %>%
  arrange(sig != "none") %>%
  ggplot(aes(x = logFC_ref, y = logFC_adj, color = sig)) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", col = "gray40") +
  scale_color_manual(values = c("none" = "grey70", "both" = "red", "unadjusted only" = "blue", "adjusted only" = "orange")) +
  theme_bw() +
  labs(x = "logFC ~ diagnosis", y = paste0("logFC ", paste(deparse(formula_adj), collapse = "")), color = paste0("SpatialFDR < ", fdr_thr))

p_fdr <- df_cmp %>%
  ggplot(aes(x = -log10(SpatialFDR_ref), y = -log10(SpatialFDR_adj))) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", col = "gray40") +
  geom_hline(yintercept = -log10(fdr_thr), linetype = "dashed", col = "red") +
  geom_vline(xintercept = -log10(fdr_thr), linetype = "dashed", col = "red") +
  theme_bw() +
  labs(x = "-log10 SpatialFDR ~ diagnosis", y = "-log10 SpatialFDR adjusted")

ggsave(paste0(folder_plot, "08c_miloR_comparison_unadjusted_vs_adjusted.pdf"), plot = p_logfc | p_fdr, width = 14, height = 6)

# neighbourhood graph on the UMAP, adjusted results -----------------------
# the neighbourhood graph (buildNhoodGraph) is already in the saved object. only the significant neighbourhoods are colored (left) and all the neighbourhoods (right, the SpatialFDR is set to 0 in a copy of the table only for this plot, plotNhoodGraphDA sets the logFC to 0 above alpha)
scale_logfc <- scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, na.value = "grey50", name = paste0("logFC\n", levels_contrast[2], " vs ", levels_contrast[1]))

p_nhood <- plotNhoodGraphDA(milo, da_adj, layout = "UMAP", alpha = fdr_thr) +
  scale_logfc +
  ggtitle(paste0("Adjusted design, significant neighbourhoods only (SpatialFDR < ", fdr_thr, ")"))

p_nhood_cont <- plotNhoodGraphDA(milo, da_adj %>% mutate(SpatialFDR = 0), layout = "UMAP", alpha = fdr_thr) +
  scale_logfc +
  ggtitle("Adjusted design, all neighbourhoods (logFC, no significance threshold)")

ggsave(paste0(folder_plot, "08c_miloR_adjusted_nhoodGraphDA_UMAP.pdf"), plot = p_nhood | p_nhood_cont, width = 20, height = 7)

# annotation of the neighbourhoods ----------------------------------------
# same steps of 08 for each level of the annotation (annotateNhoods, "Mixed" below frac_thr, summary of the significant neighbourhoods, beeswarm), plus the comparison with the unadjusted design per label
# level <- "cell_id"
list_da_annotated <- lapply(list_level, function(level){
  da_level <- annotateNhoods(milo, da_adj, coldata_col = level)

  # distribution of the fraction of cells with the label in each neighbourhood
  p_frac <- ggplot(da_level, aes(.data[[paste0(level, "_fraction")]])) +
    geom_histogram(bins = 50) +
    theme_bw() +
    labs(x = paste0(level, "_fraction"))

  # neighbourhoods with a mix of labels
  da_level <- da_level %>%
    mutate(!!level := ifelse(.data[[paste0(level, "_fraction")]] < frac_thr, "Mixed", .data[[level]]))

  da_level %>%
    write_tsv(paste0(folder_table, "08c_miloR_DA_adjusted_", contrast_label, "_", level, ".tsv"))

  # number of significant neighbourhoods per label
  da_level %>%
    group_by(.data[[level]]) %>%
    summarise(n_nhood = n(),
              n_up = sum(SpatialFDR < fdr_thr & logFC > 0),
              n_down = sum(SpatialFDR < fdr_thr & logFC < 0),
              .groups = "drop") %>%
    write_tsv(paste0(folder_table, "08c_miloR_DA_summary_adjusted_", contrast_label, "_", level, ".tsv"))

  # comparison per label: the labels are the same in the two designs (they depend only on the neighbourhoods), matched by Nhood
  da_level %>%
    dplyr::select(Nhood, label = all_of(level), logFC_adj = logFC, SpatialFDR_adj = SpatialFDR) %>%
    inner_join(df_cmp %>% dplyr::select(Nhood, logFC_ref, SpatialFDR_ref), by = "Nhood") %>%
    group_by(label) %>%
    summarise(n_nhood = n(),
              median_logFC_ref = median(logFC_ref),
              median_logFC_adj = median(logFC_adj),
              n_sig_ref = sum(SpatialFDR_ref < fdr_thr),
              n_sig_adj = sum(SpatialFDR_adj < fdr_thr),
              best_SpatialFDR_ref = min(SpatialFDR_ref),
              best_SpatialFDR_adj = min(SpatialFDR_adj),
              .groups = "drop") %>%
    arrange(dplyr::desc(median_logFC_adj)) %>%
    write_tsv(paste0(folder_table, "08c_miloR_comparison_perLabel_", contrast_label, "_", level, ".tsv"))

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
    labs(title = paste0("miloR DA ", levels_contrast[2], " vs ", levels_contrast[1], " - ", level, " - adjusted for ", paste(var_covariate, collapse = " + ")),
         subtitle = paste0("SpatialFDR < ", fdr_thr, ", positive logFC = enriched in ", levels_contrast[2]),
         y = NULL)

  # the height follows the number of labels
  ggsave(paste0(folder_plot, "08c_miloR_adjusted_beeswarm_", level, ".pdf"), plot = p_bee | p_frac,
         width = 14, height = max(6, 0.25 * length(order_level) + 2), limitsize = FALSE)

  da_level
})
names(list_da_annotated) <- list_level

message("\nAll done.")
