# AIM ---------------------------------------------------------------------
# the EBV positive raw barcodes (spinal 07 and brain WM/CX 11) are mostly unlabelled and almost empty droplets.
# plot the expression of a panel of classical marker genes for each individual barcode (one row per barcode) to see if any of them looks like a specific cell type.
# there is no real grouping and the proportion of cells expressing a gene is meaningless for a single barcode, therefore a heatmap instead of the dotplot.

# libraries -----------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(cowplot)

# parameters ------------------------------------------------------------------
dir_table <- "../../out/table/analysis_R45_pixi"
dir_plot <- "../../out/plot/analysis_R45_pixi"
dir_pipeline <- "/beegfs/scratch/ric.absinta/ric.absinta/pipeline"

# cellranger output of the two datasets (same runs as 07 and 11)
list_dataset <- list(
  "Spinal cord snRNAseq" = list(
    hits = file.path(dir_table, "07_quick_check_EBV_test_rawBarcodes_viralPositiveBarcodes_raw.tsv"),
    dir_cellranger = file.path(dir_pipeline, "scrnaseq_cosr_standard_workflow_AbsintaSpinal_snRNAseq_testEBV/results/cellranger/merged")),
  "Brain WM/CX snRNAseq" = list(
    hits = file.path(dir_table, "11_quick_check_EBV_test_rawBarcodes_AbsintaBrainWMCX_viralPositiveBarcodes_raw.tsv"),
    dir_cellranger = file.path(dir_pipeline, "scrnaseq_cosr_standard_workflow_brainFull_testEBV/results/cellranger/merged"))
)

# marker panel (from 27_subcluster_analysis_VAS_plot.R)
shortlist_features_list_long <- list(
  IMMUNE = c("LYVE1","CD163","MRC1","LINGO1","HSPA1A","MOBP","CD22","CD83","HIF1A","VEGFA","SOD1","TREM2","CX3CR1","P2RY12","C3","CSF1R", "CD74", "RUNX1","C1QB","PTPRC","AIF1","HLA-DRA","TYROBP"),
  B_CELLS = c("IGHG1", "CD38"),
  T_CELLS =  c("SKAP1", "CD8A", "CD2"),
  OLIGOLINEAGE = c("PLP1","MOG","PPP1R16B","TNS3","HMGB1","CD81","B2M","C1QL1","HLA-A","HLA-C","NLGN4X","OLIG1","OLIG2","CSPG4", "PDGFRA", "SOX6", "SOX10", "MBP","MAG"),
  ASTRO = c("AQP4", "GFAP", "CD44", "AQP1", "VIM","APOE", "VCAN", "STAT3", "ABCA1", "TNC", "SDC4","SLC1A2","S100B"),
  NEURONS = c("GAD2", "PVALB", "SV2C", "VIP", "TLE4", "CUX2", "THY1", "SLC17A7", "NRGN", "SATB2", "RORB", "SST", "STX1A", "STX1B", "SYP", "TH", "NEFL","SYT1"),
  ENDO = c("VWF","CDH5","TEK","PECAM1","FLT1","KDR","NOS3","MCAM","MMRN1","CLDN5","BMX","ANGPT2","GJA4","TIE1","ROBO4","ECSCR"),
  PERICYTE = c("PDGFRB","DES","ACTA2","ANPEP","RGS5","ABCC9","KCNJ8","CD248","DLK1","NT5E","ANGPT1"),
  SCHWANN = c("PMP22","MPZ","PRX"),
  EPENDYMA = c("CFAP299","DNAH7","DNAH9"),
  STROMAL = c("LAMA2","RBMS3","CEMIP","GPC6")
)

# make the list as a data.frame()
shortlist_features_df <- lapply(shortlist_features_list_long, function(x) {
  data.frame(gene = x)
}) %>%
  bind_rows(.id = "cat")

GOI <- unique(shortlist_features_df$gene)

# expression of the positive barcodes -----------------------------------------
# one row per positive barcode (a barcode with hits on several viral genes is kept once)
# the raw matrix has ~1M barcodes per sample: load it, keep only the positive barcodes, then log normalise on the full-transcriptome library size of the barcode (the same normalisation of NormalizeData(), scale factor 1e4)
# sample_id <- "Sample_2015_017_300"
get_expression <- function(sample_id, barcodes, dir_cellranger) {
  raw_dir <- file.path(dir_cellranger, sample_id, "outs", "raw_feature_bc_matrix")
  message("sample: ", sample_id, "  ", length(barcodes), " barcode(s)")
  raw_mat <- Read10X(data.dir = raw_dir)
  mat <- raw_mat[, barcodes, drop = FALSE]
  rm(raw_mat); gc()

  # library size over the full transcriptome, before any gene subsetting
  lib_size <- Matrix::colSums(mat)

  gene_present <- intersect(GOI, rownames(mat))
  mat_goi <- mat[gene_present, , drop = FALSE]

  # log1p(count / library size * 1e4), sparse -> dense is fine here (tens of barcodes)
  norm_goi <- log1p(sweep(as.matrix(mat_goi), 2, lib_size, "/") * 1e4)

  as.data.frame(norm_goi) %>%
    rownames_to_column("gene") %>%
    pivot_longer(names_to = "barcode", values_to = "norm_exp", -gene) %>%
    left_join(as.data.frame(as.matrix(mat_goi)) %>%
                rownames_to_column("gene") %>%
                pivot_longer(names_to = "barcode", values_to = "count", -gene),
              by = c("gene", "barcode")) %>%
    mutate(sample_id = sample_id, lib_size = lib_size[barcode])
}

# x <- list_dataset$`Brain WM/CX snRNAseq`
# dataset <- "Brain WM/CX snRNAseq"
df_expression <- imap_dfr(list_dataset, function(x, dataset) {
  df_hits <- read_tsv(x$hits, show_col_types = FALSE)

  # one row per positive barcode, the viral UMI is summed over the viral genes
  df_barcode <- df_hits %>%
    group_by(sample_id, barcode) %>%
    summarise(viral_UMI = sum(count),
              is_called_cell = any(is_called_cell),
              .groups = "drop")

  df_barcode %>%
    group_by(sample_id) %>%
    group_map(function(d, key) {
      get_expression(key$sample_id, d$barcode, x$dir_cellranger)
    }, .keep = FALSE) %>%
    bind_rows() %>%
    left_join(df_barcode, by = c("sample_id", "barcode")) %>%
    mutate(dataset = dataset)
})

# genes of the panel not found in the matrix
message("genes of the panel missing from the matrix: ", paste(setdiff(GOI, df_expression$gene), collapse = ", "))

write_tsv(df_expression, file.path(dir_table, "18_EBV_positive_barcodes_marker_expression.tsv"))

# heatmap -----------------------------------------------------------------
# the two datasets are kept separated: one plot per dataset, the scaling is computed within each dataset
# two scalings are tested:
# - gene: z-score of each gene across the barcodes of the dataset. highlights the barcodes where each gene is mostly expressed
# - barcode: z-score of each barcode across the genes of the panel. highlights the genes mostly expressed within each barcode
# row label: sample, barcode, library size (UMI of the full transcriptome) and a star if the barcode is in the cellranger called cells
# the rows are ordered by library size within sample

# scale on the unique barcode x gene values, before the gene is replicated across categories
# a constant vector (e.g. a gene never detected) gives NaN with scale(), set it to 0
df_scaled <- df_expression %>%
  mutate(row_id = paste(sample_id, barcode)) %>%
  group_by(dataset, gene) %>%
  mutate(norm_exp.scaled_gene = replace_na(scale(norm_exp)[, 1], 0)) %>%
  group_by(dataset, row_id) %>%
  mutate(norm_exp.scaled_barcode = replace_na(scale(norm_exp)[, 1], 0)) %>%
  ungroup()

# the largest barcode on top (levels are reversed because the first level is at the bottom of the y axis)
df_plot <- df_scaled %>%
  mutate(label = paste0(sample_id, " ", barcode, " (", lib_size, " UMI)", if_else(is_called_cell, " *", ""))) %>%
  left_join(shortlist_features_df, by = "gene", relationship = "many-to-many") %>%
  mutate(cat = factor(cat, levels = names(shortlist_features_list_long)),
         gene = factor(gene, levels = GOI),
         dataset = factor(dataset, levels = names(list_dataset))) %>%
  arrange(sample_id, desc(lib_size), barcode) %>%
  mutate(label = factor(label, levels = rev(unique(label))))

df_spinal <- df_plot %>% filter(dataset == "Spinal cord snRNAseq")
df_brain <- df_plot %>% filter(dataset == "Brain WM/CX snRNAseq")

# spinal cord, scaled within genes
p_spinal_gene <- df_spinal %>%
  ggplot(aes(x = gene, y = label)) +
  geom_tile(aes(fill = norm_exp.scaled_gene), width = 0.95, height = 0.95, linewidth = 0.2, color = "white") +
  facet_grid(~ cat, scales = "free", space = "free") +
  theme_cowplot() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90, size = 7),
        axis.text.y = element_text(size = 7),
        strip.text.x = element_text(angle = 90, hjust = 0),
        strip.text.y = element_text(angle = 0, hjust = 0)) +
  scale_fill_gradientn(colours = viridis::turbo(20), limits = c(-0.5, 3), oob = scales::squish, name = "exp scaled\nwithin genes") +
  labs(x = NULL, y = NULL, title = "Spinal cord snRNAseq",
       caption = "one row per EBV positive raw barcode. log1p(count / barcode UMI * 1e4), z-score gene-wise. * = barcode in the cellranger called cells")
p_spinal_gene
ggsave(file.path(dir_plot, "18_EBV_positive_barcodes_marker_heatmap_spinal_scale_gene.pdf"), p_spinal_gene,
       width = 32, height = max(6, 2 + n_distinct(df_spinal$row_id) * 0.18), limitsize = FALSE)

# spinal cord, scaled within barcodes
p_spinal_barcode <- df_spinal %>%
  ggplot(aes(x = gene, y = label)) +
  geom_tile(aes(fill = norm_exp.scaled_barcode), width = 0.95, height = 0.95, linewidth = 0.2, color = "white") +
  facet_grid(~ cat, scales = "free", space = "free") +
  theme_cowplot() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90, size = 7),
        axis.text.y = element_text(size = 7),
        strip.text.x = element_text(angle = 90, hjust = 0),
        strip.text.y = element_text(angle = 0, hjust = 0)) +
  scale_fill_gradientn(colours = viridis::turbo(20), limits = c(-0.5, 3), oob = scales::squish, name = "exp scaled\nwithin barcodes") +
  labs(x = NULL, y = NULL, title = "Spinal cord snRNAseq",
       caption = "one row per EBV positive raw barcode. log1p(count / barcode UMI * 1e4), z-score barcode-wise. * = barcode in the cellranger called cells")
p_spinal_barcode
ggsave(file.path(dir_plot, "18_EBV_positive_barcodes_marker_heatmap_spinal_scale_barcode.pdf"), p_spinal_barcode,
       width = 32, height = max(6, 2 + n_distinct(df_spinal$row_id) * 0.18), limitsize = FALSE)

# brain WM/CX, scaled within genes
p_brain_gene <- df_brain %>%
  ggplot(aes(x = gene, y = label)) +
  geom_tile(aes(fill = norm_exp.scaled_gene), width = 0.95, height = 0.95, linewidth = 0.2, color = "white") +
  facet_grid(~ cat, scales = "free", space = "free") +
  theme_cowplot() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90, size = 7),
        axis.text.y = element_text(size = 7),
        strip.text.x = element_text(angle = 90, hjust = 0),
        strip.text.y = element_text(angle = 0, hjust = 0)) +
  scale_fill_gradientn(colours = viridis::turbo(20), limits = c(-0.5, 3), oob = scales::squish, name = "exp scaled\nwithin genes") +
  labs(x = NULL, y = NULL, title = "Brain WM/CX snRNAseq",
       caption = "one row per EBV positive raw barcode. log1p(count / barcode UMI * 1e4), z-score gene-wise. * = barcode in the cellranger called cells")
p_brain_gene
ggsave(file.path(dir_plot, "18_EBV_positive_barcodes_marker_heatmap_brain_scale_gene.pdf"), p_brain_gene,
       width = 32, height = max(6, 2 + n_distinct(df_brain$row_id) * 0.18), limitsize = FALSE)

# brain WM/CX, scaled within barcodes
p_brain_barcode <- df_brain %>%
  ggplot(aes(x = gene, y = label)) +
  geom_tile(aes(fill = norm_exp.scaled_barcode), width = 0.95, height = 0.95, linewidth = 0.2, color = "white") +
  facet_grid( ~ cat, scales = "free", space = "free") +
  theme_cowplot() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90, size = 7),
        axis.text.y = element_text(size = 7),
        strip.text.x = element_text(angle = 90, hjust = 0),
        strip.text.y = element_text(angle = 0, hjust = 0)) +
  scale_fill_gradientn(colours = viridis::turbo(20), limits = c(-0.5, 3), oob = scales::squish, name = "exp scaled\nwithin barcodes") +
  labs(x = NULL, y = NULL, title = "Brain WM/CX snRNAseq",
       caption = "one row per EBV positive raw barcode. log1p(count / barcode UMI * 1e4), z-score barcode-wise. * = barcode in the cellranger called cells")
p_brain_barcode
ggsave(file.path(dir_plot, "18_EBV_positive_barcodes_marker_heatmap_brain_scale_barcode.pdf"), p_brain_barcode,
       width = 32, height = max(6, 2 + n_distinct(df_brain$row_id) * 0.18), limitsize = FALSE)
