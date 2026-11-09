# AIM ---------------------------------------------------------------------
# quick check that the sex in the metadata is assigned correctly.
# infer the sex of each sample from XIST/TSIX (female) and the Y-linked genes (male) and compare it with the sex in the metadata.
# the counts are aggregated per sample (sum of the raw counts, Seurat::AggregateExpression) on the final consensus object generated in 03_apply_filter_annotate_fullObject.R.
# NOTICE: the aggregation is run on ALL the genes, so the counts can be normalized by the real library size of each sample (log2 CPM)
# readapted from the sex assessment of the bulk RNAseq project (vst + heatmap of the sex genes)

# LIBRARIES ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# PARAMETERS ------------------------------------------------------------------
options(stringsAsFactors = FALSE)
options(Seurat.object.assay.version = "v5")

# genes used to infer the sex
sex_genes <- c(XIST = "female", TSIX = "female",
               DDX3Y = "male", RPS4Y1 = "male", USP9Y = "male", UTY = "male", EIF1AY = "male", ZFY = "male")

# sample level variable to aggregate on and the sex in the metadata
id_sample <- "sample_id"
id_sex <- "sex"

folder_table <- "../../out/table/analysis_R45_pixi/00_check_sex/"
folder_plot <- "../../out/plot/analysis_R45_pixi/00_check_sex/"

# read in the data --------------------------------------------------------
sobj <- readRDS("../../out/object/analysis_R45_pixi/03_sobj_integrated_cleanup_manualAnnotation_subclusterFiltered_v2.rds")
DefaultAssay(sobj) <- "RNA"

# sample level metadata. sample_id maps 1:1 to the clinical covariates
LUT_sample <- sobj@meta.data %>%
  distinct(.data[[id_sample]], .data[[id_sex]], location, diagnosis_short, pathological_stage_short) %>%
  set_names(c("sample_id", "sex_meta", "location", "diagnosis_short", "pathological_stage_short"))

stopifnot(nrow(LUT_sample) == n_distinct(sobj@meta.data[[id_sample]]))
table(LUT_sample$sex_meta, useNA = "ifany")

# aggregate the counts per sample -----------------------------------------
# sum of the raw counts per sample. NOTICE: AggregateExpression replaces the "_" with "-" in the group names
pbulk <- AggregateExpression(sobj, assays = "RNA", group.by = id_sample, return.seurat = FALSE)$RNA

rm(sobj)
gc()

# match the aggregated names with the sample_id
LUT_sample <- LUT_sample %>%
  mutate(sample_agg = str_replace_all(sample_id, "_", "-"))

# log2 CPM on the library size of all the genes
mat_logcpm <- log2(t(t(as.matrix(pbulk)) / colSums(pbulk) * 1e6) + 1)

# sex genes ---------------------------------------------------------------
genes_found <- names(sex_genes)[names(sex_genes) %in% rownames(mat_logcpm)]
message("sex genes found: ", paste(genes_found, collapse = ", "))
message("sex genes missing: ", paste(setdiff(names(sex_genes), genes_found), collapse = ", "))

mat_sex <- mat_logcpm[genes_found, , drop = FALSE]
sex_found <- sex_genes[genes_found]

# z-score of the log2 CPM per gene across the samples: each gene is centered on 0 and weighs the same, independently of its absolute level
mat_sex_z <- t(scale(t(mat_sex)))
mat_sex_z[is.na(mat_sex_z)] <- 0 # genes with no variation

# score per sample: mean z-score of the Y genes minus the mean z-score of the female genes (> 0 male-like, < 0 female-like)
df_sex <- tibble(sample_agg = colnames(mat_sex_z),
                 female_mean = colMeans(mat_sex_z[sex_found == "female", , drop = FALSE]),
                 male_mean = colMeans(mat_sex_z[sex_found == "male", , drop = FALSE]),
                 XIST = mat_sex["XIST", ]) %>%
  mutate(sex_score = male_mean - female_mean,
         sex_inferred = ifelse(sex_score > 0, "male", "female")) %>%
  left_join(LUT_sample, by = "sample_agg") %>%
  # harmonize the label of the metadata (F/M/female/male...)
  mutate(sex_meta_std = case_when(str_detect(str_to_lower(sex_meta), "^f") ~ "female",
                                  str_detect(str_to_lower(sex_meta), "^m") ~ "male",
                                  TRUE ~ NA_character_),
         match = sex_meta_std == sex_inferred) %>%
  arrange(match, sex_score)

df_sex %>%
  print(n = Inf)

# samples where the metadata and the expression disagree
df_sex %>%
  filter(!match | is.na(match)) %>%
  print(n = Inf)

table(metadata = df_sex$sex_meta_std, inferred = df_sex$sex_inferred, useNA = "ifany")

# plot --------------------------------------------------------------------
# score per sample
p_score <- df_sex %>%
  mutate(sample_id = fct_reorder(sample_id, sex_score)) %>%
  ggplot(aes(x = sample_id, y = sex_score, fill = sex_meta_std)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c(female = "firebrick", male = "steelblue"), na.value = "grey70") +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1, angle = 90)) +
  labs(title = "sex score per sample", subtitle = "mean z-score Y genes - mean z-score XIST/TSIX (> 0 male-like)",
       fill = "sex metadata", x = NULL, y = "sex score")
ggsave(paste0(folder_plot, "00_barplot_sex_score.pdf"), plot = p_score, width = 10, height = 5)

# heatmap: z-score of the log2 CPM per gene across the samples (mat_sex_z)
df_col <- tibble(sample_agg = colnames(mat_sex_z)) %>%
  left_join(df_sex, by = "sample_agg")
colnames(mat_sex_z) <- df_col$sample_id

column_ha <- HeatmapAnnotation(sex_metadata = df_col$sex_meta_std,
                               sex_inferred = df_col$sex_inferred,
                               location = df_col$location,
                               diagnosis = df_col$diagnosis_short,
                               col = list(sex_metadata = c(female = "firebrick", male = "steelblue"),
                                          sex_inferred = c(female = "firebrick", male = "steelblue")),
                               na_col = "grey70")

pdf(paste0(folder_plot, "00_heatmap_sex_genes.pdf"), width = 10, height = 6)
draw(Heatmap(mat_sex_z, name = "z-score\nlog2CPM",
             col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
             top_annotation = column_ha,
             row_split = factor(sex_found, levels = c("female", "male")),
             cluster_row_slices = FALSE,
             column_title = "sex genes (XIST/TSIX: female; Y-linked: male)"),
     heatmap_legend_side = "left", annotation_legend_side = "left")
dev.off()

message("\nAll done.")
