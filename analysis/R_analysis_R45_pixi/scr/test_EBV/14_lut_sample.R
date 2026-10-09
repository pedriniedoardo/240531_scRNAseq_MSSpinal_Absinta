# AIM ---------------------------------------------------------------------
# complete the summary tables from EBV test with the sample metadata.
# if availabe also add the barcode metadata

# libraries ---------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(patchwork)

# read in the data --------------------------------------------------------
dir_table <- "../../out/table/analysis_R45_pixi"
dir_project <- "/beegfs/scratch/ric.cosr/pedrini.edoardo/project_edoardo"

# per sample EBV summaries
# load only the summcires that have at least one EBV read, therefore exlcude jakel and shirmer
# summary_jakel <- read_tsv(file.path(dir_table, "08_quick_check_EBV_test_rawBarcodes_Jakel_perSampleSummary.tsv"), show_col_types = FALSE)
# summary_shirmer <- read_tsv(file.path(dir_table, "10_quick_check_EBV_test_rawBarcodes_Shirmer2019_perSampleSummary.tsv"), show_col_types = FALSE)
summary_spinal <- read_tsv(file.path(dir_table, "07_quick_check_EBV_test_rawBarcodes_perSampleSummary.tsv"), show_col_types = FALSE)
summary_brain <- read_tsv(file.path(dir_table, "11_quick_check_EBV_test_rawBarcodes_AbsintaBrainWMCX_perSampleSummary.tsv"), show_col_types = FALSE)
summary_visium <- read_tsv(file.path(dir_table, "13_quick_check_EBV_test_rawBarcodes_VisiumBrain_perSampleSummary.tsv"), show_col_types = FALSE)
summary_gse189141 <- read_tsv(file.path(dir_table, "09_quick_check_EBV_test_rawBarcodes_GSE189141_perSampleSummary.tsv"), show_col_types = FALSE)

# sample metadata
lut_spinal_full <- read_csv(file.path(dir_project,"240531_scRNAseq_MSSpinal_Absinta/data/LUT_sample_full.csv"), show_col_types = FALSE)
lut_spinal_demye <- read_csv(file.path(dir_project, "240531_scRNAseq_MSSpinal_Absinta/data/20260416_spinal_demyelination_update.csv"), show_col_types = FALSE)
lut_brain <- read_csv(file.path(dir_project, "220501_scRNAseq_MSbrain_Absinta/data/LUT_sample_WM_CX.csv"), show_col_types = FALSE)
lut_visium <- read_csv(file.path(dir_project, "260216_visium_brain_absinta/data/LUT_samples.csv"), show_col_types = FALSE)
# for summary_gse189141 the metadata is backed into the sample name

# build the spinal sample metadata as in the annotated object
# the sample level columns of the object come from three steps (R44/26_annotation.R and R44/001_update annotation_demyelination.R):
# 1. LUT_sample_full.csv (clinical metadata)
# 2. the *_short factors derived from diagnosis, pathological_stage and demyelination
# 3. the quantified demyelination (20260416 csv) joined by original_sample_name, plus the high/low classes
# the object has the column "demyelination.score" (data.frame() renamed it) and "location" (the csv header has a trailing space)
# strip the trailing spaces from the values (e.g. "no demyelination ") otherwise the factor levels below return NA
lut_spinal_clean <- lut_spinal_full %>%
  rename_with(str_trim) %>%
  rename(demyelination.score = `demyelination score`) %>%
  mutate(across(c(diagnosis, pathological_stage, demyelination), str_trim))

# classes of the quantified demyelination (cutoffs as in the R44 script, controls are "no")
lut_spinal_demye_class <- lut_spinal_demye %>%
  mutate(
    demye_WM_class = case_when(is.na(demye_WM_prop) ~ NA,
                               diagnosis == "Non-demented control" ~ "no",
                               demye_WM_prop > 0.3 ~ "high",
                               TRUE ~ "low") %>%
      factor(levels = c("no", "low", "high")),
    demye_tot_class = case_when(is.na(demye_tot_prop) ~ NA,
                                diagnosis == "Non-demented control" ~ "no",
                                demye_tot_prop > 0.3 ~ "high",
                                TRUE ~ "low") %>%
      factor(levels = c("no", "low", "high")),
    demye_GM_class = case_when(is.na(demye_GM_prop) ~ NA,
                               diagnosis == "Non-demented control" ~ "no",
                               demye_GM_prop > 0.5 ~ "high",
                               TRUE ~ "low") %>%
      factor(levels = c("no", "low", "high"))) %>%
  # drop the columns already in the clinical LUT, join on original_sample_name
  select(-c(sample_id, autopsy, sex, age, location, diagnosis))

lut_spinal <- lut_spinal_clean %>%
  mutate(diagnosis_short = factor(diagnosis, levels = c("Non-demented control", "Multiple sclerosis"), labels = c("CTRL", "MS")),
         pathological_stage_short = factor(pathological_stage, levels = c("control", "inactive", "active demyelination"), labels = c("CTRL", "IN", "ACT")),
         demyelination_short = factor(demyelination, levels = c("no demyelination", "less than half", "more than half"), labels = c("NO", "LESS50", "MORE50"))) %>%
  left_join(lut_spinal_demye_class, by = "original_sample_name")

# join the metadata -------------------------------------------------------
# Spinal (Absinta): sample_id matches the pipeline sample name 1:1
df_summary_spinal <- left_join(summary_spinal, lut_spinal, by = "sample_id")

# Brain WM/CX (Absinta): the official id (e.g. s9) is embedded in the pipeline sample name (GSM...__s9__Homo_sapiens__RNA-Seq)
# the sample_id of the LUT is the original library name, so it is dropped to avoid clashing with the pipeline sample_id
df_summary_brain <- summary_brain %>%
  mutate(official_id = str_extract(sample_id, "__(s\\d+)__",group = 1)) %>%
  left_join(select(lut_brain, -sample_id), by = "official_id")

# Visium brain
df_summary_brain_visium <- left_join(summary_visium, lut_visium, by = c("sample_id" = "library_id"))

# GSE189141: GSM_TX<donor>_<timepoint>_GEX
df_summary_gse189141 <- summary_gse189141 %>%
  mutate(gsm = str_extract(sample_id, "^GSM\\d+"),
         donor = str_extract(sample_id, "TX\\d+"),
         time_point = as.numeric(str_extract(sample_id, "_(\\d+)_",group = 1)))

# save the tables ---------------------------------------------------------
list_tab <- list(Spinal = df_summary_spinal,
     GSE189141 = df_summary_gse189141,
     AbsintaBrainWMCX = df_summary_brain,
     VisiumBrain = df_summary_brain_visium)

# save the tables
iwalk(list_tab, function(x,nm){
  write_tsv(x,file.path(dir_table, paste0("14_quick_check_EBV_test_rawBarcodes_", nm, "_perSampleSummary_metadata.tsv")))
})

# summarise the metadata --------------------------------------------------
df_group_spinal <- df_summary_spinal %>%
  group_by(diagnosis) %>%
  summarise(n_samples = n(),
            tot_barcode = sum(n_raw_barcodes, na.rm = TRUE),
            # tot_UMI = sum(total_UMI, na.rm = TRUE),
            tot_viral_hits = sum(n_viral_hits, na.rm = TRUE),
            tot_called_hits = sum(n_called_cell_hits, na.rm = TRUE))

df_group_brain <- df_summary_brain %>%
  group_by(pathology_class) %>%
  summarise(n_samples = n(),
            tot_barcode = sum(n_raw_barcodes, na.rm = TRUE),
            # tot_UMI = sum(total_UMI, na.rm = TRUE),
            tot_viral_hits = sum(n_viral_hits, na.rm = TRUE),
            tot_called_hits = sum(n_called_cell_hits, na.rm = TRUE))

df_group_visium <- df_summary_brain_visium %>%
  group_by(sample_classification2) %>%
  summarise(n_samples = n(),
            tot_barcode = sum(n_raw_barcodes, na.rm = TRUE),
            # tot_UMI = sum(total_UMI, na.rm = TRUE),
            tot_viral_hits = sum(n_viral_hits, na.rm = TRUE),
            tot_called_hits = sum(n_called_cell_hits, na.rm = TRUE))

df_group_gse189141 <- df_summary_gse189141 %>%
  group_by(time_point) %>%
  summarise(n_samples = n(),
            tot_barcode = sum(n_raw_barcodes, na.rm = TRUE),
            # tot_UMI = sum(total_UMI, na.rm = TRUE),
            tot_viral_hits = sum(n_viral_hits, na.rm = TRUE),
            tot_called_hits = sum(n_called_cell_hits, na.rm = TRUE))

# summary plot per dataset ------------------------------------------------
# one panel per dataset, bars = viral hits per million UMI of the raw matrix (the datasets have very different sequencing depth, the raw barcode number is a poor proxy: ~1M per nuclei sample vs ~5k per visium slide)
# the label is the raw number of hits (of which in called cells/spots)
df_plot <- bind_rows(
  df_group_spinal %>% rename(group = diagnosis) %>% mutate(dataset = "Spinal cord snRNAseq", group = as.character(group)),
  df_group_brain %>% rename(group = pathology_class) %>% mutate(dataset = "Brain WM/CX snRNAseq", group = as.character(group)),
  df_group_visium %>% rename(group = sample_classification2) %>% mutate(dataset = "Brain Visium", group = as.character(group)),
  df_group_gse189141 %>% rename(group = time_point) %>% mutate(dataset = "GSE189141 (time point)", group = as.character(group))
) %>%
  mutate(
    # hits_per_million = tot_viral_hits / tot_UMI * 1e6,
    label = paste0(tot_viral_hits, "\n(", tot_called_hits, ")\nn=", n_samples))

# rescale per dataset
p_summary <- df_plot %>%
  # ggplot(aes(x = group, y = hits_per_million, fill = dataset)) +
  ggplot(aes(x = group, y = tot_viral_hits, fill = dataset)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = label), vjust = -0.2, size = 3) +
  facet_wrap(~dataset, scales = "free",nrow = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(x = NULL, y = "EBV viral hits",
       caption = "label: n hits \n(n in called cells)\n n samples") +
  # labs(x = NULL, y = "EBV hits per million UMI",
  #      caption = "label: n hits \n(n in called cells)\n n samples") +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
p_summary

ggsave(file.path("../../out/plot/analysis_R45_pixi", "14_quick_check_EBV_test_summary_perDataset_scaled.pdf"), p_summary, width = 12, height = 6)

# use the same scale across all the samples
p_summary2 <- df_plot %>%
  # ggplot(aes(x = group, y = hits_per_million, fill = dataset)) +
  ggplot(aes(x = group, y = tot_viral_hits, fill = dataset)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = label), vjust = -0.2, size = 3) +
  facet_wrap(~dataset, scales = "free_x",nrow = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(x = NULL, y = "EBV viral hits",
       caption = "label: n hits \n(n in called cells)\n n samples") +
  # labs(x = NULL, y = "EBV hits per million UMI",
  #      caption = "label: n hits \n(n in called cells)\n n samples") +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
p_summary2
ggsave(file.path("../../out/plot/analysis_R45_pixi", "14_quick_check_EBV_test_summary_perDataset_unscaled.pdf"), p_summary2, width = 12, height = 6)

# split of the viral counts per gene --------------------------------------
# EBV/viral marker panel -- identical (all 86 annotated genes) to scr/test_EBV/07_quick_check_EBV_test_Spinal.R
list_markers_viral <- list(
  EBV_LATENT = c("EBNA-1", "EBNA-2", "EBNA-3A", "EBNA-3B-EBNA-3C", "EBNA-LP", "LMP-1", "RPMS1"),
  EBV_LYTIC_IE = c("BZLF1", "BRLF1"),
  EBV_LYTIC_E = c("BMRF1", "BALF2", "BSLF2-BMLF1", "BALF5", "BHRF1"),
  EBV_LYTIC_L = c("BLLF1", "BcLF1", "BLRF2", "BALF4"),
  EBV_OTHER = c(
    "BNLF2b", "BNLF2a", "BARF1", "BALF1", "BARF0", "BALF3", "A73", "LF2", "LF1", "BILF1", "LF3",
    "BILF2", "BdRF1", "BVRF2", "BVLF1", "BVRF1", "BXRF1", "BXLF2", "BXLF1", "BTRF1", "BcRF1",
    "BDLF3", "BDLF2", "BDLF1", "BGLF2", "BGLF1", "BDLF4", "BDLF3.5", "BGRF1-BDRF1", "BBLF1",
    "BGLF5", "BGLF4", "BGLF3.5", "BGLF3", "BBRF3", "BBLF2-BBLF3", "BBRF2", "BBRF1", "BBLF4",
    "BKRF4", "BKRF3", "BKRF2", "BRRF2", "BRRF1", "BZLF2", "BLLF2", "BLRF1", "BLLF3", "BSRF1",
    "BSLF1", "BMRF2", "BaRF1", "BORF2", "BORF1", "BPLF1", "BOLF1", "BFRF3", "BFRF2", "BFRF1",
    "BFRF1A", "BFLF2", "BFLF1", "BHLF1", "BWRF1", "BCRF1",
    "rna-NC-007605.1:6956..7128", "rna-NC-007605.1:6629..6795", "BNRF1"
  )
)

df_gene_class <- stack(list_markers_viral) %>%
  dplyr::rename(gene = values, viral_class = ind)

# the viral positive barcodes (one row per barcode x gene) of each dataset
df_hits_spinal <- read_tsv(file.path(dir_table, "07_quick_check_EBV_test_rawBarcodes_viralPositiveBarcodes_raw.tsv"))
df_hits_brain <- read_tsv(file.path(dir_table, "11_quick_check_EBV_test_rawBarcodes_AbsintaBrainWMCX_viralPositiveBarcodes_raw.tsv"))
df_hits_visium <- read_tsv(file.path(dir_table, "13_quick_check_EBV_test_rawBarcodes_VisiumBrain_viralPositiveBarcodes_raw.tsv"))
df_hits_gse189141 <- read_tsv(file.path(dir_table, "09_quick_check_EBV_test_rawBarcodes_GSE189141_viralPositiveBarcodes_raw.tsv"))

# attach the same group used in the summary plot (the dataset names are the ones of df_plot) visium has is_called_tissue instead of is_called_cell
df_hits_group_spinal <- df_hits_spinal %>%
  select(sample_id, gene, barcode, count, is_called_cell) %>%
  left_join(df_summary_spinal %>% select(sample_id, group = diagnosis), by = "sample_id") %>%
  mutate(dataset = "Spinal cord snRNAseq", group = as.character(group))

df_hits_group_brain <- df_hits_brain %>%
  select(sample_id, gene, barcode, count, is_called_cell) %>%
  left_join(df_summary_brain %>% select(sample_id, group = pathology_class), by = "sample_id") %>%
  mutate(dataset = "Brain WM/CX snRNAseq", group = as.character(group))

df_hits_group_visium <- df_hits_visium %>%
  select(sample_id, gene, barcode, count, is_called_cell = is_called_tissue) %>%
  left_join(df_summary_brain_visium %>% select(sample_id, group = sample_classification2), by = "sample_id") %>%
  mutate(dataset = "Brain Visium", group = as.character(group))

df_hits_group_gse189141 <- df_hits_gse189141 %>%
  select(sample_id, gene, barcode, count, is_called_cell) %>%
  left_join(df_summary_gse189141 %>% select(sample_id, group = time_point), by = "sample_id") %>%
  mutate(dataset = "GSE189141 (time point)", group = as.character(group))

df_hits_all <- bind_rows(df_hits_group_spinal, df_hits_group_brain, df_hits_group_visium, df_hits_group_gse189141)

# viral UMI and number of positive barcodes per gene, dataset and group
df_gene <- df_hits_all %>%
  group_by(dataset, gene) %>%
  summarise(n_barcodes = n_distinct(sample_id, barcode),
            viral_UMI = sum(count),
            n_called = sum(is_called_cell, na.rm = TRUE),
            .groups = "drop") %>%
  left_join(df_gene_class, by = "gene")

# every gene has to be in the panel and the UMI have to match the per sample summaries

df_gene %>%
  group_by(dataset) %>%
  summarise(viral_UMI = sum(viral_UMI)) %>%
  left_join(bind_rows(
    summary_spinal %>% summarise(total_viral_UMI = sum(total_viral_UMI)) %>% mutate(dataset = "Spinal cord snRNAseq"),
    summary_brain %>% summarise(total_viral_UMI = sum(total_viral_UMI)) %>% mutate(dataset = "Brain WM/CX snRNAseq"),
    summary_visium %>% summarise(total_viral_UMI = sum(total_viral_UMI)) %>% mutate(dataset = "Brain Visium"),
    summary_gse189141 %>% summarise(total_viral_UMI = sum(total_viral_UMI)) %>% mutate(dataset = "GSE189141 (time point)")),
    by = "dataset") %>%
  mutate(match = viral_UMI == total_viral_UMI)

# order the genes by class (panel order) and then by total UMI across the datasets
order_gene <- df_gene %>%
  group_by(viral_class, gene) %>%
  summarise(tot = sum(viral_UMI), .groups = "drop") %>%
  arrange(factor(viral_class, levels = names(list_markers_viral)), desc(tot)) %>%
  pull(gene)

df_gene <- df_gene %>%
  mutate(gene = factor(gene, levels = order_gene),
         viral_class = factor(viral_class, levels = names(list_markers_viral)))

write_tsv(df_gene, file.path(dir_table, "14_quick_check_EBV_test_perGene_summary.tsv"))

# one plot per dataset (free x scale), only the genes with at least a hit
# the panel heights are proportional to the number of genes of each dataset
n_gene <- df_gene %>%
  group_by(dataset) %>%
  summarise(n = n_distinct(gene))

# fixed color per viral class, defined a priori so it is the same in every panel
col_class <- c(EBV_LATENT = "#E41A1C",
               EBV_LYTIC_IE = "#377EB8",
               EBV_LYTIC_E = "#4DAF4A",
               EBV_LYTIC_L = "#984EA3",
               EBV_OTHER = "grey60")

list_p_gene <- map(n_gene$dataset, function(x) {
  df_gene %>%
    filter(dataset == x) %>%
    ggplot(aes(y = gene, x = viral_UMI, fill = viral_class)) +
    geom_col() +
    # first gene on top
    scale_y_discrete(limits = rev) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_fill_manual(values = col_class) +
    labs(x = "EBV UMI", y = NULL, fill = "viral class", title = x) +
    theme_bw() +
    theme(strip.background = element_blank())
})

# keep the legend only of the dataset with the most viral classes, drop it from the others with a single legend left, plot_layout(guides = "collect") places it on the side of the whole figure
n_class <- df_gene %>%
  group_by(dataset) %>%
  summarise(n_class = n_distinct(viral_class))
i_legend <- which.max(n_class$n_class)

# loop only the plot that do not have to top legend
list_p_gene[-i_legend] <- map(list_p_gene[-i_legend], function(x){
  x + theme(legend.position = "none")
})

p_gene <- wrap_plots(list_p_gene, ncol = 1, heights = n_gene$n) +
  plot_layout(guides = "collect")
p_gene

ggsave(file.path("../../out/plot/analysis_R45_pixi", "14_quick_check_EBV_test_perGene_perDataset_scaled.pdf"), p_gene, width = 18, height = 14)

# also plot the unscaled version
p_gene2 <- df_gene %>%
  ggplot(aes(y = gene, x = viral_UMI, fill = viral_class)) +
  geom_col() +
  # first gene on top
  scale_y_discrete(limits = rev) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = col_class) +
  # wrap the long dataset names so the horizontal strips stay narrow
  facet_grid(dataset ~ ., scales = "free_y", space = "free_y", labeller = labeller(dataset = label_wrap_gen(15))) +
  theme_bw() +
  # horizontal strip text (the default on the right is rotated by -90 degrees)
  theme(strip.background = element_blank(),
        strip.text.y = element_text(angle = 0))
ggsave(file.path("../../out/plot/analysis_R45_pixi", "14_quick_check_EBV_test_perGene_perDataset_unscaled.pdf"), plot = p_gene2, width = 18, height = 14)

# barcode level annotation (spinal and brain nuclei, visium brain) ------------------------
# the viral positive barcodes (one row per barcode x gene, read above) are matched to the fully annotated objects (for visium the manual annotation of the spots).
# only the barcodes that also passed the QC of the annotated object (different pipeline run, same FASTQs) get an annotation, the others stay NA

# keep only the metadata of the annotated objects and free the objects right away
sobj_spinal <- readRDS("../../out/object/analysis_R45_pixi/03_sobj_integrated_cleanup_manualAnnotation_subclusterFiltered_v2.rds")
DimPlot(sobj_spinal)
meta_spinal <- sobj_spinal@meta.data %>% rownames_to_column("cell_name")
rm(sobj_spinal); gc()

# meta_spinal %>%
#   group_by(sample_id,GM_prop) %>%
#   summarise() %>%
#   left_join(lut_spinal %>%
#               select(sample_id,GM_prop),by = "sample_id") %>%
#   mutate(test = GM_prop.x == GM_prop.y ) %>%
#   filter(test ==F)

sobj_brain <- readRDS("/beegfs/scratch/ric.cosr/pedrini.edoardo/project_edoardo/220501_scRNAseq_MSbrain_Absinta/out/object/revision/120_WMCX_ManualClean4_harmonySkipIntegration_AllSoupX_4000_AnnotationSCType_manualAnnotation.rds")
DimPlot(sobj_brain)
meta_brain <- sobj_brain@meta.data %>% rownames_to_column("cell_name")
rm(sobj_brain); gc()

# Spinal: cell names are <sample_id>_<barcode> (e.g. Sample_1_AAACCCACACATAGCT-1), same sample_id as the pipeline sample
df_hits_spinal_annotated <- df_hits_spinal %>%
  mutate(cell_name = paste0(sample_id, "_", barcode)) %>%
  left_join(meta_spinal, by = "cell_name", suffix = c("", "_annot")) %>%
  dplyr::select(sample_id:in_postQC,cell_name,scDblFinder.class,orig.ident,diagnosis_short,pathological_stage_short,demye_WM_class,demye_GM_class,demye_tot_class,cell_id_subcluster,cell_id_subcluster2,cell_id_subcluster3)

# Brain: the barcode alone is not unique across samples, so match on the sample key + the 16-mer barcode (with the -1 suffix) the sample key is the official id (s1, s18, ...):
# - in the metadata of the object it is the prefix of the cell name (e.g. s18_<barcode>)
# - in the hits it is baked in the pipeline sample name (e.g. GSM5470490__s18__Homo_sapiens__RNA-Seq, GSM8522356__Cortex__s35__Homo_sapiens__RNA-Seq)
# n_match flags the barcodes that are still ambiguous within a sample (n_match > 1)
meta_brain_barcode <- meta_brain %>%
  mutate(sample_key = str_extract(cell_name, "^s\\d+"),
         barcode = str_extract(cell_name, "[ACGT]{16}-*."))

# make sure all the barcode are available
meta_brain_barcode %>%
  filter(is.na(sample_key))
meta_brain_barcode %>%
  filter(is.na(barcode))

df_hits_brain_annotated <- df_hits_brain %>%
  mutate(sample_key = str_extract(sample_id, "__(s\\d+)__", group = 1)) %>%
  left_join(meta_brain_barcode, by = c("sample_key", "barcode"), suffix = c("", "_annot")) %>%
  dplyr::select(sample_id:in_postQC,cell_name,scDblFinder.class,orig.ident,pathology_class,expertAnno.l1,expertAnno.l2)

# Visium brain: the manual annotation of the spots is one csv per slide in the "fix" folder (columns Barcode, manual_anno, manula_anno2) the files are named after the sample_id of the LUT (e.g. 02_SP1) while the hits use the library_id (e.g. MA7678_2), so the LUT is used for the conversion S18_044_manual_annotation.csv is not a slide of the visium brain LUT, therefore only the *_SP1_ files are read
files_visium <- list.files(file.path(dir_project, "260216_visium_brain_absinta/data/manual_annotation/fix"),
                           pattern = "_SP1_manual_annotation\\.csv$", full.names = TRUE)
names(files_visium) <- str_remove(basename(files_visium), "_manual_annotation\\.csv$")

meta_visium <- map_dfr(files_visium, read_csv, .id = "sample_id_lut") %>%
  # fix the typo of the column name
  rename(barcode = Barcode, manual_anno2 = manula_anno2) %>%
  left_join(select(lut_visium, library_id, sample_id), by = c("sample_id_lut" = "sample_id"))

# the sample_id of the hits is the library_id. sample_id_lut is not NA only for the spots found in the annotation files (manual_anno can be NA for a spot that is in the file)
df_hits_visium_annotated <- df_hits_visium %>%
  left_join(meta_visium, by = c("sample_id" = "library_id", "barcode"), relationship = "many-to-one") %>%
  left_join(select(lut_visium, library_id, sample_classification, sample_classification2), by = c("sample_id" = "library_id")) %>%
  dplyr::select(sample_id:is_called_tissue, sample_id_lut, manual_anno, manual_anno2, sample_classification, sample_classification2)

# how many hits could be annotated
df_hits_spinal_annotated %>% summarise(n_hits = n(), n_annotated = sum(cell_name %in% meta_spinal$cell_name))
df_hits_brain_annotated %>% summarise(n_hits = n(), n_annotated = sum(!is.na(cell_name)))
df_hits_visium_annotated %>% summarise(n_hits = n(), n_in_annotation = sum(!is.na(sample_id_lut)), n_with_label = sum(!is.na(manual_anno)))

# save the tables
write_tsv(df_hits_spinal_annotated, file.path(dir_table, "14_quick_check_EBV_test_rawBarcodes_Spinal_viralPositiveBarcodes_annotated.tsv"))
write_tsv(df_hits_brain_annotated, file.path(dir_table, "14_quick_check_EBV_test_rawBarcodes_AbsintaBrainWMCX_viralPositiveBarcodes_annotated.tsv"))
write_tsv(df_hits_visium_annotated, file.path(dir_table, "14_quick_check_EBV_test_rawBarcodes_VisiumBrain_viralPositiveBarcodes_annotated.tsv"))
