# AIM ---------------------------------------------------------------------
# same readout as 14_lut_sample.R (EBV hits per dataset/group) but normalised by the sequencing depth of the library.
# the depth term is total_UMI of the *_perSampleSummary.tsv: sum of nCount_RNA over ALL the raw barcodes (empty droplets/ambient included), so it is a per library depth proxy that does not depend on the cell calling.
# the readout is the viral UMI per million raw UMI, computed per sample. group level values are the mean of the per sample rates, so every sample weighs the same regardless of its depth.
# the samples are grouped as control or case (group2) within each dataset, and the difference is tested per dataset with a negative binomial glm (viral UMI count, log(total_UMI) as offset).
# the called cell hits can not be normalised by the raw UMI (different denominator), they are only reported as counts in the tables.
# jakel and shirmer are not included, as in 14 (no EBV reads).
# the datasets are different technologies/depths: the comparisons are made within each dataset, no comparison across datasets is meant.

# libraries ---------------------------------------------------------------
library(tidyverse)
library(patchwork)

# read in the data --------------------------------------------------------
dir_table <- "../../out/table/analysis_R45_pixi"
dir_plot <- "../../out/plot/analysis_R45_pixi"
dir_project <- "/beegfs/scratch/ric.cosr/pedrini.edoardo/project_edoardo"

# per sample EBV summaries
summary_spinal <- read_tsv(file.path(dir_table, "07_quick_check_EBV_test_rawBarcodes_perSampleSummary.tsv"))
summary_brain <- read_tsv(file.path(dir_table, "11_quick_check_EBV_test_rawBarcodes_AbsintaBrainWMCX_perSampleSummary.tsv"))
summary_visium <- read_tsv(file.path(dir_table, "13_quick_check_EBV_test_rawBarcodes_VisiumBrain_perSampleSummary.tsv"))
summary_gse189141 <- read_tsv(file.path(dir_table, "09_quick_check_EBV_test_rawBarcodes_GSE189141_perSampleSummary.tsv"))

# the normalisation is meaningless without the depth
stopifnot(all(c("total_UMI", "total_viral_UMI") %in% colnames(summary_spinal)),
          !anyNA(c(summary_spinal$total_UMI, summary_brain$total_UMI, summary_visium$total_UMI, summary_gse189141$total_UMI)))

# sample metadata
lut_spinal_full <- read_csv(file.path(dir_project,"240531_scRNAseq_MSSpinal_Absinta/data/LUT_sample_full.csv"))
lut_spinal_demye <- read_csv(file.path(dir_project, "240531_scRNAseq_MSSpinal_Absinta/data/20260416_spinal_demyelination_update.csv"))
lut_brain <- read_csv(file.path(dir_project, "220501_scRNAseq_MSbrain_Absinta/data/LUT_sample_WM_CX.csv"))
lut_visium <- read_csv(file.path(dir_project, "260216_visium_brain_absinta/data/LUT_samples.csv"))
# for summary_gse189141 the metadata is backed into the sample name

# build the spinal sample metadata as in 14_lut_sample.R (and in the annotated object)
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

# Brain WM/CX (Absinta): the official id (e.g. s9) is embedded in the pipeline sample name (GSM..._s9)
# the sample_id of the LUT is the original library name, so it is dropped to avoid clashing with the pipeline sample_id
df_summary_brain <- summary_brain %>%
  mutate(official_id = str_extract(sample_id, "_(s\\d+)", group = 1)) %>%
  left_join(select(lut_brain, -sample_id), by = "official_id")

# Visium brain
df_summary_brain_visium <- left_join(summary_visium, lut_visium, by = c("sample_id" = "library_id"))

# GSE189141: GSM_TX<donor>_<timepoint>_GEX
df_summary_gse189141 <- summary_gse189141 %>%
  mutate(gsm = str_extract(sample_id, "^GSM\\d+"),
         donor = str_extract(sample_id, "TX\\d+"),
         time_point = as.numeric(str_extract(sample_id, "_(\\d+)_", group = 1)))

# per sample table --------------------------------------------------------
# one row per sample with the detailed group (group) and the control/case class (group2); same dataset names as 14
df_sample <- bind_rows(
  df_summary_spinal %>% mutate(dataset = "Spinal cord snRNAseq", group = as.character(diagnosis)),
  df_summary_brain %>% mutate(dataset = "Brain WM/CX snRNAseq", group = as.character(pathology_class)),
  df_summary_brain_visium %>% mutate(dataset = "Brain Visium", group = as.character(sample_classification2)),
  df_summary_gse189141 %>% mutate(dataset = "GSE189141 (time point)", group = as.character(time_point))
) %>%
  # group2: control vs case regardless of the level (brain WM and CX together; GSE189141 time point 0 is the baseline)
  mutate(group2 = case_when(group %in% c("Non-demented control", "WM_Ctrl", "CX_Ctrl", "non-lesion", "0") ~ "control",
                            TRUE ~ "case") %>%
           factor(levels = c("control", "case"))) %>%
  select(dataset, group, group2, sample_id, n_raw_barcodes, total_UMI, n_viral_hits, total_viral_UMI, n_called_cell_hits) %>%
  mutate(viral_UMI_per_million = total_viral_UMI / total_UMI * 1e6,
         hits_per_million = n_viral_hits / total_UMI * 1e6)

write_tsv(df_sample, file.path(dir_table, "15_quick_check_EBV_test_perSample_UMInorm.tsv"))

# summarise per group -----------------------------------------------------
# per dataset and control/case: the avg_* columns are the mean of the per sample rates (equal weight per sample), the tot_* columns are raw sums
df_plot <- df_sample %>%
  group_by(dataset,group2) %>%
  summarise(n_samples = n(),
            # samples with at least one viral UMI
            n_pos_samples = sum(total_viral_UMI > 0),
            frac_pos_samples = n_pos_samples / n_samples,
            tot_UMI = sum(total_UMI),
            tot_viral_hits = sum(n_viral_hits),
            tot_viral_UMI = sum(total_viral_UMI),
            tot_called_hits = sum(n_called_cell_hits),
            avg_viral_UMI_per_million = mean(viral_UMI_per_million),
            avg_hits_per_million = mean(hits_per_million),
            .groups = "drop")

write_tsv(df_plot, file.path(dir_table, "15_quick_check_EBV_test_perGroup_UMInorm.tsv"))

# summary plot per dataset ------------------------------------------------
# one panel per dataset, top: bars = mean of the per sample viral UMI per million raw UMI; bottom: the single samples (the datasets have very different sequencing depth)
cap_summary <- "label: n samples"

p_summary <- df_plot %>%
  ggplot(aes(x = group2, y = avg_viral_UMI_per_million, fill = dataset)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = paste0("n=", n_samples)), vjust = -0.2, size = 3) +
  facet_wrap(~dataset, scales = "free", nrow = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.35))) +
  labs(x = NULL, y = "EBV viral UMI per million UMI", caption = cap_summary) +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

p_summary2 <- df_sample %>%
  ggplot(aes(x = group2, y = viral_UMI_per_million, colour = dataset)) +
  geom_point(position = position_jitter(width = 0.2, height = 0), shape = 1) +
  # group mean (same value as the bars above), the median would sit at 0 for most groups
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, linewidth = 0.3) +
  facet_wrap(~dataset, scales = "free", nrow = 1) +
  # a small lower margin so the points at 0 are not cut by the axis
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(x = NULL, y = "EBV viral UMI per million UMI") +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# fraction of samples with at least one viral UMI, label = positive samples / total samples
p_summary3 <- df_plot %>%
  ggplot(aes(x = group2, y = frac_pos_samples, fill = dataset)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = paste0(n_pos_samples, "/", n_samples)), vjust = -0.2, size = 3) +
  facet_wrap(~dataset, nrow = 1) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "fraction of samples with EBV UMI", caption = "label: positive samples / total samples") +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(dir_plot, "15_quick_check_EBV_test_summary_perDataset_UMInorm_scaled.pdf"), plot = (p_summary/p_summary2/p_summary3), width = 12, height = 12)

# case vs control test per dataset ----------------------------------------
# negative binomial glm (MASS::glm.nb) on the viral UMI count (poisson was overdispersed, e.g. GSE189141), the depth (total_UMI) enters as log offset, so the coefficient is the log rate ratio case vs control (exp(coef) = fold change of viral UMI per UMI)
# no random effect: the samples are treated as independent (the repeated donors/patients are ignored, this is exploratory)
# only the datasets with at least 3 samples with a viral UMI are tested, with fewer events the model is not informative
df_glm <- df_sample %>%
  mutate(group2 = factor(group2, levels = c("control", "case"))) %>%
  group_by(dataset) %>%
  mutate(n_pos_sample = sum(total_viral_UMI > 0)) %>%
  ungroup() %>%
  filter(n_pos_sample >= 3) %>%
  split(f = .$dataset) %>%
  imap(function(x, key) {
    fit <- MASS::glm.nb(total_viral_UMI ~ group2 + offset(log(total_UMI)), data = x)
    broom::tidy(fit, conf.int = TRUE) %>%
      filter(term == "group2case") %>%
      mutate(rate_ratio = exp(estimate),
             rate_ratio_low = exp(conf.low),
             rate_ratio_high = exp(conf.high),
             # pearson dispersion as a fit check: about 1 if the negative binomial variance fits (>> 1 = still overdispersed)
             dispersion = sum(residuals(fit, type = "pearson")^2) / fit$df.residual,
             n_control = sum(x$group2 == "control"),
             n_case = sum(x$group2 == "case"))
  }) %>%
  bind_rows(.id = "dataset")
df_glm

write_tsv(df_glm, file.path(dir_table, "15_quick_check_EBV_test_glm_caseVsControl.tsv"))
