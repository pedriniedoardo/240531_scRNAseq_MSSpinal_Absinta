# AIM ---------------------------------------------------------------------
# this is just an update for the script @26_propeller.R to inplement the suggestion from Aletta to remove the control samples form the analysis of the proportion for the pathology
# notice that to work on the pathology cavariate we could even work with the previous version 26_sobj_integrated_cleanup_manualAnnotation.rds

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(speckle)
library(limma)
library(statmod)
library(cowplot)
library(ggrepel)
library(finalfit)
library(glmmTMB)
library(emmeans)

# read in the data --------------------------------------------------------
scobj <- readRDS("../../out/object/analysis_R44/001_sobj_integrated_cleanup_manualAnnotation_updateDemyelination.rds")

# confirm the object and annotations
DimPlot(scobj,group.by="cell_id")

meta_ref <- scobj@meta.data %>%
  rownames_to_column()

# Aletta suggested removing the CTRL sample from the analysis
meta_ref_testAletta <- meta_ref %>%
  filter(pathological_stage_short != "CTRL") %>%
  # drop the level CTRL
  mutate(pathological_stage_short = fct_drop(pathological_stage_short))

meta_ref_testAletta$pathological_stage_short %>% summary()

# diagnosis ---------------------------------------------------------------
# confirm the numbers from the tissue dataset
meta_ref_testAletta %>% 
  group_by(sample_id,RNA_snn_res.0.1,diagnosis) %>% 
  summarise(n = n())

# run the proportion test pathology stage ---------------------------------
# Run propeller testing for cell type proportion differences between the groups.
# cluster is the cluster/celltype id
# sample is the reference id of the biological replicates
# group is the grouping id
out_pathology <- propeller(clusters = meta_ref_testAletta$cell_id,
                           sample = meta_ref_testAletta$sample_id,
                           group = meta_ref_testAletta$pathological_stage_short)

out_pathology %>%
  rownames_to_column("cell_id") %>%
  write_tsv("../../out/table/analysis_R44/002_propeller_out_PatStage_cellid_testAletta.tsv")

# plotting diagnosis ------------------------------------------------------
# implementation to account for missing cells combinations
df_summary_pathology <- meta_ref_testAletta %>%
  group_by(cell_id, sample_id, pathological_stage_short) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = cell_id, values_from = n, values_fill = 0) %>%
  pivot_longer(names_to = "cell_id", values_to = "n", -c(sample_id, pathological_stage_short)) %>%
  group_by(sample_id) %>% 
  mutate(tot = sum(n), prop = n/tot) %>%
  ungroup()

# save the aggregated table of values
df_summary_pathology %>%
  write_tsv("../../out/table/analysis_R44/002_df_summary_pathology_testAletta.tsv")

# plot 01
df_summary_pathology %>%
  ggplot(aes(x=pathological_stage_short,y=prop))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),shape=1,alpha =0.7)+
  facet_wrap(~cell_id,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/plot/analysis_R44/002_propeller_plot01_PatStage_cellid_testAletta.pdf",width = 10,height = 10)

# attempt modelling -------------------------------------------------------
# as we tried in the @005_scPorportionTest_alletta.R, for modelling the diagnosis, I am trying to model the proportion data using the mixed model.
# since after removign the CTRL variable, the covarite has become binary, this is quite straightforward to implement
# for the beta family modelling tweak the absolute 0 and 1 into sometihng slightly below
# pull the minimumnb value

min_prop <- df_summary_pathology %>%
  ungroup() %>%
  filter(prop > 0) %>%
  summarise(min_prop = min(prop)) %>%
  pull(min_prop)

df_summary_pathology <- df_summary_pathology %>%
  mutate(prop_fix = case_when(prop == 0 ~ (0 + min_prop/100),
                              prop == 1 ~ (1 - min_prop/100),
                              T ~ prop))

# compare prop and prop_fix
df_summary_pathology %>%
  ggplot(aes(x = prop,y = prop_fix)) + geom_point()

# how many where changed
df_summary_pathology %>%
  filter(prop == 0 | prop == 1) %>%
  summarise(n = n())

model <- glmmTMB(prop_fix ~ cell_id * pathological_stage_short + (1 | sample_id),
                 data = df_summary_pathology,
                 family = beta_family(link = "logit"))

summary(model)

# make the stiamates base on the model
emm <- emmeans(model, ~ pathological_stage_short | cell_id, type = "response")

# make the pairs split by the comparisons of interest
res <- pairs(emm,reverse = T)

# add the confidence interval
df_res <- left_join(
  res %>%
    # Adds lower.CL and upper.CL columns
    confint() %>%
    data.frame(),
  res %>%
    data.frame() %>%
    select(contrast,cell_id,null,z.ratio,p.value), by = c("contrast","cell_id")) %>%
  arrange(p.value)

# make the plot
df_res %>%
  mutate(cell_id = fct_reorder(cell_id,odds.ratio)) %>%
  ggplot(aes(y=cell_id,x=odds.ratio)) +
  geom_point() +
  geom_errorbar(aes(xmin = asymp.LCL,
                    xmax = asymp.UCL),width = 0.15, linewidth = 0.5) +
  geom_vline(col="red",linetype = "dashed",xintercept = 1) +
  theme_bw() +
  scale_x_log10()

# calculate the estimates and the odds ration in the same table
emm %>%
  data.frame() %>%
  dplyr::select(-df) %>%
  pivot_wider(names_from = pathological_stage_short,
              values_from = c(response,SE,asymp.LCL,asymp.UCL)) %>%
  left_join(res %>% data.frame(),by = c("cell_id")) %>%
  # manual calculation fo the odds ration: 
  # OR = [response_ACT / (1 - response_ACT)] / [response_IN / (1 - response_IN)]
  # mutate(OR_calc = (response_ACT/(1-response_ACT))/(response_IN/(1-response_IN))) %>%
  mutate(delta_prop_abs = response_ACT - response_IN,
         delta_prop_rel = delta_prop_abs / response_IN) %>%
  arrange(desc(abs(delta_prop_rel)))
