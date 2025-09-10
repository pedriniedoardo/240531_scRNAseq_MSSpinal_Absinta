# AIM ---------------------------------------------------------------------
# run a statistical test on the proportion differences across clusters, between conditions

# libraries ---------------------------------------------------------------
library(tidyverse)
library(speckle)
library(limma)
library(statmod)
library(cowplot)
library(ggrepel)
library(finalfit)

# read in the data --------------------------------------------------------
scobj <- readRDS("../../out/object/analysis_R44/26_sobj_integrated_cleanup_manualAnnotation.rds")

# # add the cell_id based on martina's annotation
# scobj$cell_id <- scobj@meta.data %>%
#   mutate(cell_id = case_when(seurat_clusters %in% c(0,4,14)~"OLIGO",
#                              seurat_clusters %in% c(9)~"OPC",
#                              seurat_clusters %in% c(3,12)~"ASTRO",
#                              seurat_clusters %in% c(5)~"IMMUNE",
#                              seurat_clusters %in% c(13)~"LYM",
#                              seurat_clusters %in% c(11)~"VAS",
#                              seurat_clusters %in% c(1, 2, 10,6)~"EXC NEU",
#                              seurat_clusters %in% c(7,8)~"INH NEU")) %>%
#   pull(cell_id)


# meta <- scobj@meta.data %>%
#   rownames_to_column()
# write_tsv(meta,"../../out/table/meta_data.combined_harmonySkipIntegration_AllSoupX_00500_07000_05_AnnotationSCType.tsv")

meta_ref <- scobj@meta.data %>%
  rownames_to_column() %>%
  # mutate(RNA_snn_res.0.2 = fct_relevel(as.character(RNA_snn_res.0.2),as.character(0:15))) %>%
  # mutate(diagnosis = fct_relevel(diagnosis,c("Non-demented control","Multiple sclerosis"))) %>%
  mutate(diagnosis = fct_relevel(diagnosis_short,c("CTRL","MS"))) %>%
  mutate(location = fct_relevel(location,c("cervical","thoracic","lumbar"))) %>%
  # mutate(pathological_stage = fct_relevel(pathological_stage,c("control","inactive","active demyelination"))) %>%
  # mutate(demyelination = fct_relevel(demyelination,c("no demyelination","less than half","more than half")))
  mutate(pathological_stage = fct_relevel(pathological_stage_short,c("CTRL","IN","ACT"))) %>%
  mutate(demyelination = fct_relevel(demyelination_short,c("NO","LESS50","MORE50")))

# diagnosis ---------------------------------------------------------------
# confirm the numbers from the tissue dataset
meta_ref %>% 
  group_by(sample_id,RNA_snn_res.0.1,diagnosis) %>% 
  summarise(n = n())

# run the proportion test diagnosis ---------------------------------------
# Run propeller testing for cell type proportion differences between the groups.
# cluster is the cluster/celltype id
# sample is the reference id of the biological replicates
# group is the grouping id
out_diagnosis <- propeller(clusters = meta_ref$cell_id,
                           sample = meta_ref$sample_id,
                           group = meta_ref$diagnosis)

out_diagnosis %>%
  rownames_to_column("cell_id") %>%
  write_tsv("../../out/table/analysis_R44/26_propeller_out_diagnosis_cellid.tsv")

# plotting diagnosis ------------------------------------------------------
df_summary_diagnosis <- meta_ref %>% 
  group_by(cell_id,
           sample_id,
           diagnosis) %>% 
  summarise(n=n()) %>% 
  ungroup() %>% 
  group_by(sample_id) %>% 
  mutate(tot = sum(n),
         prop = n/tot)

# plot 01
df_summary_diagnosis %>%
  ggplot(aes(x=diagnosis,y=prop))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),shape=1,alpha =0.7)+
  facet_wrap(~cell_id,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/plot/analysis_R44/26_propeller_plot01_diagnosis_cellid.pdf",width = 10,height = 10)

# plot 02
df_summary_diagnosis %>%
  ggplot() +
  geom_boxplot(aes(x=cell_id,y=prop,color=diagnosis),outlier.shape = NA) +
  geom_point(aes(x=cell_id,y=prop,color=diagnosis),position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.8),alpha=0.7) +
  theme_cowplot()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  scale_y_sqrt()
ggsave("../../out/plot/analysis_R44/26_propeller_plot02_diagnosis_cellid.pdf",width = 8,height = 5)

# location ----------------------------------------------------------------
# confirm the numbers from the tissue dataset
meta_ref %>% 
  group_by(sample_id,cell_id,location) %>% 
  summarise(n = n())

# run the proportion test location ----------------------------------------
# Run propeller testing for cell type proportion differences between the groups.
# cluster is the cluster/celltype id
# sample is the reference id of the biological replicates
# group is the grouping id
out_location <- propeller(clusters = meta_ref$cell_id,
                          sample = meta_ref$sample_id,
                          group = meta_ref$location)

out_location %>%
  rownames_to_column("cell_id") %>%
  write_tsv("../../out/table/analysis_R44/26_propeller_out_location_cellid.tsv")

# plotting diagnosis ------------------------------------------------------
df_summary_location <- meta_ref %>% 
  group_by(cell_id,
           sample_id,
           location) %>% 
  summarise(n=n()) %>% 
  ungroup() %>% 
  group_by(sample_id) %>% 
  mutate(tot = sum(n),
         prop = n/tot)

# plot 01
df_summary_location %>%
  ggplot(aes(x=location,y=prop))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),shape=1,alpha =0.7)+
  facet_wrap(~cell_id,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/plot/analysis_R44/26_propeller_plot01_location_cellid.pdf",width = 10,height = 10)

# plot 02
df_summary_location %>%
  ggplot() +
  geom_boxplot(aes(x=cell_id,y=prop,color=location),outlier.shape = NA) +
  geom_point(aes(x=cell_id,y=prop,color=location),position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.8),alpha=0.7) +
  theme_cowplot()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  scale_y_sqrt()
ggsave("../../out/plot/analysis_R44/26_propeller_plot02_location_cellid.pdf",width = 8,height = 5)

df_summary_location2 <- meta_ref %>% 
  group_by(cell_id,
           sample_id,
           location,
           diagnosis) %>% 
  summarise(n=n()) %>% 
  ungroup() %>% 
  group_by(sample_id) %>% 
  mutate(tot = sum(n),
         prop = n/tot)

# plot 01
df_summary_location2 %>%
  ggplot(aes(x=location,y=prop,color = diagnosis))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.8),shape=1,alpha =0.7)+
  facet_wrap(~cell_id,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/plot/analysis_R44/propeller_plot01alt_location_cellid.pdf",width = 10,height = 10)

# run the proportion test pathology stage ---------------------------------
# Run propeller testing for cell type proportion differences between the groups.
# cluster is the cluster/celltype id
# sample is the reference id of the biological replicates
# group is the grouping id
out_pathology <- propeller(clusters = meta_ref$cell_id,
                           sample = meta_ref$sample_id,
                           group = meta_ref$pathological_stage)

out_pathology %>%
  rownames_to_column("cell_id") %>%
  write_tsv("../../out/table/analysis_R44/26_propeller_out_PatStage_cellid.tsv")

# plotting diagnosis ------------------------------------------------------
df_summary_pathology <- meta_ref %>% 
  group_by(cell_id,
           sample_id,
           pathological_stage) %>% 
  summarise(n=n()) %>% 
  ungroup() %>% 
  group_by(sample_id) %>% 
  mutate(tot = sum(n),
         prop = n/tot)

# plot 01
df_summary_pathology %>%
  ggplot(aes(x=pathological_stage,y=prop))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),shape=1,alpha =0.7)+
  facet_wrap(~cell_id,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/plot/analysis_R44/26_propeller_plot01_PatStage_cellid.pdf",width = 10,height = 10)

# plot 02
df_summary_pathology %>%
  ggplot() +
  geom_boxplot(aes(x=cell_id,y=prop,color=pathological_stage),outlier.shape = NA) +
  geom_point(aes(x=cell_id,y=prop,color=pathological_stage),position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.8),alpha=0.7) +
  theme_cowplot()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  scale_y_sqrt()
ggsave("../../out/plot/analysis_R44/26_propeller_plot02_PatStage_cellid.pdf",width = 8,height = 5)

df_summary_pathology2 <- meta_ref %>% 
  group_by(cell_id,
           sample_id,
           pathological_stage,
           location) %>% 
  summarise(n=n()) %>% 
  ungroup() %>% 
  group_by(sample_id) %>% 
  mutate(tot = sum(n),
         prop = n/tot)

# plot 01
df_summary_pathology2 %>%
  ggplot(aes(x=pathological_stage,y=prop,color = location))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.8),shape=1,alpha =0.7)+
  facet_wrap(~cell_id,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/plot/analysis_R44/26_propeller_plot01alt_PatStage_cellid.pdf",width = 10,height = 10)

# run the proportion test demyelination -----------------------------------
# Run propeller testing for cell type proportion differences between the groups.
# cluster is the cluster/celltype id
# sample is the reference id of the biological replicates
# group is the grouping id
out_dem <- propeller(clusters = meta_ref$cell_id,
                     sample = meta_ref$sample_id,
                     group = meta_ref$demyelination)

out_dem %>%
  rownames_to_column("cell_id") %>%
  write_tsv("../../out/table/analysis_R44/26_propeller_out_dem_cellid.tsv")

# plotting diagnosis ------------------------------------------------------
df_summary_dem <- meta_ref %>% 
  group_by(cell_id,
           sample_id,
           demyelination) %>% 
  summarise(n=n()) %>% 
  ungroup() %>% 
  group_by(sample_id) %>% 
  mutate(tot = sum(n),
         prop = n/tot)

# plot 01
df_summary_dem %>%
  ggplot(aes(x=demyelination,y=prop))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),shape=1,alpha =0.7)+
  facet_wrap(~cell_id,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/plot/analysis_R44/26_propeller_plot01_dem_cellid.pdf",width = 10,height = 10)

# plot 02
df_summary_dem %>%
  ggplot() +
  geom_boxplot(aes(x=cell_id,y=prop,color=demyelination),outlier.shape = NA) +
  geom_point(aes(x=cell_id,y=prop,color=demyelination),position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.8),alpha=0.7) +
  theme_cowplot()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  scale_y_sqrt()
ggsave("../../out/plot/analysis_R44/26_propeller_plot02_dem_cellid.pdf",width = 8,height = 5)

df_summary_dem2 <- meta_ref %>% 
  group_by(cell_id,
           sample_id,
           demyelination,
           location) %>% 
  summarise(n=n()) %>% 
  ungroup() %>% 
  group_by(sample_id) %>% 
  mutate(tot = sum(n),
         prop = n/tot)

# plot 01
df_summary_dem2 %>%
  ggplot(aes(x=demyelination,y=prop,color = location))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.8),shape=1,alpha =0.7)+
  facet_wrap(~cell_id,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/plot/analysis_R44/26_propeller_plot01alt_dem_cellid.pdf",width = 10,height = 10)

# try to use a linear model to account for more metadata
test <- meta_ref %>% 
  group_by(cell_id,
           sample_id,
           diagnosis_short,
           pathological_stage_short,
           demyelination_short,
           location) %>% 
  summarise(n=n()) %>% 
  ungroup() %>% 
  group_by(sample_id) %>% 
  mutate(tot = sum(n),
         prop = n/tot,
         log_ratio = log(prop/(1-prop)))

# confirm that in B cells there is a significant difference
# use the logratio to avoid issue with the proportions

# test with wilcoxon rank sum test
test %>%
  filter(cell_id == "B CELLS") %>%
  data.frame() %>%
  wilcox.test(prop ~ diagnosis_short, data = .)

# test the diagnosis
# Using a linear model directly on the proportions (prop ~ diagnosis_short) is not recommended. The main reason is that linear regression assumes that the errors (residuals) are normally distributed and have constant variance. Proportion data is bounded between 0 and 1, so these assumptions are often violated. This can lead to unreliable p-values and confidence intervals.

test %>%
  filter(cell_id == "B CELLS") %>%
  data.frame() %>%
  finalfit(formula = log_ratio ~ diagnosis_short)

# account for the location
test %>%
  filter(cell_id == "B CELLS") %>%
  data.frame() %>%
  finalfit(formula = log_ratio ~ location)

test %>%
  filter(cell_id == "B CELLS") %>%
  data.frame() %>%
  finalfit(formula = log_ratio ~ location + diagnosis_short)

# -------------------------------------------------------------------------
test %>%
  filter(cell_id == "IMMUNE") %>%
  data.frame() %>%
  finalfit(formula = log_ratio ~ pathological_stage_short)

test %>%
  filter(cell_id == "IMMUNE") %>%
  data.frame() %>%
  lm(formula = log_ratio ~ pathological_stage) %>%
  summary()

test %>%
  filter(cell_id == "IMMUNE") %>%
  data.frame() %>%
  lm(formula = log_ratio ~ pathological_stage) %>%
  summary()

test %>%
  filter(cell_id == "IMMUNE") %>%
  data.frame() %>%
  finalfit(formula = log_ratio ~ pathological_stage)

test %>%
  filter(cell_id == "IMMUNE") %>%
  data.frame() %>%
  finalfit(formula = log_ratio ~ pathological_stage+location)

test %>%
  filter(cell_id == "IMMUNE") %>%
  data.frame() %>%
  lm(formula = log_ratio ~ pathological_stage+location) %>%
  summary()

test %>%
  filter(cell_id == "IMMUNE") %>%
  data.frame() %>%
  finalfit(formula = log_ratio ~ demyelination+location)

test %>%
  filter(cell_id == "IMMUNE") %>%
  data.frame() %>%
  lm(formula = log_ratio ~ demyelination+location+pathological_stage) %>%
  model.matrix()

test %>%
  filter(cell_id == "IMMUNE") %>%
  data.frame() %>%
  lm(formula = log_ratio ~ demyelination+location+pathological_stage) %>%
  summary()

test %>%
  filter(cell_id == "IMMUNE") %>%
  data.frame() %>%
  group_by(demyelination,location,pathological_stage) %>%
  summarise(n = n())

# try finalfit with categorical outcome variable
meta_ref %>%
  # filter(cell_id == "IMMUNE") %>%
  summary_factorlist(dependent = "pathological_stage",explanatory = "cell_id",p=T)
