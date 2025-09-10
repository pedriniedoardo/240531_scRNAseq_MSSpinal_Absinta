# AIM ---------------------------------------------------------------------
# run a statistical test on the proportion differences across clusters, between conditions

# libraries ---------------------------------------------------------------
library(tidyverse)
library(speckle)
library(limma)
library(statmod)
library(cowplot)
library(ggrepel)

# read in the data --------------------------------------------------------
# scobj <- readRDS("../../out/object/data.combined_harmonySkipIntegration_AllSoupX_00500_07000_05_AnnotationSCType.rds")
# 
# # # add the cell_id based on martina's annotation
# # scobj$cell_id <- scobj@meta.data %>%
# #   mutate(cell_id = case_when(seurat_clusters %in% c(0,4,14)~"OLIGO",
# #                              seurat_clusters %in% c(9)~"OPC",
# #                              seurat_clusters %in% c(3,12)~"ASTRO",
# #                              seurat_clusters %in% c(5)~"IMMUNE",
# #                              seurat_clusters %in% c(13)~"LYM",
# #                              seurat_clusters %in% c(11)~"VAS",
# #                              seurat_clusters %in% c(1, 2, 10,6)~"EXC NEU",
# #                              seurat_clusters %in% c(7,8)~"INH NEU")) %>%
# #   pull(cell_id)
# 
# 
# meta <- scobj@meta.data %>%
#   rownames_to_column()
# write_tsv(meta,"../../out/table/meta_data.combined_harmonySkipIntegration_AllSoupX_00500_07000_05_AnnotationSCType.tsv")
# 
meta_ref <- read_tsv(file = "../../out/table/meta_data.combined_harmonySkipIntegration_AllSoupX_00500_07000_05_AnnotationSCType.tsv") %>%
  mutate(RNA_snn_res.0.1 = fct_relevel(as.character(RNA_snn_res.0.1),as.character(0:12))) %>%
  mutate(diagnosis = fct_relevel(diagnosis,c("Non-demented control","Multiple sclerosis"))) %>%
  mutate(location = fct_relevel(location,c("cervical","thoracic","lumbar")))

# diagnosis ---------------------------------------------------------------
# confirm the numbers from the tissue dataset
meta_ref %>% 
  group_by(orig.ident,RNA_snn_res.0.1,diagnosis) %>% 
  summarise(n = n())

# run the proportion test diagnosis ---------------------------------------
# Run propeller testing for cell type proportion differences between the groups.
# cluster is the cluster/celltype id
# sample is the reference id of the biological replicates
# group is the grouping id
out_diagnosis <- propeller(clusters = meta_ref$RNA_snn_res.0.1,
                           sample = meta_ref$orig.ident,
                           group = meta_ref$diagnosis)

out_diagnosis %>%
  rownames_to_column("cluster") %>%
  write_tsv("../../out/table/propeller_out_diagnosis_cluster.tsv")

# plotting diagnosis ------------------------------------------------------
df_summary_diagnosis <- meta_ref %>% 
  group_by(RNA_snn_res.0.1,
           orig.ident,
           diagnosis) %>% 
  summarise(n=n()) %>% 
  ungroup() %>% 
  group_by(orig.ident) %>% 
  mutate(tot = sum(n),
         prop = n/tot)

# plot 01
df_summary_diagnosis %>%
  ggplot(aes(x=diagnosis,y=prop))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),shape=1,alpha =0.7)+
  facet_wrap(~RNA_snn_res.0.1,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/plot/propeller_plot01_diagnosis_cluster.pdf",width = 6,height = 10)

# plot 02
df_summary_diagnosis %>%
  ggplot() +
  geom_boxplot(aes(x=RNA_snn_res.0.1,y=prop,color=diagnosis),outlier.shape = NA) +
  geom_point(aes(x=RNA_snn_res.0.1,y=prop,color=diagnosis),position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.8),alpha=0.7) +
  theme_cowplot()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  scale_y_sqrt()
ggsave("../../out/plot/propeller_plot02_diagnosis_cluster.pdf",width = 8,height = 5)

# location ----------------------------------------------------------------
# confirm the numbers from the tissue dataset
meta_ref %>% 
  group_by(orig.ident,RNA_snn_res.0.1,location) %>% 
  summarise(n = n())

# run the proportion test location ----------------------------------------
# Run propeller testing for cell type proportion differences between the groups.
# cluster is the cluster/celltype id
# sample is the reference id of the biological replicates
# group is the grouping id
out_location <- propeller(clusters = meta_ref$RNA_snn_res.0.1,
                          sample = meta_ref$orig.ident,
                          group = meta_ref$location)

out_location %>%
  rownames_to_column("cluster") %>%
  write_tsv("../../out/table/propeller_out_location_cluster.tsv")

# plotting diagnosis ------------------------------------------------------
df_summary_location <- meta_ref %>% 
  group_by(RNA_snn_res.0.1,
           orig.ident,
           location) %>% 
  summarise(n=n()) %>% 
  ungroup() %>% 
  group_by(orig.ident) %>% 
  mutate(tot = sum(n),
         prop = n/tot)

# plot 01
df_summary_location %>%
  ggplot(aes(x=location,y=prop))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),shape=1,alpha =0.7)+
  facet_wrap(~RNA_snn_res.0.1,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/plot/propeller_plot01_location_cluster.pdf",width = 6,height = 10)

# plot 02
df_summary_location %>%
  ggplot() +
  geom_boxplot(aes(x=RNA_snn_res.0.1,y=prop,color=location),outlier.shape = NA) +
  geom_point(aes(x=RNA_snn_res.0.1,y=prop,color=location),position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.8),alpha=0.7) +
  theme_cowplot()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  scale_y_sqrt()
ggsave("../../out/plot/propeller_plot02_location_cluster.pdf",width = 8,height = 5)

df_summary_location2 <- meta_ref %>% 
  group_by(RNA_snn_res.0.1,
           orig.ident,
           location,
           diagnosis) %>% 
  summarise(n=n()) %>% 
  ungroup() %>% 
  group_by(orig.ident) %>% 
  mutate(tot = sum(n),
         prop = n/tot)

# plot 01
df_summary_location2 %>%
  ggplot(aes(x=location,y=prop,color = diagnosis))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.8),shape=1,alpha =0.7)+
  facet_wrap(~RNA_snn_res.0.1,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/plot/propeller_plot01alt_location_cluster.pdf",width = 10,height = 10)

