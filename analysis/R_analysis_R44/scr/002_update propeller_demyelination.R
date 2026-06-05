# AIM ---------------------------------------------------------------------
# run a statistical test on the proportion differences across clusters, between conditions
# this is specifically to update the demyeliantion status covariate

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(speckle)
library(limma)
library(statmod)
library(cowplot)
library(ggrepel)
library(finalfit)
library(ggpubr)
library(broom)

# read in the data --------------------------------------------------------
scobj <- readRDS("../../out/object/analysis_R44/001_sobj_integrated_cleanup_manualAnnotation_updateDemyelination.rds")

# confirm the object and annotations
DimPlot(scobj,group.by="cell_id")

meta_ref <- scobj@meta.data %>%
  rownames_to_column()

# wrangling ---------------------------------------------------------------
# compare the new classifications
meta_ref %>%
  group_by(sample_id,demye_tot_class,demye_WM_class) %>%
  summarise() %>%
  ungroup() %>%
  group_by(demye_tot_class,demye_WM_class) %>%
  summarise(n = n())

meta_ref %>%
  group_by(sample_id,demye_GM_class,demye_WM_class) %>%
  summarise() %>%
  ungroup() %>%
  group_by(demye_GM_class,demye_WM_class) %>%
  summarise(n = n())

# find the samples that swithc labels
meta_ref %>%
  group_by(sample_id,demye_GM_class,demye_WM_class) %>%
  summarise() %>%
  filter(demye_GM_class != demye_WM_class)

# compare the old with the new classifcation
df_summary_dem_WMnew <- meta_ref %>% 
  group_by(cell_id,
           sample_id,
           demye_WM_class) %>% 
  summarise(n=n()) %>% 
  ungroup() %>% 
  filter(!is.na(demye_WM_class)) %>%
  group_by(sample_id) %>% 
  mutate(tot = sum(n),
         prop = n/tot)

df_summary_dem_GMnew <- meta_ref %>% 
  group_by(cell_id,
           sample_id,
           demye_GM_class) %>% 
  summarise(n=n()) %>% 
  ungroup() %>% 
  filter(!is.na(demye_GM_class)) %>%
  group_by(sample_id) %>% 
  mutate(tot = sum(n),
         prop = n/tot)

df_summary_dem_old <- meta_ref %>% 
  group_by(cell_id,
           sample_id,
           demyelination_short) %>% 
  summarise(n=n()) %>% 
  ungroup() %>% 
  filter(!is.na(demyelination_short)) %>%
  group_by(sample_id) %>% 
  mutate(tot = sum(n),
         prop = n/tot)

# plot by group
p01 <- df_summary_dem_WMnew %>%
  filter(cell_id %in% c("B CELLS","T CELLS")) %>%
  ggplot(aes(x=demye_WM_class,y=prop))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),shape=1,alpha =0.7)+
  facet_wrap(~cell_id,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))

p02 <- df_summary_dem_GMnew %>%
  filter(cell_id %in% c("B CELLS","T CELLS")) %>%
  ggplot(aes(x=demye_GM_class,y=prop))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),shape=1,alpha =0.7)+
  facet_wrap(~cell_id,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))

p03 <- df_summary_dem_old %>%
  filter(cell_id %in% c("B CELLS","T CELLS")) %>%
  ggplot(aes(x=demyelination_short,y=prop))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),shape=1,alpha =0.7)+
  facet_wrap(~cell_id,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))

p01 / p02 / p03

# identify the sample swithching labels
sample_switch_WM <- meta_ref %>%
  group_by(sample_id,demyelination_short,demye_WM_class) %>%
  summarise() %>%
  mutate(demyelination_short = fct_recode(demyelination_short,no = "NO",low = "LESS50",high = "MORE50")) %>%
  filter(demyelination_short != demye_WM_class) %>%
  pull(sample_id)

sample_switch_GM <- meta_ref %>%
  group_by(sample_id,demyelination_short,demye_GM_class) %>%
  summarise() %>%
  mutate(demyelination_short = fct_recode(demyelination_short,no = "NO",low = "LESS50",high = "MORE50")) %>%
  filter(demyelination_short != demye_GM_class) %>%
  pull(sample_id)

# show the values for the demyelination scores for the missclassified samples
meta_ref %>%
  group_by(sample_id,demyelination_short,demye_WM_class,demye_WM_prop) %>%
  summarise() %>%
  filter(sample_id %in% sample_switch_WM)

meta_ref %>%
  group_by(sample_id,demyelination_short,demye_GM_class,demye_GM_prop) %>%
  summarise() %>%
  filter(sample_id %in% sample_switch_GM)

# run the proportion test demyelination -----------------------------------
# Run propeller testing for cell type proportion differences between the groups.
# cluster is the cluster/celltype id
# sample is the reference id of the biological replicates

# group is the grouping id
out_dem_GM <- propeller(clusters = meta_ref$cell_id,
                        sample = meta_ref$sample_id,
                        group = meta_ref$demye_GM_class)

out_dem_GM %>%
  rownames_to_column("cell_id") %>%
  write_tsv("../../out/table/analysis_R44/002_propeller_out_dem_cellid_GMclass.tsv")

out_dem_WM <- propeller(clusters = meta_ref$cell_id,
                        sample = meta_ref$sample_id,
                        group = meta_ref$demye_WM_class)

out_dem_WM %>%
  rownames_to_column("cell_id") %>%
  write_tsv("../../out/table/analysis_R44/002_propeller_out_dem_cellid_WMclass.tsv")

# not needed on the tot as the classification is the same as WM
# out_dem_tot <- propeller(clusters = meta_ref$cell_id,
#                          sample = meta_ref$sample_id,
#                          group = meta_ref$demye_tot_class)
# 
# out_dem_tot %>%
#   rownames_to_column("cell_id") %>%
#   write_tsv("../../out/table/analysis_R44/002_propeller_out_dem_cellid_Totclass.tsv")

# plotting diagnosis ------------------------------------------------------
# x <- "demye_WM_class"
lapply(c("demye_WM_class","demye_GM_class"), function(x){
  df_summary_dem <- meta_ref %>% 
    group_by(cell_id,
             sample_id,
             .data[[x]]) %>% 
    summarise(n=n()) %>% 
    ungroup() %>% 
    filter(!is.na(.data[[x]])) %>%
    group_by(sample_id) %>% 
    mutate(tot = sum(n),
           prop = n/tot)
  
  # plot by group split
  df_summary_dem %>%
    ggplot(aes(x=.data[[x]],y=prop))+
    geom_boxplot(outlier.shape = NA)+
    geom_point(position = position_jitter(width = 0.1),shape=1,alpha =0.7)+
    facet_wrap(~cell_id,scales = "free")+
    theme_bw()+
    theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
  
  ggsave(paste0("../../out/plot/analysis_R44/002_propeller_plot01_",x,"_cellid.pdf"),width = 10,height = 10)
  
  # plot by group color
  df_summary_dem %>%
    ggplot() +
    geom_boxplot(aes(x=cell_id,y=prop,color=.data[[x]]),outlier.shape = NA) +
    geom_point(aes(x=cell_id,y=prop,color=.data[[x]]),position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.8),alpha=0.7) +
    theme_cowplot()+
    theme(axis.text.x = element_text(hjust = 1,angle = 90))+
    scale_y_sqrt()
  ggsave(paste0("../../out/plot/analysis_R44/002_propeller_plot02_",x,"_cellid.pdf"),width = 10,height = 10)
})

# with the point estimates we can also plot the correlations values
# x <- "demye_WM_prop"
lapply(c("demye_WM_prop","demye_GM_prop","demye_tot_prop"), function(x){
  df_summary_dem <- meta_ref %>% 
    group_by(cell_id,
             sample_id,
             .data[[x]]) %>% 
    summarise(n=n()) %>% 
    ungroup() %>% 
    filter(!is.na(.data[[x]])) %>%
    group_by(sample_id) %>% 
    mutate(tot = sum(n),
           prop = n/tot)
  
  # plot correlatoin value
  # 1. Compute p-values for the slope per cell_id
  p_values <- df_summary_dem %>%
    group_by(cell_id) %>%
    do(tidy(lm(prop ~ .data[[x]], data = .))) %>%
    filter(term != "(Intercept)") %>% 
    mutate(p_label = paste0("p = ", format.pval(p.value, digits = 3))) %>%
    select(cell_id, p_label)
  
  # 2. Join the labels back to your main data and create a combined header string
  df_plot <- df_summary_dem %>%
    left_join(p_values, by = "cell_id") %>%
    # Create a multi-line string: Name on top, p-value on bottom
    mutate(facet_header = paste0(cell_id, "\n", p_label))
  
  # 3. Plot using the new facet_header variable
  df_plot %>%
    ggplot(aes(x = .data[[x]], y = prop)) +
    geom_point(shape = 1) +
    geom_smooth(method = "lm") +
    # Facet by the new header string instead of just cell_id
    facet_wrap(~facet_header, scales = "free") +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      # Optional: Center or left-align the two-line header text
      strip.text = element_text(hjust = 0.5, face = "bold"), 
      axis.text.x = element_text(hjust = 1, angle = 45)
    )
  
  ggsave(paste0("../../out/plot/analysis_R44/002_scatterProp_plot01_",x,"_cellid.pdf"),width = 10,height = 10)

})

