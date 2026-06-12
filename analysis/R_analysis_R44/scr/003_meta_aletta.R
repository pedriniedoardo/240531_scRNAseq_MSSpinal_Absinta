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
# aggregate the data per cell numbers
meta_01 <- meta_ref %>%
  distinct(project_id,
           original_sample_name,
           sample_number,
           sample_id,
           nbb,
           autopsy,
           cohort,
           sex,
           age,
           braak,
           amyloid,
           braaklb,
           pmd,
           ph,
           weight,
           csf,
           apoe,
           iduit,
           datumuit,
           uitvraag,
           recipient,
           dcodeprot,
           dcode,
           diagnosis,
           dcodewk,
           wcode,
           region,
           location,
           pathological_stage,
           specific,
           ocode,
           storage,
           diagnosis_short,
           pathological_stage_short,
           sequencing,
           GM_prop,
           demye_WM_prop,
           demye_WM_class,
           demye_GM_prop,
           demye_GM_class,
           demye_tot_prop,
           demye_tot_class)

# if missing combination but 0 count of cells  
# meta_02 <- meta_ref %>%
#   group_by(sample_id,
#            cell_id) %>%
#   summarise(n_cell = n())

# df_meta_full <- crossing(sample_id = meta_02$sample_id %>% unique(),
#          cell_id = meta_02$cell_id %>% unique()) %>%
#   left_join(meta_02,by = c("sample_id","cell_id")) %>%
#   # filter(is.na(n_cell)) %>%
#   mutate(n_cell = case_when(is.na(n_cell) ~ 0,
#                             T ~ n_cell)) %>%
#   left_join(meta_01,by = "sample_id")

# better implementation
meta_02 <- meta_ref %>%
  count(sample_id, cell_id, name = "n_cell") %>%
  complete(sample_id, cell_id, fill = list(n_cell = 0))

df_meta_full <- meta_02 %>%
  left_join(meta_01, by = "sample_id")

# meta_02 %>%
#   filter(n_cell == 0)

df_meta_full %>%
  # filter(is.na(n_cell))
  filter(n_cell == 0)
  
# save
df_meta_full %>%
  write_tsv(file = "../../out/table/analysis_R44/003_meta_spinal_full_aletta.tsv")
