# AIM ---------------------------------------------------------------------
# plot some genes of interest

# libraries ---------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(ComplexHeatmap)
library(UpSetR)
library(gplots)

# plot genes dotplot ------------------------------------------------------
data <- readRDS("../../out/object/analysis_R44/30_ddsHTSeq_pseudobulk_filterExp_IMMUNE.rds")

lut <- colData(data) %>%
  as.data.frame() %>%
  mutate(sample.id = str_replace_all(sample.id,pattern = "-","."))

# GOI <- read_tsv("../../out/table/res_DIAvsCTRL.txt") %>% 
#   # add a clor variable in case significant
#   mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1),yes = 1,no = 0)) %>%
#   dplyr::filter(col==1) %>%
#   pull(symbol)

GOI <- c("HIST2H2BE", "FAM135B", "HIST2H2BF", "CYCS", "PRKCA", "EME2", "APBB3", "LMO7", "SBF2-AS1", "AC023282.1", "RAMP1")


MR <- counts(data,normalized=T)%>%
  data.frame()%>%
  rownames_to_column("symbol") %>%
  pivot_longer(names_to = "sample",values_to = "exp",-symbol)%>%
  group_by(sample)%>%
  summarise(MR = sum(exp)/10^6)

# plot the data following the methods implemented in the plotCounts funciton from DESeq2
# Normalized counts plus a pseudocount of 0.5 are shown by default.
# sample plot in split by diagnosis
counts(data,normalized=T)%>%
  data.frame()%>%
  rownames_to_column("symbol") %>%
  filter(symbol %in% GOI) %>%
  pivot_longer(names_to = "sample",values_to = "count",-symbol) %>%
  # add the milion reads per sample
  left_join(MR,by = "sample") %>%
  left_join(lut,by = c("sample" = "sample.id")) %>%
  mutate(count_norm_adj = count + 0.5)%>%
  ggplot(aes(x=diagnosis_short,y = count_norm_adj))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha=0.6)+facet_wrap(~symbol,scales = "free")+scale_y_log10()+ theme_bw()+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
ggsave(paste0("../../out/plot/analysis_R44/31_boxplot_GOI01_","IMMUNE",".pdf"),width = 8,height = 6)

# sample plot split by location and diagnosis
counts(data,normalized=T)%>%
  data.frame()%>%
  rownames_to_column("symbol") %>%
  filter(symbol %in% GOI) %>%
  pivot_longer(names_to = "sample",values_to = "count",-symbol) %>%
  # add the milion reads per sample
  left_join(MR,by = "sample") %>%
  left_join(lut,by = c("sample" = "sample.id")) %>%
  mutate(count_norm_adj = count + 0.5)%>%
  ggplot(aes(x=diagnosis_short,y = count_norm_adj,colour = location))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.8),alpha=0.6)+facet_wrap(~symbol,scales = "free")+scale_y_log10()+ theme_bw()+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
ggsave(paste0("../../out/plot/analysis_R44/31_boxplot_GOI01_2_","IMMUNE",".pdf"),width = 8,height = 6)
