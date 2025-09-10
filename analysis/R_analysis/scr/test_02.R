# LIBRARIES ---------------------------------------------------------------
library(scater)
library(Seurat)
library(tidyverse)
library(robustbase)
# library(SeuratData)
library(patchwork)
library(harmony)
library(ggrepel)
library(ComplexHeatmap)

# read in the data --------------------------------------------------------
# load the LUT
# LUT <- read_csv("../../data/LUT_samples.csv")

# read in the list of objects. use the filtered dataset for the singlets only
data.list <- readRDS("../../out/object/list_datasc_fix_filter_norm_doubletSinglet_SoupX_00500_07000_05.rds")

list_plot <- lapply(names(data.list),function(x){
  data <- data.list[[x]]
  plot <- DimPlot(data,raster = T)+ggtitle(x)
  return(plot)
})

wrap_plots(list_plot)

lapply(data.list, function(x){
  gene <- dim(x)[1]
  cell <- dim(x)[2]
  
  data.frame(gene,cell)
}) %>%
  bind_rows(.id = "sample") %>%
  summarise(tot = sum(cell))

read_tsv("../../out/table/meta_datasc_beforeQC.tsv") %>%
  summarise(n = n(),
            tot_UMI = sum(nCount_RNA),
            avg_UMI = mean(nCount_RNA),
            median_UMI = median(nCount_RNA),
            avg_mito_prt = mean(percent.mt),
            median_mito_prt = median(percent.mt))

lapply(data.list, function(x){
  x@meta.data
}) %>%
  bind_rows(.id = "sample") %>%
  summarise(n = n(),
          tot_UMI = sum(nCount_RNA),
          avg_UMI = mean(nCount_RNA),
          median_UMI = median(nCount_RNA),
          avg_mito_prt = mean(percent.mt),
          median_mito_prt = median(percent.mt))


read_tsv("../../out/table/meta_datasc_beforeQC.tsv")
