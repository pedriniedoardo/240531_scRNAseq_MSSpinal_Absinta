# AIM ---------------------------------------------------------------------
# sample routine for the pseudobulk analysis.
# the following is a readaptation of different resources
# https://github.com/sib-swiss/single-cell-training/
# https://hbctraining.github.io/scRNA-seq/lessons/pseudobulk_DESeq2_scrnaseq.html

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(GGally)
library(cowplot)
library(ComplexHeatmap)
library(scales)
library(circlize)
library(DESeq2)
library(RNAseqQC)
library(limma)
library(ashr)
library(magick)
library(UpSetR)

# read in the data --------------------------------------------------------
# read in the sample dataset
data.combined <- readRDS("../../out/object/analysis_R44/26_sobj_integrated_cleanup_manualAnnotation.rds")
DimPlot(data.combined, group.by = "cell_id")

# For this test I would focus on the subset only the presumed cells of interest cells for the test
# in this case use all the cells
scobj_subset <- subset(data.combined,subset = cell_id %in% c("IMMUNE"))

# scobj_subset <- data.combined
# DimPlot(scobj_subset_subset,label = T,raster = T,group.by = "stim")

# summarise the metadata per sample
meta_summary <- scobj_subset@meta.data %>%
  group_by(sample_id,age,weight,csf,apoe,pathological_stage_short,demyelination_short) %>%
  summarise(avg_percent.mt = mean(percent.mt), avg_percent.ribo = mean(percent.ribo))

# wrangling ---------------------------------------------------------------
# explore the full meta
scobj_subset@meta.data %>%
  group_by(sample_id,diagnosis_short) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = diagnosis_short,values_from = n)

# aggregate the expression per sample per donor.id and stimulation
cts_sample_all <- AggregateExpression(object = scobj_subset,
                                      group.by = c("sample_id","sex","diagnosis_short","location"),
                                      assays = 'RNA',
                                      slot = "counts",
                                      return.seurat = FALSE)

# processing --------------------------------------------------------------
# extract
# 1. Get counts matrix
counts <- cts_sample_all$RNA %>%
  as.data.frame()

# build a metadata to add to the deseq2 object to track the total number of cells of reads per sample
meta_reads <- colSums(counts) %>% data.frame(tot_reads = .) %>% mutate(M_reads = tot_reads/10^6) %>% rownames_to_column("sample") %>%
  separate(sample,into = c("sample_id","sex","diagnosis_short","location"),sep = "_") %>%
  dplyr::select(-c(tot_reads,sex,diagnosis_short,location)) %>%
  mutate(sample_id = str_replace_all(sample_id,pattern = "-","_"))

# 2. generate sample level metadata
LUT_sample <- scobj_subset@meta.data %>%
  group_by(sample_id,sex,diagnosis_short,location) %>%
  summarise() %>%
  # match the sample_id with the one in the LUT
  mutate(sample_id_fix = str_replace_all(sample_id,pattern = "_",replacement = "-")) %>%
  # harmonyze the naming of the donor_id
  # mutate(sample_disease = paste(sample_id_fix,diagnosis_short,sep = "_")) %>%
  mutate(sample_disease_origin = paste(sample_id_fix,sex,diagnosis_short,location,sep = "_")) %>%
  ungroup() %>%
  # add the summary info
  left_join(meta_summary,by = "sample_id") %>%
  # add total number of reads per sample
  left_join(meta_reads,by = "sample_id")

# match the LUT with the expression matrix
colData <- data.frame(sample.id = colnames(counts)) %>%
  # separate(samples,into = c(c("orig.ident","origin_01","pathology_class","origin_02")),remove = F,sep = "_") %>%
  left_join(LUT_sample,by = c("sample.id" = "sample_disease_origin"))

# save matrix and metadata
saveRDS(counts,file = paste0("../../out/object/analysis_R44/30_counts_mono_pBulk_","IMMUNE",".rds"))
saveRDS(colData,file = paste0("../../out/object/analysis_R44/30_colData_mono_pBulk_","IMMUNE",".rds"))

# perform DESeq2 ----------------------------------------------------------
# build the model
treat <- colData$diagnosis
# gender <- colData$sex
location <- colData$location

# design <- model.matrix(~ block_donor + treat)
# design <- model.matrix(~ gender + location + treat)
design <- model.matrix(~ location + treat)
colnames(design)[1] <- c("intercept")

# save the disign
saveRDS(design,paste0("../../out/object/analysis_R44/30_design_mono_pBulk","IMMUNE",".rds"))

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = design)

# filter low aboundant features
feat_keep <- edgeR::filterByExpr(counts(dds), group = colData$disease)
dds_filter <- dds[feat_keep,]

# plot the raw number of reads per sample
colSums(counts(dds_filter)) %>%
  data.frame(tot_reads = .) %>%
  rownames_to_column("sample") %>%
  ggplot(aes(x=sample,y=tot_reads)) + geom_col()+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))
ggsave(paste0("../../out/plot/analysis_R44/30_pBulk_totReads_","IMMUNE",".pdf"),width = 10,height = 6)

# scale the data
vds_filter <- vst(dds_filter, blind = T)

# clustering samples ------------------------------------------------------
# set seed to control random annotation colors
pdf(paste0("../../out/plot/analysis_R44/30_heatmap_SampleCluster_filterExp_","IMMUNE",".pdf"),width = 12,height = 6)
set.seed(21)
hm <- plot_sample_clustering(vds_filter,
                             anno_vars = c("diagnosis_short","location","sex","pathological_stage_short","demyelination_short","age","M_reads"),
                             distance = "euclidean")
draw(hm,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(2, 2, 2, 30), "mm"))
dev.off()

# PCA ---------------------------------------------------------------------
plot_vsd <- plotPCA(vds_filter,
                    intgroup = c("diagnosis_short","sex","location","sample_id_fix","age","weight","csf","apoe","pathological_stage_short","demyelination_short","avg_percent.mt","avg_percent.ribo","M_reads")) +
  theme_bw()

plot_vsd$data %>%
  separate(name,into = c("sample","sex","diagnosis","location"),sep = "_",remove = F) %>%
  ggplot(aes(x=PC1,y=PC2,col=location)) +
  geom_point(size =3,alpha=0.6) +
  # ggrepel::geom_text_repel(show.legend = F)+
  scale_x_continuous(expand = expansion(mult = 0.1))+
  theme_bw() + ylab(plot_vsd$labels[1]) + xlab(plot_vsd$labels[2])
ggsave(paste0("../../out/plot/analysis_R44/30_PCA_pseudobulk_filterExp_","IMMUNE",".pdf"),width = 12,height = 10)

plot_vsd$data %>%
  separate(name,into = c("sample","sex","diagnosis","location"),sep = "_",remove = F) %>%
  ggplot(aes(x=PC1,y=PC2,col=diagnosis,label = sample)) +
  geom_point(size =3,alpha=0.6) +
  ggrepel::geom_text_repel(show.legend = F)+
  scale_x_continuous(expand = expansion(mult = 0.1))+
  theme_bw() + ylab(plot_vsd$labels[1]) + xlab(plot_vsd$labels[2])
ggsave(paste0("../../out/plot/analysis_R44/30_PCA_pseudobulk_filterExp2_","IMMUNE",".pdf"),width = 12,height = 10)

# pull more PC
rv <- rowVars(assay(vds_filter),useNames = TRUE)

# select the ntop genes by variance
select_var <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]

# perform a PCA on the data in assay(x) for the selected genes
test <- prcomp(t(assay(vds_filter)[select_var,]))$x %>% 
  data.frame() %>% 
  rownames_to_column("sample")

# plot more PC by condition
left_join(plot_vsd$data %>% dplyr::select(-c("PC1","PC2")),test,by = c("name"="sample")) %>%
  # ggpairs(columns = 5:14,ggplot2::aes(colour=condition),upper = "blank")+
  ggpairs(columns = 16:40,ggplot2::aes(colour=treat),upper = "blank") +
  theme_bw()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave(paste0("../../out/plot/analysis_R44/30_PCA_pseudobulk_panel_","IMMUNE",".pdf"),width = 20,height = 20)

# explore pc score by metadata for the samples cat
test2 <- left_join(plot_vsd$data %>% dplyr::select(-c("PC1","PC2")),test,by = c("name"="sample"))

test_df1 <- test2 |> 
  dplyr::select(name,diagnosis_short,sex,location,pathological_stage_short,demyelination_short) |> 
  # dplyr::rename(approx.time = Approx.Time.between.collection.and.processing) |> 
  pivot_longer(names_to = "var_1",values_to = "value_1",-c(name))

test_df2 <- test2 |> 
  dplyr::select(name,PC1:PC10) |>
  pivot_longer(names_to = "var_2",values_to = "value_2",-c(name))

left_join(test_df1,test_df2,by=c("name")) |> 
  mutate(var_2 = str_remove(var_2,"PC") %>% str_pad(side = "left",width = 2,pad = "0") %>% paste0("PC",.)) |> 
  mutate(comparison = paste0(var_1,"_vs_",var_2)) |> 
  ggplot(aes(x=value_1,y=value_2)) +
  facet_wrap(~comparison,scales = "free",ncol=10) +
  # facet_grid(var_1~var_2,scales = "free_x",drop = T) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.2))+
  theme_bw()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave(paste0("../../out/plot/analysis_R44/30_panel_metadata_PC_pbulk_cat_","IMMUNE",".pdf"),width = 30,height = 15)

# continuos variables
# test2 <- left_join(plot_vsd$data %>% dplyr::select(-c("PC1","PC2")),test,by = c("name"="sample"))

test_df3 <- test2 |> 
  dplyr::select(name,age,weight,csf,avg_percent.mt,avg_percent.ribo,M_reads) |> 
  # dplyr::rename(approx.time = Approx.Time.between.collection.and.processing) |> 
  pivot_longer(names_to = "var_1",values_to = "value_1",-c(name))

left_join(test_df3,test_df2,by=c("name")) |> 
  mutate(var_2 = str_remove(var_2,"PC") %>% str_pad(side = "left",width = 2,pad = "0") %>% paste0("PC",.)) |> 
  mutate(comparison = paste0(var_1,"_vs_",var_2)) |> 
  ggplot(aes(x=value_1,y=value_2)) +
  facet_wrap(~comparison,scales = "free",ncol=10) +
  # facet_grid(var_1~var_2,scales = "free_x",drop = T) +
  # geom_boxplot(outlier.shape = NA)+
  geom_smooth(method = "lm") +
  geom_point(position = position_jitter(width = 0.2))+
  theme_bw()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave(paste0("../../out/plot/analysis_R44/30_panel_metadata_PC_pbulk_cont_","IMMUNE",".pdf"),width = 30,height = 15)

# cell composition --------------------------------------------------------
# load the full dataset
# scobj <- readRDS("../../out/object/129_MG_subcluster_HarmonySample_martinaCluster.rds")
DimPlot(data.combined,group.by = "cell_id")

# check the relative proportion of different cell types from different sample
sample_prop_wide <- data.combined@meta.data %>%
  group_by(sample_id,sex,location,diagnosis_short,cell_id) %>%
  summarise(n = n(),.groups = "drop") %>%
  # match the sample_id with the one in the LUT
  mutate(sample_id_fix = str_replace_all(sample_id,pattern = "_",replacement = "-")) %>%
  mutate(sample_disease_location = paste(sample_id_fix,sex,diagnosis_short,location,sep = "_")) %>%
  # ensure all the combinations are presents
  complete(cell_id,sample_disease_location, fill = list(n = 0)) %>%
  group_by(sample_disease_location) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  # make it wider
  dplyr::select(cell_id,sample_disease_location,prop) %>%
  pivot_wider(names_from = sample_disease_location,values_from = prop) %>%
  # pivot_wider(names_from = stim_donor,values_from = prop,values_fill = 0) %>%
  column_to_rownames("cell_id") %>%
  # select only the samples in the dataset
  dplyr::select(colData$sample.id)

# confim the dimensions
colSums(sample_prop_wide)

# plot the data as heatmap
meta_sample_prop <- data.frame(colname = colnames(sample_prop_wide)) %>%
  left_join(colData,by=c("colname"="sample.id"))

column_meta_sample_prop <- HeatmapAnnotation(diagnosis = meta_sample_prop$diagnosis_short,
                                             location = meta_sample_prop$location,
                                             col = list(diagnosis = c("CTRL" = "blue",
                                                                    "MS" = "red"),
                                                        location = c("cervical" = "gray",
                                                                     "thoracic" = "gray30",
                                                                     "lumbar" = "black")))

ht2_shr <- Heatmap(sample_prop_wide, show_column_names = T,raster_by_magick = T,show_row_dend = F, use_raster = T,
                   name = "prop_cell_type",
                   column_title = "sample",
                   col = viridis::viridis(option = "turbo",n = 10),
                   
                   # row_names_gp = gpar(fontsize = 3),
                   top_annotation = column_meta_sample_prop,show_row_names = T
                   # cluster_rows = F,
                   # right_annotation = row_ha,
                   # row_split = rep(c(1,2,3,4),c(2,3,4,7))
)

pdf(paste0("../../out/plot/analysis_R44/30_heatmap_cellType_pseudobulk_","IMMUNE",".pdf"),width = 10,height = 8)
draw(ht2_shr,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(30,2,2, 2), "mm"))
dev.off()

# scale the proportions by cell type
sample_prop_wide2 <- data.combined@meta.data %>%
  group_by(sample_id,sex,location,diagnosis_short,cell_id) %>%
  summarise(n = n(),.groups = "drop") %>%
  # match the sample_id with the one in the LUT
  mutate(sample_id_fix = str_replace_all(sample_id,pattern = "_",replacement = "-")) %>%
  mutate(sample_disease_location = paste(sample_id_fix,sex,diagnosis_short,location,sep = "_")) %>%
  # ensure all the combinations are presents
  complete(cell_id,sample_disease_location, fill = list(n = 0)) %>%
  group_by(sample_disease_location) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  ungroup() %>%
  group_by(cell_id) %>%
  mutate(scale_prop = scale(prop) %>% .[,1]) %>%
  # make it wider
  dplyr::select(cell_id,sample_disease_location,scale_prop) %>%
  pivot_wider(names_from = sample_disease_location,values_from = scale_prop) %>%
  # pivot_wider(names_from = stim_donor,values_from = prop,values_fill = 0) %>%
  column_to_rownames("cell_id") %>%
  # select only the samples in the dataset
  dplyr::select(colData$sample.id)

# confirm the scaling
apply(sample_prop_wide2,1,mean)
apply(sample_prop_wide2,1,sd)

# plot the data as heatmap
ht3_shr <- Heatmap(sample_prop_wide2, show_column_names = T,raster_by_magick = T,show_row_dend = F, use_raster = T,
                   name = "prop_cell_type",
                   column_title = "sample",
                   col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")),
                   
                   # row_names_gp = gpar(fontsize = 3),
                   top_annotation = column_meta_sample_prop,show_row_names = T
                   # cluster_rows = F,
                   # right_annotation = row_ha,
                   # row_split = rep(c(1,2,3,4),c(2,3,4,7))
)

pdf(paste0("../../out/plot/analysis_R44/30_heatmap_cellType_pseudobulk2_","IMMUNE",".pdf"),width = 10,height = 8)
draw(ht3_shr,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(30,2,2, 2), "mm"))
dev.off()

# run DE ------------------------------------------------------------------
# run DESeq2
ddsHTSeq_filter <- DESeq(dds_filter)

# Check the coefficients for the comparison
resultsNames(ddsHTSeq_filter)

# save the filtered object
saveRDS(ddsHTSeq_filter,paste0("../../out/object/analysis_R44/30_ddsHTSeq_pseudobulk_filterExp_","IMMUNE",".rds"))

# print the contrast
resultsNames(ddsHTSeq_filter)

# save the resut of the contrast of interest. notice that is also possible to set the alpha value for the significance (adj-p value)
contrast <- makeContrasts(MS_vs_CTRL = treatMS,
                          levels = design)

# pull the results table
res_MS <- results(ddsHTSeq_filter, contrast=contrast[,"MS_vs_CTRL"],alpha = 0.05)
summary(res_MS)

# add the gene symbols
df_res <-
  list(res_MS = res_MS) %>%
  lapply(function(x){
    df <- x %>%
      data.frame() %>%
      rownames_to_column(var = "symbol") %>%
      arrange(pvalue)
    # left_join(LUT_gene,by = c("symbol"="gene_name"))
    df
  }) %>%
  bind_rows(.id = "condVsCTRL")

# save the table
df_res %>%
  mutate(cluster = "IMMUNE") %>%
  write_tsv(paste0("../../out/table/analysis_R44/30_DE_pseudobulk_filterExp_","IMMUNE",".tsv"))

# check some of the genes that were significant in the sc DE analysis
df_res %>%
  filter(symbol %in% c("IFITM3", "BST2", "MT2A", "PLSCR1", "BTG1", "MX2", "MX1", "FPR1", "CASP1"))

# shrink ------------------------------------------------------------------
res_MS_shr <- lfcShrink(ddsHTSeq_filter, res = res_MS, type = "ashr")

# save the table of DEGs
df_res_shr <-
  list(res_MS_shr = res_MS_shr) %>%
  lapply(function(x){
    df <- x %>%
      data.frame() %>%
      rownames_to_column(var = "symbol") %>%
      arrange(pvalue)
    # left_join(LUT_gene,by = c("symbol"="gene_name"))
    df
  }) %>%
  bind_rows(.id = "condVsCTRL")

# save the table
df_res_shr %>%
  mutate(cluster = "IMMUNE") %>%
  write_tsv(paste0("../../out/table/analysis_R44/30_DE_pseudobulk_filterExp_shr_","IMMUNE",".tsv"))


# Another useful diagnostic plot is the histogram of the p values (figure below). This plot is best formed by excluding genes with very small counts, which otherwise generate spikes in the histogram.
df_res %>%
  data.frame()%>%
  dplyr::filter(baseMean>1) %>%
  ggplot(aes(x=pvalue))+geom_histogram(breaks = 0:20/20) +
  facet_wrap(~condVsCTRL)+
  theme_bw()+
  theme(strip.background = element_blank())
ggsave(paste0("../../out/plot/analysis_R44/30_histogram_pvalue_pseudobulk_filterExp_","IMMUNE",".pdf"),width = 6,height = 4)

# volcano -----------------------------------------------------------------
# add the info of the genename
plot_volcano <- df_res %>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>0.5&!is.na(symbol)),yes = 1,no = 0)) %>%
  filter(!is.na(col))

plot_volcano %>%
  ggplot(aes(x=log2FoldChange,y=-log(padj)))+
  # geom_point()
  geom_point(data = plot_volcano[plot_volcano$col==0,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.05)+
  geom_point(data = plot_volcano[plot_volcano$col==1,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.5)+
  geom_vline(xintercept = c(-0.5,0.5),col="red",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="red",linetype="dashed")+
  scale_color_manual(values = c("black","red"))+theme(legend.position = "none")+
  # ggrepel::geom_text_repel(
  #   data = plot_volcano[plot_volcano$col==1,][1:1000,],
  #   aes(label = symbol),max.overlaps = 1,segment.alpha=0.4,
  #   size = 2,
  #   box.padding = unit(0.35, "lines"),
  #   point.padding = unit(0.3, "lines")) +
  ggrepel::geom_text_repel(
    data = plot_volcano %>% 
      group_by(condVsCTRL) %>% 
      arrange(padj) %>% 
      dplyr::slice(1:30) %>% 
      dplyr::filter(abs(log2FoldChange)>0.5) %>% 
      dplyr::filter(padj<0.05),
    aes(label = symbol),segment.alpha=0.4) +
  facet_wrap(~condVsCTRL)+
  theme_bw()+
  theme(strip.background = element_blank())+
  theme(legend.position = "none")
ggsave(paste0("../../out/plot/analysis_R44/30_vulcano_plot_pseudobulk_filterExp_","IMMUNE",".pdf"),width = 10,height = 10)

#
plot_volcano_shr <- df_res_shr %>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>0.5&!is.na(symbol)),yes = 1,no = 0)) %>%
  filter(!is.na(col))

plot_volcano_shr %>%
  ggplot(aes(x=log2FoldChange,y=-log(padj)))+
  # geom_point()
  geom_point(data = plot_volcano_shr[plot_volcano_shr$col==0,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.05)+
  geom_point(data = plot_volcano_shr[plot_volcano_shr$col==1,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.5)+
  geom_vline(xintercept = c(-0.5,0.5),col="red",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="red",linetype="dashed")+
  scale_color_manual(values = c("black","red"))+theme(legend.position = "none")+
  ggrepel::geom_text_repel(
    data = plot_volcano_shr %>% group_by(condVsCTRL) %>% arrange(padj) %>% dplyr::slice(1:30)%>% filter(abs(log2FoldChange)>0.5),
    aes(label = symbol),segment.alpha=0.4) +
  facet_wrap(~condVsCTRL)+
  theme_bw()+
  theme(strip.background = element_blank())+
  theme(legend.position = "none")
ggsave(paste0("../../out/plot/analysis_R44/30_vulcano_plot_pseudobulk_filterExp_shr_","IMMUNE",".pdf"),width = 10,height = 10)

# MA plot -----------------------------------------------------------------
df_res %>%
  filter(!is.na(padj)) %>%
  data.frame()%>%
  mutate(color = ifelse(padj<0.05,yes = 1,no = 0))%>%
  mutate(color = factor(ifelse(is.na(color),yes = 0,no = color)))%>%
  arrange(color) %>%
  ggplot(aes(baseMean,y = log2FoldChange,col=color)) + geom_point(alpha=0.2) +
  scale_x_log10() + scale_color_manual(values = c("gray","red")) +
  facet_wrap(~condVsCTRL)+
  theme_bw()+
  theme(strip.background = element_blank()) +
  theme(legend.position = "none")
# ggsave("../../out/image/02_MA_plot_pseudobulk_filterExp.pdf",width = 10,height = 10)

# shr
df_res_shr %>%
  filter(!is.na(padj)) %>%
  data.frame()%>%
  mutate(color = ifelse(padj<0.05,yes = 1,no = 0))%>%
  mutate(color = factor(ifelse(is.na(color),yes = 0,no = color)))%>%
  arrange(color) %>%
  ggplot(aes(baseMean,y = log2FoldChange,col=color)) + geom_point(alpha=0.2) +
  scale_x_log10() + scale_color_manual(values = c("gray","red")) +
  facet_wrap(~condVsCTRL)+
  theme_bw()+
  theme(strip.background = element_blank()) +
  theme(legend.position = "none")
# ggsave("../../out/image/02_MA_plot_pseudobulk_filterExp_shr.pdf",width = 10,height = 10)

# heatmaps ----------------------------------------------------------------
# pull the scaled values
mat_filter <- assay(vds_filter) %>%
  as.data.frame() %>%
  # dplyr::select(contains(c("_0_","_6_"))) %>%
  as.matrix()

# the DEGs plot
DEG_2 <- df_res_shr %>%
  as.data.frame()%>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>0.5),yes = 1,no = 0)) %>%
  dplyr::filter(col==1) %>%
  pull(symbol) %>%
  unique()

# mat_filter <- assay(vds_filter) %>%
#   data.frame() %>%
#   # dplyr::select(contains(c("_0_","_6_"))) %>%
#   as.matrix()

mat_shr <- mat_filter[rownames(vds_filter) %in% DEG_2, ]
mat2_shr <- (mat_shr - rowMeans(mat_shr))/rowSds(mat_shr,useNames = TRUE)
#
meta_sample <- data.frame(colname = colnames(mat2_shr)) %>%
  left_join(colData,by=c("colname"="sample.id")) %>%
  mutate(location = factor(location,levels = c("cervical","thoracic","lumbar")))

# make the column of the matrix more readable
colnames(mat2_shr) <- meta_sample$colname

column_ha_shr <- HeatmapAnnotation(disease = meta_sample$diagnosis_short,
                                   origin = meta_sample$location,
                                   col = list(disease = c("CTRL" = "blue",
                                                          "MS" = "red"),
                                              origin = c("lumbar" = "gray",
                                                         "thoracic" = "gray30",
                                                         "cervical" = "black")))

ht2_shr <- Heatmap(mat2_shr, show_column_names = T,raster_by_magick = T,show_row_dend = F, use_raster = T,
                   name = "exp",
                   column_title = "test_shr",
                   # row_names_gp = gpar(fontsize = 3),
                   top_annotation = column_ha_shr,show_row_names = T,
                   # cluster_rows = F,
                   # right_annotation = row_ha,
                   # row_split = rep(c(1,2,3,4),c(2,3,4,7))
)
pdf(paste0("../../out/plot/analysis_R44/30_heatmap_DEG_plot_pseudobulk_filterExp_shr_","IMMUNE",".pdf"),width = 10,height = 8)
draw(ht2_shr,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(30,2,2, 2), "mm"))
dev.off()

# # upset plot --------------------------------------------------------------
# # read in the table of DEGs
# df_res_shr <- read_tsv("../../out/table/DE_pseudobulk_ASTRO_ctrl_refCX_shr.tsv")
# 
# # build a list of common elements belonging to each set of fegs
# list_DE_up <- df_res_shr %>%
#   split(f = .$conditionVsCX) %>%
#   lapply(function(x){
#     x %>%
#       filter(padj < 0.05,log2FoldChange>1) %>%
#       pull(symbol) %>%
#       unique()
#   })
# 
# glimpse(list_DE_up)
# 
# list_DE_down <- df_res_shr %>%
#   split(f = .$conditionVsCX) %>%
#   lapply(function(x){
#     x %>%
#       filter(padj < 0.05,log2FoldChange<(-1)) %>%
#       pull(symbol) %>%
#       unique()
#   })
# glimpse(list_DE_down)
# 
# # try the upset plot version
# # library(UpSetR)
# pdf("../../out/image/upset_DEG_UP_plot_pseudobulk_ASTRO_refBASELINE_shr.pdf",width = 14,height = 7)
# upset(fromList(list_DE_up), order.by = "freq",nsets = 7)
# dev.off()
# 
# pdf("../../out/image/upset_DEG_DOWN_plot_pseudobulk_ASTRO_refBASELINE_shr.pdf",width = 14,height = 7)
# upset(fromList(list_DE_down), order.by = "freq",nsets = 7)
# dev.off()
# 
# # pull the intersections
# df1_UP <- lapply(list_DE_up,function(x){
#   data.frame(gene = x)
# }) %>%
#   bind_rows(.id = "path")
# 
# df1_DOWN <- lapply(list_DE_down,function(x){
#   data.frame(gene = x)
# }) %>%
#   bind_rows(.id = "path")
# 
# head(df1_UP)
# head(df1_DOWN)
# 
# df2_UP <- data.frame(gene=unique(unlist(list_DE_up)))
# df2_DOWN <- data.frame(gene=unique(unlist(list_DE_down)))
# 
# head(df2_UP)
# head(df2_DOWN)
# 
# df_int_UP <- lapply(df2_UP$gene,function(x){
#   # pull the name of the intersections
#   intersection <- df1_UP %>%
#     dplyr::filter(gene==x) %>%
#     arrange(path) %>%
#     pull("path") %>%
#     paste0(collapse = "|")
# 
#   # build the dataframe
#   data.frame(gene = x,int = intersection)
# }) %>%
#   bind_rows()
# 
# df_int_DOWN <- lapply(df2_DOWN$gene,function(x){
#   # pull the name of the intersections
#   intersection <- df1_DOWN %>%
#     dplyr::filter(gene==x) %>%
#     arrange(path) %>%
#     pull("path") %>%
#     paste0(collapse = "|")
# 
#   # build the dataframe
#   data.frame(gene = x,int = intersection)
# }) %>%
#   bind_rows()
# 
# df_int_UP %>%
#   write_tsv("../../out/table/upset_DEG_UP_plot_pseudobulk_ASTRO_refBASELINE_shr.tsv")
# 
# df_int_DOWN %>%
#   write_tsv("../../out/table/upset_DEG_DOWN_plot_pseudobulk_ASTRO_refBASELINE_shr.tsv")
# 
# head(df_int_UP,n=20)
# head(df_int_DOWN,n=20)
# 
# df_int_UP %>%
#   group_by(int) %>%
#   summarise(n=n()) %>%
#   arrange(desc(n))
# 
# df_int_DOWN %>%
#   group_by(int) %>%
#   summarise(n=n()) %>%
#   arrange(desc(n))

# PLOT DISPERSION ---------------------------------------------------------
pdf(paste0("../../out/plot/analysis_R44/30_dispersion_plot_pseudobulk_filterExp_","IMMUNE",".pdf"),width = 5,height = 5)
plotDispEsts(ddsHTSeq_filter)
dev.off()

# plot some specific genes ------------------------------------------------


