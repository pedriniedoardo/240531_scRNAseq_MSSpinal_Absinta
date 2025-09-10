# AIM ---------------------------------------------------------------------
# run GSEA on the table of DE from pbulk analysis

# library -----------------------------------------------------------------
library(tidyverse)
library(fgsea)
# install.packages('msigdbdf', repos = 'https://igordot.r-universe.dev')
# library(msigdbdf)
library(msigdbr)
library(GSEABase)
library(patchwork)

# prepare the dataset with all the annoration needed ---------------------- 
file <- "30_DE_pseudobulk_filterExp_shr_IMMUNE.tsv"

# load the results 
results <- lapply(paste0("../../out/table/analysis_R44/",file),function(x){
  read_tsv(x) %>%
    filter(cluster == "IMMUNE")
}) %>%
  setNames(str_remove_all(file,pattern = ".tsv"))

# GSEA -------------------------------------------------------------------- 
# use the FC dataset to create the ranked list of genes 
# Symbol or Entrez? 
# x <- results$res_MutvsWT_shr
list_ranks <- lapply(results, function(x){
  
  x <- dplyr::filter(x,!is.na(symbol)) %>%
    dplyr::filter(!is.na(pvalue)) %>%
    group_by(symbol) %>%
    # average pvalue in case of duplicated genenames
    summarise(pvalue = mean(pvalue),
              log2FC = mean(log2FoldChange)) %>%
    ungroup() %>%
    # arrange(pvalue) %>%
    # change the inf to big numbers
    # mutate(signed_rank_stats = sign(log2FC) * -log10(pvalue)) %>%
    mutate(negative_log10pvalue = -log10(pvalue)) %>%
    # multiply the 1000 by the log2FC to break the ranking ties
    # mutate(negative_log10pvalue = ifelse(is.infinite(negative_log10pvalue), 1000 * log2FC, negative_log10pvalue)) %>%
    mutate(negative_log10pvalue = ifelse(is.infinite(negative_log10pvalue), 1000, negative_log10pvalue)) %>%
    mutate(signed_rank_stats = sign(log2FC) * negative_log10pvalue)
  
  
  ranks <- setNames(x$signed_rank_stats, x$symbol)
  ranks
}) 
glimpse(list_ranks)

# define the signature object ---------------------------------------------
# library("msigdbr")
msigdbr_collections() %>% print(n = 30)

# define a list of parameters to run a set of annotations
list_annotation <- read_csv("../../data/LUT_GSEA.csv") %>%
  split(f = .$annotation)

# loop the processing
list_gene_sets_01 <- lapply(list_annotation,function(x){
  # track the progress
  print(x$annotation)
  
  # load the annotation
  gene_sets <- msigdbr(species = "Homo sapiens", collection = x$collection,subcollection = x$subcollection)
  
  return(gene_sets)
})

list_gene_sets_02 <- list("HALLMARK" = msigdbr(species = "Homo sapiens", collection = "H"))

list_gene_sets <- c(list_gene_sets_01,
                    list_gene_sets_02)


# format in order to be accepted by GSEA

list_pathways <- lapply(list_gene_sets, function(x){
  split(x = x$gene_symbol,
        f = x$gs_name)
})
# head(pathways)

# RUN GSEA ----------------------------------------------------------------
# loop the process in the pathway list
# name_anno <- names(list_pathways[1])
# pathways <- list_pathways[[1]]

pmap(list(names(list_pathways),list_pathways),function(name_anno,pathways){
  # track the progress
  print(name_anno)
  
  # run the enrichment
  # add a seed to fix the GSEA result
  set.seed(123)
  
  df_tables_GSEA_all <- lapply(list_ranks, function(x){
    fgsea(pathways, x, minSize=10, maxSize=500)  
  }) %>%
    bind_rows(.id = "dataset") %>% 
    # the ladingEdge columns has to be re-arranged in order to save the file as a table (originally is a list) 
    mutate(leadingEdge = unlist(lapply(.$leadingEdge, function(x){
      paste0(x,collapse = "|")
    }))) %>%
    arrange(padj,-abs(NES)) 
  
  dim(df_tables_GSEA_all)
  
  head(df_tables_GSEA_all,n=20) 
  
  # save the whole table
  df_tables_GSEA_all %>%
    write_tsv(paste0("../../out/table/analysis_R44/33_df_table_GSEA_",name_anno,"_rankSigLogPval_","IMMUNE",".tsv"))
  
  # COLLAPSE REDUNDANT ------------------------------------------------------
  # collapsing the similar pathways 
  
  # split the dataset per type
  list_tables_GSEA_all <- split(df_tables_GSEA_all,f = df_tables_GSEA_all$dataset)
  
  names(list_tables_GSEA_all)
  
  list_collapsedPathways <- lapply(names(list_tables_GSEA_all),function(x){
    collapsePathways(list_tables_GSEA_all[[x]], pathways, list_ranks[[x]])
  }) %>%
    setNames(names(list_tables_GSEA_all))
  
  str(list_collapsedPathways)
  
  list_mainPathways <- pmap(list(list_tables_GSEA_all,list_collapsedPathways),function(x,y){
    x %>%
      dplyr::filter(pathway %in% y$mainPathways) %>%
      arrange(padj,-abs(NES)) %>%
      pull(pathway) 
  })
  
  str(list_mainPathways)
  
  # save list of non redundant terms
  # chackt the order of the names is the same
  sum(!names(list_tables_GSEA_all) == names(list_mainPathways))
  
  # filter only the non redundant fro each comparison
  df_tables_GSEA_all_non_redundant <- 
    pmap(list(list_tables_GSEA_all,list_mainPathways),function(x,y){
      x %>%
        dplyr::filter(pathway %in% y)
    }) %>%
    bind_rows()
  
  # save the table
  df_tables_GSEA_all_non_redundant %>%
    write_tsv(file = paste0("../../out/table/analysis_R44/33_df_table_GSEA_",name_anno,"_nonredundant_rankSigLogPval_","IMMUNE",".tsv"))
  
  test <- df_tables_GSEA_all_non_redundant %>%
    group_by(dataset) %>%
    top_n(wt = padj*(-1),n = 5)
  
  # test plot to show the main terms in each dataset
  library(ggrepel)
  df_tables_GSEA_all_non_redundant %>%
    # shorten the label of the pathway
    mutate(pathway2 = str_remove(pathway,pattern = paste0(name_anno,"_")) %>%
             str_sub(start = 1,end = 35)) %>%
    # mutate(min_log10_padj = -log10(padj)) %>%
    ggplot(aes(y = -log(padj),x = NES,label = pathway2)) + geom_point(aes(size = size),alpha = 0.2) + facet_wrap(~dataset) + theme_bw() +
    theme(strip.background = element_blank())+
    geom_text_repel(size = 2,box.padding = 0.5,segment.alpha = 0.6,max.overlaps = 10)+
    geom_hline(yintercept = -log(0.05),col="gray",linetype="dashed")
  ggsave(paste0("../../out/plot/analysis_R44/33_GSEA_unbiased_",name_anno,"_nonredundant_rankSigLogPval_","IMMUNE",".pdf"),width = 15,height = 10)
  
  # library(ggrepel)
  df_tables_GSEA_all %>%
    # shorten the label of the pathway
    mutate(pathway2 = str_remove(pathway,pattern = paste0(name_anno,"_")) %>%
             str_sub(start = 1,end = 35)) %>%
    # mutate(min_log10_padj = -log10(padj)) %>%
    ggplot(aes(y = -log(padj),x = NES,label = pathway2)) + geom_point(aes(size = size),alpha = 0.2) + facet_wrap(~dataset) + theme_bw() +
    theme(strip.background = element_blank())+
    geom_text_repel(size = 2,box.padding = 0.5,segment.alpha = 0.6,max.overlaps = 10)+
    geom_hline(yintercept = -log(0.05),col="gray",linetype="dashed")
  ggsave(paste0("../../out/plot/analysis_R44/33_GSEA_unbiased_",name_anno,"_rankSigLogPval_","IMMUNE",".pdf"),width = 15,height = 9)
})

# merge all the individual tables in one
df_GSEA_rankPvalue_all <- lapply(names(list_pathways),function(x){
  # check the progress
  print(x)
  
  # read in the table
  df_tables_GSEA <- read_tsv(paste0("../../out/table/analysis_R44/33_df_table_GSEA_",x,"_rankSigLogPval_","IMMUNE",".tsv")) %>%
    mutate(annotation = x)
  
  return(df_tables_GSEA)
}) %>%
  bind_rows()

df_GSEA_rankPvalue_all %>%
  write_tsv("../../out/table/analysis_R44/33_df_table_GSEA_all_rankSigLogPval_IMMUNE.tsv")

# compare the two results from the two ranking systems --------------------

# read all the terms from the different annotaitons
df_GSEA_rankPvalue <- lapply(names(list_pathways),function(x){
  # check the progress
  print(x)
  
  # read in the table
  df_tables_GSEA_rankPvalue <- read_tsv(paste0("../../out/table/analysis_R44/33_df_table_GSEA_",x,"_rankSigLogPval_","IMMUNE",".tsv")) %>%
    dplyr::select(dataset,pathway,NES,padj) %>%
    mutate(annotation = x)
  
  return(df_tables_GSEA_rankPvalue)
}) %>%
  bind_rows()

df_GSEA_rankLogFC <- lapply(names(list_pathways),function(x){
  # check the progress
  print(x)
  
  # read in the table
  df_tables_GSEA_rankLogFC <- read_tsv(paste0("../../out/table/analysis_R44/33_df_table_GSEA_",x,"_ranklogfc_","IMMUNE",".tsv")) %>%
    dplyr::select(dataset,pathway,NES,padj) %>%
    mutate(annotation = x)
  
  return(df_tables_GSEA_rankLogFC)
}) %>%
  bind_rows()

# join the two tables
test <-  df_GSEA_rankPvalue %>%
  left_join(df_GSEA_rankLogFC,by = c("annotation","dataset","pathway"),suffix = c(".pvalue",".logfc"))

# plot the comparison for the two ranking systems
test %>%  
  ggplot(aes(x=NES.pvalue, y = NES.logfc)) + geom_point(shape = 1) + theme_bw() + facet_wrap(~annotation) + theme(strip.background = element_blank())

# 

# arrange the plot for the top terms in GSEA
# anno<-"CP"
list_plot <- lapply(names(list_pathways),function(anno){
  # check the progress
  print(anno)
  
  # pull the range of the NES to scale the dots
  global_min <- min(-log(c(test %>%
                             filter(annotation %in% anno) %>%
                             slice_max(abs(NES.logfc),n = 10) %>%
                             pull(padj.logfc),
                           test %>%
                             filter(annotation %in% anno) %>%
                             slice_max(abs(NES.pvalue),n = 10) %>%
                             pull(padj.pvalue))))
  global_max <- max(-log(c(test %>%
                             filter(annotation %in% anno) %>%
                             slice_max(abs(NES.logfc),n = 10) %>%
                             pull(padj.logfc),
                           test %>%
                             filter(annotation %in% anno) %>%
                             slice_max(abs(NES.pvalue),n = 10) %>%
                             pull(padj.pvalue))))
  
  p1 <- test %>%
    filter(annotation %in% anno) %>%
    slice_max(abs(NES.pvalue),n = 10) %>%
    # mutate(Term = str_remove_all(pathway,pattern = "KEGG_MEDICUS_") %>% str_sub(start = 1,end = 50)) %>%
    mutate(Term = str_remove_all(pathway,pattern = paste0(anno,"_")) %>% str_sub(start = 1,end = 50)) %>%
    mutate(Term = fct_reorder(Term, NES.pvalue,.desc = F)) %>%
    mutate(direction = factor(sign(NES.pvalue),labels = c(1,-1),levels = c(1,-1))) %>%
    ggplot(aes(x = NES.pvalue,y=Term)) + 
    geom_point(aes(size = -log(padj.pvalue),col=direction)) +
    theme_bw() +
    geom_vline(xintercept = 0,col="gray",linetype = "dashed") +
    scale_color_manual(values = c("red","blue")) +
    scale_size_continuous(limits = c(global_min, global_max)) +
    coord_cartesian(xlim = c(-3,3)) + 
    ggtitle(paste0("GSEA ", anno," signed NegLogPval ranking"))
  
  p2 <- test %>%
    filter(annotation %in% anno) %>%
    slice_max(abs(NES.logfc),n = 10) %>%
    # mutate(Term = str_remove_all(pathway,pattern = "KEGG_MEDICUS_") %>% str_sub(start = 1,end = 50)) %>%
    mutate(Term = str_remove_all(pathway,pattern = paste0(anno,"_")) %>% str_sub(start = 1,end = 50)) %>%
    mutate(Term = fct_reorder(Term, NES.logfc,.desc = F)) %>%
    mutate(direction = factor(sign(NES.logfc),labels = c(1,-1),levels = c(1,-1))) %>%
    ggplot(aes(x = NES.logfc,y=Term)) + 
    geom_point(aes(size = -log(padj.logfc),col=direction)) +
    theme_bw() +
    geom_vline(xintercept = 0,col="gray",linetype = "dashed") +
    scale_color_manual(values = c("red","blue")) +
    scale_size_continuous(limits = c(global_min, global_max)) +
    coord_cartesian(xlim = c(-3,3)) + 
    ggtitle(paste0("GSEA ", anno," logFC ranking"))
  
  # arrange the plot
  # p_tot <- (p1 / p2) + plot_annotation(
  #   title = paste0("GSEA_", anno), 
  #   theme = theme(plot.title = element_text(hjust = 0.5))
  # )
  
  p_tot <- (p1 / p2)
  
  return(p_tot)
})

#
wrap_plots(list_plot,nrow = 1)
ggsave("../../out/plot/analysis_R44/33_GSEA_top10_all.pdf",width = 40,height = 10)

# focus on interferon response
read_tsv(paste0("../../out/table/analysis_R44/33_df_table_GSEA_HALLMARK_rankSigLogPval_IMMUNE.tsv")) %>%
  filter(str_detect(pathway,pattern = "INTERFERON_ALPHA_RESPONSE")) %>%
  pull(leadingEdge) %>%
  str_split(pattern = "\\|") %>%
  unlist() %>%
  saveRDS("../../out/object/analysis_R44/33_leadingEges_HALLMARK_INTERFERON_ALPHA_RESPONSE_rankSigLogPval.rds")
  
# plot volcano like plot for all the annoations together

test %>%
  # shorten the label of the pathway
  mutate(pathway2 = str_remove(pathway,pattern = paste0(paste0(names(list_pathways),collapse = "|"),"_")) %>%
           str_sub(start = 1,end = 50)) %>%
  # mutate(min_log10_padj = -log10(padj)) %>%
  ggplot(aes(y = -log(padj.pvalue),x = NES.pvalue,label = pathway2)) + geom_point(alpha = 0.2) + facet_wrap(~annotation) + theme_bw() +
  theme(strip.background = element_blank())+
  geom_text_repel(size = 2,box.padding = 0.5,segment.alpha = 0.6,max.overlaps = 10)+
  geom_hline(yintercept = -log(0.05),col="gray",linetype="dashed")
ggsave("../../out/plot/analysis_R44/33_GSEA_volcano_rankSigLogPval_all.pdf",width = 30,height = 20)

test %>%
  # shorten the label of the pathway
  mutate(pathway2 = str_remove(pathway,pattern = paste0(paste0(names(list_pathways),collapse = "|"),"_")) %>%
           str_sub(start = 1,end = 50)) %>%
  # mutate(min_log10_padj = -log10(padj)) %>%
  ggplot(aes(y = -log(padj.logfc),x = NES.logfc,label = pathway2)) + geom_point(alpha = 0.2) + facet_wrap(~annotation) + theme_bw() +
  theme(strip.background = element_blank())+
  geom_text_repel(size = 2,box.padding = 0.5,segment.alpha = 0.6,max.overlaps = 10)+
  geom_hline(yintercept = -log(0.05),col="gray",linetype="dashed")
ggsave("../../out/plot/analysis_R44/33_GSEA_volcano_rankLogFC_all.pdf",width = 30,height = 20)