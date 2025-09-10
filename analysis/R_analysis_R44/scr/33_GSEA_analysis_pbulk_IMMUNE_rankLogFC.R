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
    # average logFC in case of duplicated genenames
    group_by(symbol) %>%
    summarise(logFC = mean(log2FoldChange))
  
  ranks <- setNames(x$logFC, x$symbol)
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
    write_tsv(paste0("../../out/table/analysis_R44/33_df_table_GSEA_",name_anno,"_ranklogfc_","IMMUNE",".tsv"))
  
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
    write_tsv(file = paste0("../../out/table/analysis_R44/33_df_table_GSEA_",name_anno,"_nonredundant_ranklogf_","IMMUNE",".tsv"))
  
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
    ggplot(aes(y = -log10(padj),x = NES,label = pathway2)) + geom_point(aes(size = size),alpha = 0.2) + facet_wrap(~dataset) + theme_bw() +
    theme(strip.background = element_blank())+
    geom_text_repel(size = 2,box.padding = 0.5,segment.alpha = 0.6,max.overlaps = 10)+
    geom_hline(yintercept = -log(0.05),col="gray",linetype="dashed")
  ggsave(paste0("../../out/plot/analysis_R44/33_GSEA_unbiased_",name_anno,"_nonredundant_ranklogfc_","IMMUNE",".pdf"),width = 15,height = 10)
  
  # library(ggrepel)
  df_tables_GSEA_all %>%
    # shorten the label of the pathway
    mutate(pathway2 = str_remove(pathway,pattern = paste0(name_anno,"_")) %>%
             str_sub(start = 1,end = 35)) %>%
    # mutate(min_log10_padj = -log10(padj)) %>%
    ggplot(aes(y = -log10(padj),x = NES,label = pathway2)) + geom_point(aes(size = size),alpha = 0.2) + facet_wrap(~dataset) + theme_bw() +
    theme(strip.background = element_blank())+
    geom_text_repel(size = 2,box.padding = 0.5,segment.alpha = 0.6,max.overlaps = 10)+
    geom_hline(yintercept = -log(0.05),col="gray",linetype="dashed")
  ggsave(paste0("../../out/plot/analysis_R44/33_GSEA_unbiased_",name_anno,"_ranklogfc_","IMMUNE",".pdf"),width = 15,height = 9)
})
