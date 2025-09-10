# AIM ---------------------------------------------------------------------
# run EnrichR on the table of DE from sc analysis

# librarues ---------------------------------------------------------------
library(tidyverse)
library(enrichR)
library(scales)

# DB selection ------------------------------------------------------------
dbs <- listEnrichrDbs()
#filter fo the db of interest
dbs %>%
  filter(str_detect(libraryName,pattern = "Atlas"))

dbs %>%
  filter(str_detect(libraryName,pattern = "Azimuth"))

dbs_db <- c("KEGG_2021_Human","MSigDB_Hallmark_2020","Reactome_Pathways_2024")

# split by direction ------------------------------------------------------

res_file <- "../../out/table/analysis_R44/30_response_MS_vs_CTRL_seurat_sc.tsv"

list_genes_UP <- read_tsv(res_file) %>%
  filter(p_val_adj<0.05 & avg_log2FC > 0.5 & !is.na(gene)) %>%
  split(f = .$annotation) %>%
  map(function(x){
    x %>%
      pull(gene)
  })

list_genes_DOWN <- read_tsv(res_file) %>%
  filter(p_val_adj<0.05 & avg_log2FC < (-0.5)&!is.na(gene)) %>%
  split(f = .$annotation) %>%
  map(function(x){
    x %>%
      pull(gene)
  })

# define the background
# background <- df_modules$feature

# x <- list_res_tot_UP_filter$`DeMye_vs_Norm|clust_5`
list_enrichr_UP <- lapply(list_genes_UP,function(x){
  genes <- x
  # out_enrich <- enrichr(genes, dbs_db,background = background,include_overlap=T)
  out_enrich <- enrichr(genes, dbs_db)
  
  # filter out the annotations without an output
  filter_out <- lapply(out_enrich,function(x){dim(x)[1]}) %>% 
    unlist()
  
  out_enrich[filter_out>0] %>%
    bind_rows(.id = "annotation")
}) %>%
  bind_rows(.id = "comparison")

list_enrichr_UP %>%
  write_tsv("../../out/table/analysis_R44/30_enrichR_DE_UP_sc.tsv")

list_enrichr_DOWN <- lapply(list_genes_DOWN,function(x){
  genes <- x
  # out_enrich <- enrichr(genes, dbs_db,background = background,include_overlap=T)
  out_enrich <- enrichr(genes, dbs_db)
  
  # filter out the annotations without an output
  filter_out <- lapply(out_enrich,function(x){dim(x)[1]}) %>% 
    unlist()
  
  out_enrich[filter_out>0] %>%
    bind_rows(.id = "annotation")
}) %>%
  bind_rows(.id = "comparison")

list_enrichr_DOWN %>%
  write_tsv("../../out/table/analysis_R44/30_enrichR_DE_DOWN_sc.tsv")

# build the plots
plot_list_UP <- list_enrichr_UP %>%
  split(f = .$comparison)

plot_list_DOWN <- list_enrichr_DOWN %>%
  split(f = .$comparison)

# library(scales)
list_plot_UP <- pmap(list(plot_list_UP,names(plot_list_UP)), function(x,y){
  x %>%
    group_by(annotation) %>%
    arrange(P.value) %>%
    dplyr::slice(1:10) %>%
    mutate(Term = str_sub(Term,start = 1,end = 30)) %>%
    mutate(Term = fct_reorder(Term, Combined.Score,.desc = F)) %>%
    # ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
    ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
    scale_color_gradientn(colors = c("red","blue"),
                          values = rescale(c(0,1)),
                          limits = c(0,0.2))+
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))+
    ggtitle(y)
  # scale_color_gradient(low = "red",high = "blue")
  
  #ggsave(paste0("image/enrichR_out_",y,".pdf"),width = 7,height = 15)
})

wrap_plots(list_plot_UP,nrow = 2)
ggsave("../../out/plot/analysis_R44/30_enrichR_DE_UP_sc.pdf",width = 40,height = 20,limitsize = FALSE)

list_plot_DOWN <- pmap(list(plot_list_DOWN,names(plot_list_DOWN)), function(x,y){
  x %>%
    group_by(annotation) %>%
    arrange(P.value) %>%
    dplyr::slice(1:10) %>%
    mutate(Term = str_sub(Term,start = 1,end = 30)) %>%
    mutate(Term = fct_reorder(Term, Combined.Score,.desc = F)) %>%
    # ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
    ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
    scale_color_gradientn(colors = c("red","blue"),
                          values = rescale(c(0,1)),
                          limits = c(0,0.2))+
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))+
    ggtitle(y)
  # scale_color_gradient(low = "red",high = "blue")
  
  #ggsave(paste0("image/enrichR_out_",y,".pdf"),width = 7,height = 15)
})

wrap_plots(list_plot_DOWN,nrow = 2)
ggsave("../../out/plot/analysis_R44/30_enrichR_DE_DOWN_sc.pdf",width = 40,height = 20,limitsize = FALSE)


# check some genes --------------------------------------------------------
read_tsv("../../out/table/analysis_R44/30_enrichR_DE_UP_sc.tsv") %>%
  filter(comparison == "IMMUNE_MS_vs_CTRL") %>%
  filter(annotation == "MSigDB_Hallmark_2020")

# pull the genes of interest
genes <- c("IFITM3", "BST2", "MT2A", "PLSCR1", "BTG1", "MX2", "MX1", "FPR1", "CASP1")

# plot the stats
read_tsv(res_file) %>%
  filter(cluster == "IMMUNE") %>%
  filter(gene %in% genes)

read_tsv("../../out/table/analysis_R44/30_enrichR_DE_UP_sc.tsv") %>%
  filter(str_detect(Genes,pattern = "IGKC")) %>%
  filter(comparison == "IMMUNE_MS_vs_CTRL")
  filter(annotation == "MSigDB_Hallmark_2020")
