scipts <- c("scr/17_subcluster_AST.R",
            # "scr/17_subcluster_IMMUNE.R"
            "scr/17_subcluster_LYM.R",
            "scr/17_subcluster_NEU.R",
            "scr/17_subcluster_OLIGO.R",
            "scr/17_subcluster_OPC.R",
            "scr/17_subcluster_OTHER.R",
            "scr/17_subcluster_VAS.R")

lapply(scipts, function(sc){
  source(sc)
  rm(list = ls())
  gc()
})