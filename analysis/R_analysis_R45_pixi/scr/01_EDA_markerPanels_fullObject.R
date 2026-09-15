# AIM ---------------------------------------------------------------------
# EDA on the final subcluster-filtered full object (output of R44's 29_apply_filter_annotate_fullObject.R):
# 1) DC vs microglia vs BAM markers in IMMUNE clusters
# 2) ciliated genes in ASTRO clusters
# 3) lymphocyte subtype markers in LYM clusters
# 4) IFN-signature validity across OPC + OLIGO + IMMUNE clusters (OPC6/OLIGO9/IMMUNE10)
# 5) OPALIN/oligodendrocyte maturation panel (Jakel 2019, Castelo-Branco lab) in the full object,
#    OLIGO shown at subcluster resolution against all other cell types pooled as reference

# added also the panel propose by Sofia fro the LYM and the NEU subclusters

# LIBRARIES -----------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(patchwork)
library(cowplot)

# PARAMETERS ------------------------------------------------------------------
options(stringsAsFactors = FALSE)
options(Seurat.object.assay.version = "v5")

# read in the data --------------------------------------------------------
sobj <- readRDS("../../out/object/analysis_R44/29_sobj_integrated_cleanup_manualAnnotation_subclusterFiltered.rds")
DefaultAssay(sobj) <- "RNA"

# check the dataset
DimPlot(sobj)

# save the full metadata with the current draft anntoation from Martina
sobj@meta.data %>% group_by(cell_id_subcluster,cell_id_subcluster2) %>% summarise(n = n()) %>%
  write_tsv("../../out/table/analysis_R45_pixi/01_meta_29_sobj_integrated_cleanup_manualAnnotation_subclusterFiltered.tsv")

out_plot_dir <- "../../out/plot/analysis_R45_pixi"

# helper A: DotPlot a named gene-panel list for the target populations only, at fine resolution (group_col) -- the original focal-only view --------
plot_marker_dotplot_focal <- function(obj,
                                      gene_list,
                                      populations,
                                      group_col = "cell_id_subcluster2",
                                      pop_col = "cell_id",
                                      title) {
  genes_present <- lapply(gene_list, function(g) intersect(g, rownames(obj)))
  missing <- setdiff(unlist(gene_list), unlist(genes_present))
  if (length(missing) > 0) message("Missing from object, dropped: ", paste(missing, collapse = ", "))

  sub <- subset(obj, cells = colnames(obj)[obj[[pop_col, drop = TRUE]] %in% populations])

  dp <- DotPlot(sub, features = unique(unlist(genes_present)), group.by = group_col,
                dot.scale = 8, cluster.idents = TRUE)

  df <- lapply(genes_present, function(g) dp$data %>% filter(features.plot %in% g)) %>%
    bind_rows(.id = "gene_group") %>%
    mutate(features.plot = factor(features.plot, levels = unique(unlist(genes_present)))) %>%
    mutate(gene_group = factor(gene_group, levels = names(gene_list)))

  p <- ggplot(df, aes(x = features.plot, y = id)) +
    geom_point(aes(size = pct.exp, colour = avg.exp.scaled)) +
    scale_radius(range = c(0, 8)) +
    facet_grid(~ gene_group, scales = "free", space = "free") +
    scale_color_gradient(low = "lightgrey", high = "blue") +
    theme_cowplot() +
    labs(title = title, x = NULL, y = group_col) +
    theme(strip.background = element_blank(),
          strip.text = element_text(angle = 90),
          axis.text.x = element_text(angle = 90, hjust = 1))

  return(p)
}

# helper B: DotPlot a named gene-panel list for the target populations at fine resolution (group_col), plus every non-target cell_id pooled as "Ref_<cell_id>" below a dashed separator -- so the panel is judged against general background too, not just within the target subclusters (same focal-vs-reference idiom as task 5) --------
plot_marker_dotplot_vs_ref <- function(obj,
                                       gene_list,
                                       populations,
                                       group_col = "cell_id_subcluster2",
                                       pop_col = "cell_id",
                                       title) {
  genes_present <- lapply(gene_list, function(g) intersect(g, rownames(obj)))
  missing <- setdiff(unlist(gene_list), unlist(genes_present))
  if (length(missing) > 0) message("Missing from object, dropped: ", paste(missing, collapse = ", "))

  is_focal <- obj[[pop_col, drop = TRUE]] %in% populations
  obj$plot_id <- case_when(
    is_focal ~ as.character(obj[[group_col, drop = TRUE]]),
    TRUE ~ paste0("Ref_", obj[[pop_col, drop = TRUE]]))

  focal_ids <- sort(unique(obj$plot_id[is_focal]))
  ref_ids   <- sort(unique(obj$plot_id[!is_focal]))
  id_levels <- c(focal_ids, ref_ids)
  obj$plot_id <- factor(obj$plot_id, levels = id_levels)
  Idents(obj) <- "plot_id"

  sep_y <- length(focal_ids) + 0.5

  dp <- DotPlot(obj, features = unique(unlist(genes_present)), dot.scale = 8, cluster.idents = FALSE)

  df <- lapply(genes_present, function(g) dp$data %>% filter(features.plot %in% g)) %>%
    bind_rows(.id = "gene_group") %>%
    mutate(features.plot = factor(features.plot, levels = unique(unlist(genes_present))),
           id = factor(id, levels = id_levels)) %>%
    mutate(gene_group = factor(gene_group, levels = names(gene_list)))

  p <- ggplot(df, aes(x = features.plot, y = id)) +
    geom_point(aes(size = pct.exp, colour = avg.exp.scaled)) +
    scale_radius(range = c(0, 8)) +
    facet_grid(~ gene_group, scales = "free", space = "free") +
    geom_hline(yintercept = sep_y, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    scale_color_gradient(low = "lightgrey", high = "blue") +
    theme_cowplot() +
    labs(title = title, x = NULL, y = "plot_id",
         caption = paste0(paste(populations, collapse = "/"),
                           " subclusters below the line | Ref_* = pooled other cell_id above")) +
    theme(strip.background = element_blank(),
          strip.text = element_text(angle = 90),
          axis.text.x = element_text(angle = 90, hjust = 1),
          plot.caption = element_text(hjust = 0, colour = "grey40", size = 7))

  return(p)
}

# 1. DC vs microglia vs BAM -- IMMUNE clusters -------------------------------
dc_panel <- list(
  cDC1                       = c("CLEC9A", "XCR1", "BATF3"),
  `cDC2/conventional_DC`     = c("CD1C", "FCER1A", "CLEC10A", "CD1E"),
  `DC_activation/migratory`  = c("CD83", "COTL1", "CCR7"),
  homeostatic_microglia      = c("P2RY12", "TMEM119", "SALL1", "GPR34", "CX3CR1"),
  `BAM/perivascular_mac`     = c("MRC1", "LYVE1", "SELENOP", "CD163", "STAB1")
)

p1 <- plot_marker_dotplot_focal(
  sobj, dc_panel, populations = "IMMUNE",
  title = "DC vs microglia vs BAM markers -- IMMUNE subclusters (focal only)"
)
ggsave(plot = p1, filename = file.path(out_plot_dir, "01_dotplot_DC_markers_IMMUNE_focalOnly.pdf"),
       width = 12, height = 6)

p2 <- plot_marker_dotplot_vs_ref(
  sobj, dc_panel, populations = "IMMUNE",
  title = "DC vs microglia vs BAM markers -- IMMUNE subclusters vs reference"
)
ggsave(plot = p2, filename = file.path(out_plot_dir, "01_dotplot_DC_markers_IMMUNE_vsRef.pdf"),
       width = 12, height = 8)

# 2. Ciliated astrocyte genes -- ASTRO clusters ------------------------------
ciliated_panel <- list(
  ciliated_astrocyte = c("FRMPD2", "ADGB", "SPAG17", "FOXJ1", "DNAH11", "DNAH6", "CFAP54", "CFAP299")
)

p3 <- plot_marker_dotplot_focal(
  sobj, ciliated_panel, populations = "ASTRO",
  title = "Ciliated astrocyte genes -- ASTRO subclusters (focal only)"
)
ggsave(plot = p3, filename = file.path(out_plot_dir, "01_dotplot_ciliatedGenes_ASTRO_focalOnly.pdf"),
       width = 8, height = 5)

p4 <- plot_marker_dotplot_vs_ref(
  sobj, ciliated_panel, populations = "ASTRO",
  title = "Ciliated astrocyte genes -- ASTRO subclusters vs reference"
)
ggsave(plot = p4, filename = file.path(out_plot_dir, "01_dotplot_ciliatedGenes_ASTRO_vsRef.pdf"),
       width = 8, height = 7)

# 3. Lymphocyte subtype markers -- LYM clusters ------------------------------
# note: cell_id_subcluster(2) uses the pooled "LYM|..." prefix for these cells even though cell_id itself is split into "B CELLS"/"T CELLS"
lym_panel <- list(
  `MAIT/mucosal`  = c("TRAV1-2", "KLRB1", "SLC4A10", "DPP4"),
  pan_T           = c("TRBC1", "TRAC", "IL7R", "LTB"),
  `CD4/CD8`       = c("CD4", "CD8A", "CD8B"),
  trafficking     = c("CCR6", "CXCR6"),
  `NK/cytotoxic`  = c("KLRD1", "NCR1", "FCGR3A", "TYROBP")
)

# confirm the population
# group_col = "cell_id_subcluster2",
# pop_col = "cell_id",
sobj@meta.data %>%
  filter(str_detect(cell_id_subcluster2,pattern = "LYM")) %>%
  group_by(cell_id) %>%
  summarise(n = n())

sobj@meta.data %>%
  filter(cell_id %in% c("B CELLS","T CELLS")) %>%
  dim()

sobj@meta.data %>%
  filter(str_detect(cell_id_subcluster2,pattern = "LYM")) %>%
  dim()

p5 <- plot_marker_dotplot_focal(
  sobj, lym_panel, populations = c("B CELLS","T CELLS"),
  title = "Lymphocyte subtype markers -- LYM subclusters (focal only)"
)
ggsave(plot = p5, filename = file.path(out_plot_dir, "01_dotplot_lymphocyte_markers_LYM_focalOnly.pdf"),
       width = 10, height = 6)

p6 <- plot_marker_dotplot_vs_ref(
  sobj, lym_panel, populations = c("B CELLS","T CELLS"),
  title = "Lymphocyte subtype markers -- LYM subclusters vs reference"
)
ggsave(plot = p6, filename = file.path(out_plot_dir, "01_dotplot_lymphocyte_markers_LYM_vsRef.pdf"),
       width = 10, height = 6)

# add the panel suggested by Sofia to check the CD4 resident cells
memory_tcell_panel <- list(
  `Central memory (TCM / TSM)` =
    c("CCR7", "CD27", "PTPRC", "SELL"),
  # CCR7+, CD27+, CD45RO+/CD45RA- (PTPRC isoforms), CD62L+ (SELL)
  
  `Effector memory (TEM(RA) / TDB/WM)` =
    c("CCR7", "PTPRC", "SELL"),
  # CCR7-, CD45RA+/- (PTPRC isoform), CD62L- (SELL)
  
  `Resident memory (TRM / TDRM)` =
    c("CCR7", "ITGA1", "SELL", "CD69", "ITGAE", "CXCR6")
  # CCR7-, CD49a+ (ITGA1), CD62L- (SELL), CD69+, CD103+/- (ITGAE), CXCR6+
)

p7 <- plot_marker_dotplot_focal(
  sobj, memory_tcell_panel, populations = c("B CELLS","T CELLS"),
  title = "Lymphocyte subtype markers -- LYM subclusters (focal only)"
)
ggsave(plot = p7, filename = file.path(out_plot_dir, "01_dotplot_lymphocyte_markers_sofi_LYM_focalOnly.pdf"),
       width = 10, height = 6)

p8 <- plot_marker_dotplot_vs_ref(
  sobj, memory_tcell_panel, populations = c("B CELLS","T CELLS"),
  title = "Lymphocyte subtype markers -- LYM subclusters vs reference"
)
ggsave(plot = p8, filename = file.path(out_plot_dir, "01_dotplot_lymphocyte_markers_sofi_LYM_vsRef.pdf"),
       width = 10, height = 8)

# show also the subcluster of the LYM cells
sobj_lym <- readRDS("../../out/object/analysis_R44/27_LYM_subcluster_HarmonySample.rds")
sobj_lym <- AddMetaData(sobj_lym,metadata = sobj@meta.data[,c("cell_id_subcluster",
                                                              "cell_id_subcluster2")])
DimPlot(sobj_lym,group.by = "cell_id_subcluster2",label = T, repel = T) + NoLegend()

# 4. IFN signature -- OPC + OLIGO + IMMUNE together, one DotPlot -------------
# cell_id_subcluster2 already carries the population prefix (e.g. "OPC|OPC-IFN|6"). Rather than just ordering clusters back-to-back, pull DotPlot's own long-format data (id/features.plot/pct.exp/avg.exp.scaled) and rebuild the plot with facet_grid(cell_type ~ ., scales = "free", space = "free") so OPC/OLIGO/IMMUNE render as visually separated row-blocks -- the candidate IFN clusters (OPC-IFN/Oligo-IFN/MIMS-IFN) compare directly against their own population's other clusters, with a clear visual break between populations.
ifn_panel <- list(
  IFN_signature = c(
    "ISG15", "MX1", "MX2", "IFIT1", "IFIT2", "IFIT3", "IFITM1", "IFITM3", "RSAD2",
    "OAS1", "OAS2", "OAS3", "OASL", "STAT1", "STAT2", "IRF1", "IRF7", "JAK1", "JAK2",
    "GBP1", "GBP5", "EPSTI1", "PARP14", "XAF1", "APOL1", "APOL2", "APOL3", "NLRC5", "B2M"
  )
)

ifn_genes_present <- intersect(ifn_panel$IFN_signature, rownames(sobj))
missing_ifn <- setdiff(ifn_panel$IFN_signature, ifn_genes_present)
if (length(missing_ifn) > 0) message("Missing from object, dropped: ", paste(missing_ifn, collapse = ", "))

sobj_ifn <- subset(sobj, cells = colnames(sobj)[sobj$cell_id %in% c("OPC", "OLIGO", "IMMUNE")])
Idents(sobj_ifn) <- "cell_id_subcluster2"

dp_ifn <- DotPlot(sobj_ifn, features = ifn_genes_present)

id_to_pop_ifn <- sobj_ifn@meta.data %>%
  distinct(cell_id_subcluster2, cell_id) %>%
  mutate(cluster_num = suppressWarnings(as.integer(str_extract(cell_id_subcluster2, "[0-9]+$")))) %>%
  rename(id = cell_id_subcluster2, cell_type = cell_id) %>%
  mutate(cell_type = factor(cell_type, levels = c("OPC", "OLIGO", "IMMUNE"))) %>%
  arrange(cell_type, cluster_num)

df_plot_ifn <- dp_ifn$data %>%
  left_join(id_to_pop_ifn, by = "id") %>%
  mutate(id = factor(id, levels = id_to_pop_ifn$id))

p_ifn <- ggplot(df_plot_ifn, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, col = avg.exp.scaled)) +
  scale_size(range = c(0, 6)) +
  facet_grid(cell_type ~ ., scales = "free", space = "free") +
  theme_cowplot() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(hjust = 1, angle = 90),
        strip.text.y = element_text(angle = 0)) +
  scale_color_gradient(low = "lightgrey", high = "blue") +
  labs(title = "IFN signature across OPC / OLIGO / IMMUNE subclusters",
       x = "IFN signature genes", y = "cell_id_subcluster2")

ggsave(plot = p_ifn, filename = file.path(out_plot_dir, "01_dotplot_IFNsignature_OPC_OLIGO_IMMUNE.pdf"),
       width = 16, height = 8)

# 5. OPALIN / oligodendrocyte maturation (Jakel 2019, Castelo-Branco lab) ----
# focal-vs-reference calibrated DotPlot: OLIGO cells shown at subcluster resolution, everything else pooled as "Ref_<cell_id>" -- run twice, once per fine-label column. Literature check: OPALIN marks the mature "Oligo6" state in Jakel et al. 2019 Nature (human MS white matter oligodendrocyte heterogeneity, Castelo-Branco lab); KLK6 marks "Oligo5"; RASGRF1 marks "Oligo1". Full Oligo1-6 subtype assignment would need the paper's supplementary marker table for high confidence -- this panel is a maturation/heterogeneity screen, not a claim of exact subtype identity.
opalin_panel <- list(
  oligo_maturation = c(
    "PDGFRA", "CSPG4", "GPR17", "BCAS1", "OLIG2", "SOX10", "MYRF", "PLP1", "MBP",
    "MOBP", "MAG", "MOG", "OPALIN", "TPPP", "CNP", "GFAP", "DPP10",
    "KLK6", "RASGRF1", "RBFOX1", "ANLN"
  )
)

plot_focal_vs_ref_dotplot <- function(obj,
                                      focal_pop,
                                      focal_label_col,
                                      panel,
                                      title) {
  obj$plot_id <- case_when(
    obj$cell_id == focal_pop ~ as.character(obj[[focal_label_col, drop = TRUE]]),
    TRUE ~ paste0("Ref_", obj$cell_id)
  )

  focal_ids <- sort(unique(obj$plot_id[obj$cell_id == focal_pop]))
  ref_ids   <- sort(unique(obj$plot_id[obj$cell_id != focal_pop]))
  id_levels <- c(focal_ids, ref_ids)
  obj$plot_id <- factor(obj$plot_id, levels = id_levels)
  Idents(obj) <- "plot_id"

  panel_filt <- lapply(panel, function(gs) intersect(gs, rownames(obj)))
  missing <- setdiff(unlist(panel), unlist(panel_filt))
  if (length(missing) > 0) message("Missing from object, dropped: ", paste(missing, collapse = ", "))
  panel_filt <- panel_filt[lengths(panel_filt) > 0]

  sep_y <- length(focal_ids) + 0.5

  dp <- DotPlot(obj, features = unique(unlist(panel_filt)), dot.scale = 6, cluster.idents = FALSE)

  df <- lapply(panel_filt, function(g) dp$data %>% filter(features.plot %in% g)) %>%
    bind_rows(.id = "gene_group") %>%
    mutate(features.plot = factor(features.plot, levels = unique(unlist(panel_filt))),
           id = factor(id, levels = id_levels))

  p <- ggplot(df, aes(x = features.plot, y = id)) +
    geom_point(aes(size = pct.exp, colour = avg.exp.scaled)) +
    scale_radius(range = c(0, 6)) +
    facet_grid(~ gene_group, scales = "free", space = "free") +
    geom_hline(yintercept = sep_y, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    scale_color_gradient(low = "lightgrey", high = "blue") +
    theme_cowplot() +
    labs(title = title, x = NULL, y = "plot_id",
         caption = paste0(focal_pop, " subclusters below the line | Ref_* = pooled other cell_id above")) +
    theme(strip.background = element_blank(),
          strip.text = element_text(angle = 90),
          axis.text.x = element_text(angle = 90, hjust = 1),
          plot.caption = element_text(hjust = 0, colour = "grey40", size = 7))

  return(p)
}

p9 <- plot_focal_vs_ref_dotplot(
  sobj, "OLIGO", "cell_id_subcluster", opalin_panel,
  title = "OLIGO subclusters vs reference cell types -- OPALIN/maturation panel"
)
ggsave(plot = p9, filename = file.path(out_plot_dir, "01_dotplot_OPALIN_oligoMaturation_calibrated_subclusterLabel.pdf"),
       width = 14, height = 10)

p10 <- plot_focal_vs_ref_dotplot(
  sobj, "OLIGO", "cell_id_subcluster2", opalin_panel,
  title = "OLIGO subclusters (with res0.4 id) vs reference cell types -- OPALIN/maturation panel"
)
ggsave(plot = p10, filename = file.path(out_plot_dir, "01_dotplot_OPALIN_oligoMaturation_calibrated_subclusterLabel2.pdf"),
       width = 14, height = 14)


# 6. NEU subset -----------------------------------------------------------
# Explore the panel suggested by Sofia for the neurons subset
yadav_gene_class_panel <- list(
  # Dorsal excitatory (Ex-Dorsal-1..12)
  `Ex-Dorsal`  = c("MAFA", "COL6A2", "MYO10", "CBLN1", "NR2F1", "NTS", "GABRB2",
                   "MAF", "PDE11A", "RELN", "IQGAP2", "NMUR2", "TAC1", "SOX5",
                   "IGF1", "CARTPT", "LMO3", "NMU", "ISM1", "CALCB"),
  
  # Dorsal inhibitory (Inh-Dorsal-1..9)
  `Inh-Dorsal` = c("CAPN8", "NR2F2", "LINC01197", "PENK", "STK32B", "RORB",
                   "P2RY1", "CBLN4", "ADARB2", "CDHR3", "NFATC1", "CDH7",
                   "PDYN", "SULF1", "LEF1", "CCK", "PROX1", "NPY", "SAMD3",
                   "HMCN1"),
  
  # Mid/ventral excitatory (Ex-M-1..4, Ex-V-1..3)
  `Ex-M`       = c("TAC1", "LMX1B", "ZFHX3", "SATB2", "TAC3", "NFIX", "CDH23",
                   "PRDM8", "ISL1", "ONECUT2", "POU6F2", "RNF220", "FOXP2",
                   "NR4A2", "CDH9"),
  
  # Mid/ventral inhibitory (Inh-M-1..2, Inh-V-1..3)
  `Inh-M`      = c("GAD2", "TFAP2B", "POU6F2", "DRD2", "PROX1", "NFIX",
                   "GATA3", "FOXP2", "ZFHX3", "ASAH2"),
  
  # Motoneurons (not part of the dendrogram; markers from Fig. 3D/6A and text)
  Motoneurons  = c("CHAT", "SLC5A7", "SOD1", "TUBA4A", "NEFL", "NEFH", "NEFM",
                   "STMN2", "PRPH", "SPP1")
)

yadav_gene_class_panel_full <- list(
  # Dorsal excitatory (Ex-Dorsal-1..12)
  `Ex-Dorsal`  = c("MAFA", "COL6A2", "MYO10", "CBLN1", "NR2F1", "NTS", "GABRB2",
                   "MAF", "PDE11A", "RELN", "IQGAP2", "NMUR2", "TAC1", "SOX5",
                   "IGF1", "CARTPT", "LMO3", "NMU", "ISM1", "CALCB"),
  
  # Dorsal inhibitory (Inh-Dorsal-1..9)
  `Inh-Dorsal` = c("CAPN8", "NR2F2", "LINC01197", "PENK", "STK32B", "RORB",
                   "P2RY1", "CBLN4", "ADARB2", "CDHR3", "NFATC1", "CDH7",
                   "PDYN", "SULF1", "LEF1", "CCK", "PROX1", "NPY", "SAMD3",
                   "HMCN1"),
  
  # Mid/ventral excitatory (Ex-M-1..4, Ex-V-1..3)
  `Ex-M`       = c("TAC1", "LMX1B", "ZFHX3", "SATB2", "TAC3", "NFIX", "CDH23",
                   "PRDM8", "ISL1", "ONECUT2", "POU6F2", "RNF220", "FOXP2",
                   "NR4A2", "CDH9"),
  
  # Mid/ventral inhibitory (Inh-M-1..2, Inh-V-1..3)
  `Inh-M`      = c("GAD2", "TFAP2B", "POU6F2", "DRD2", "PROX1", "NFIX",
                   "GATA3", "FOXP2", "ZFHX3", "ASAH2"),
  
  # Motoneurons (not part of the dendrogram; markers from Fig. 3D/6A and text)
  Motoneurons  = c("CHAT", "SLC5A7", "SOD1", "TUBA4A", "NEFL", "NEFH", "NEFM",
                   "STMN2", "PRPH", "SPP1"),
  # add also generic markers
  IMMUNE = c("CX3CR1","P2RY12","C3","CSF1R", "CD74","C1QB"),
  LYM = c("IGHG1", "CD38","SKAP1", "CD8A", "CD2"),
  OL = c("MOG","MBP","MAG","NLGN4X","OLIG1","OLIG2"),
  ASTRO = c("AQP4", "GFAP", "CD44", "AQP1"),
  VAS = c("VWF","FLT1","CLDN5","PDGFRB"),
  SCHWANN = c("PMP22","MPZ","PRX"),
  EPENDYMA = c("CFAP299","DNAH7","DNAH9"),
  STROMAL = c("LAMA2","RBMS3","CEMIP","GPC6")
)

zhang_gene_class_panel <- list(
  # Excitatory (C0, C3, C4, C5, C7, C10, C13, C14, C16, C17, C19)
  Excitatory = c("NEFM", "NEFH", "NEFL", "AC018742.1", "KLHL1", "DPP10",
                 "TAFA1", "SOX5", "ADAMTSL1", "PBX3", "LINC02742", "EBF2",
                 "TLN2", "KIRREL3", "AL355612.1", "GABRB2", "DOCK4", "MCOLN3",
                 "CNTNAP3B", "TLE1", "WDR11", "ATP2B4", "BNC2", "SRPX2",
                 "LINC00378", "SCHLAP1", "IQGAP2", "PLPP2", "MORC1", "CA8",
                 "XIST"),
  
  # Inhibitory (C2, C6, C8, C9, C11, C12, C15, C18)
  Inhibitory = c("NXPH1", "TLL1", "TCF4", "CCBE1", "GAD1", "GRIK2", "SEMA3A",
                 "MAN1A2", "RHOBTB2", "ADARB2", "ADAMTS17", "GRM3", "SAMD3",
                 "ADAMTS16", "AL356737.12", "NR2F2", "PAX5", "AC093765.2",
                 "IRAK3", "AC105916.1", "CAMKMT", "CCR3", "MSC-AS1", "RECK"),
  
  # Motor (C20)
  Motor      = c("SLC5A7", "AC010967.1", "PCA3"),
  
  # Mixed functional status (C1)
  Mixed      = c("HTR2C", "PCDH11X", "ZFHX3")
)


zhang_gene_class_panel_full <- list(
  # Excitatory (C0, C3, C4, C5, C7, C10, C13, C14, C16, C17, C19)
  Excitatory = c("NEFM", "NEFH", "NEFL", "AC018742.1", "KLHL1", "DPP10",
                 "TAFA1", "SOX5", "ADAMTSL1", "PBX3", "LINC02742", "EBF2",
                 "TLN2", "KIRREL3", "AL355612.1", "GABRB2", "DOCK4", "MCOLN3",
                 "CNTNAP3B", "TLE1", "WDR11", "ATP2B4", "BNC2", "SRPX2",
                 "LINC00378", "SCHLAP1", "IQGAP2", "PLPP2", "MORC1", "CA8",
                 "XIST"),
  
  # Inhibitory (C2, C6, C8, C9, C11, C12, C15, C18)
  Inhibitory = c("NXPH1", "TLL1", "TCF4", "CCBE1", "GAD1", "GRIK2", "SEMA3A",
                 "MAN1A2", "RHOBTB2", "ADARB2", "ADAMTS17", "GRM3", "SAMD3",
                 "ADAMTS16", "AL356737.12", "NR2F2", "PAX5", "AC093765.2",
                 "IRAK3", "AC105916.1", "CAMKMT", "CCR3", "MSC-AS1", "RECK"),
  
  # Motor (C20)
  Motor      = c("SLC5A7", "AC010967.1", "PCA3"),
  
  # Mixed functional status (C1)
  Mixed      = c("HTR2C", "PCDH11X", "ZFHX3"),
  
  # add also generic markers
  IMMUNE = c("CX3CR1","P2RY12","C3","CSF1R", "CD74","C1QB"),
  LYM = c("IGHG1", "CD38","SKAP1", "CD8A", "CD2"),
  OL = c("MOG","MBP","MAG","NLGN4X","OLIG1","OLIG2"),
  ASTRO = c("AQP4", "GFAP", "CD44", "AQP1"),
  VAS = c("VWF","FLT1","CLDN5","PDGFRB"),
  SCHWANN = c("PMP22","MPZ","PRX"),
  EPENDYMA = c("CFAP299","DNAH7","DNAH9"),
  STROMAL = c("LAMA2","RBMS3","CEMIP","GPC6")
)


# confirm the population
# group_col = "cell_id_subcluster2",
# pop_col = "cell_id",
sobj@meta.data %>%
  filter(str_detect(cell_id_subcluster2,pattern = "NEU")) %>%
  group_by(cell_id) %>%
  summarise(n = n())

# show also the subcluster of the LYM cells
sobj_neu <- readRDS("../../out/object/analysis_R44/27_NEU_subcluster_HarmonySample.rds")
DimPlot(sobj_neu,group.by = "RNA_snn_res.0.1",label = T, repel = T) + NoLegend()

# add the subcluster clustering info from the neuron to the full object
sobj_test <- sobj
df_sobj_neu <- data.frame(NEU_res0.1 = sobj_neu@meta.data[,c("RNA_snn_res.0.1")],
                          row.names = rownames(sobj_neu@meta.data))
sobj_test <- AddMetaData(sobj_test,metadata = df_sobj_neu)

p11 <- plot_marker_dotplot_focal(
  sobj_test, yadav_gene_class_panel,
  pop_col = "cell_id",
  group_col = "NEU_res0.1",
  populations = c("NEU"),
  title = "Neuron markers -- NEU subclusters (focal only)")
ggsave(plot = p11, filename = file.path(out_plot_dir, "01_dotplot_yadav_markers_NEU_focalOnly.pdf"),
       width = 15, height = 8)

p12 <- plot_marker_dotplot_vs_ref(
  sobj_test, yadav_gene_class_panel_full,
  pop_col = "cell_id",
  group_col = "NEU_res0.1",
  populations = c("NEU"),
  title = "Neuron + generic markers -- NEU subclusters vs reference")
ggsave(plot = p12, filename = file.path(out_plot_dir, "01_dotplot_yadav_markers_NEU_vsRef.pdf"),
       width = 25, height = 8)

p13 <- plot_marker_dotplot_focal(
  sobj_test, zhang_gene_class_panel,
  pop_col = "cell_id",
  group_col = "NEU_res0.1",
  populations = c("NEU"),
  title = "Neuron markers -- NEU subclusters (focal only)")
ggsave(plot = p13, filename = file.path(out_plot_dir, "01_dotplot_zhang_markers_NEU_focalOnly.pdf"),
       width = 15, height = 8)

p14 <- plot_marker_dotplot_vs_ref(
  sobj_test, zhang_gene_class_panel_full,
  pop_col = "cell_id",
  group_col = "NEU_res0.1",
  populations = c("NEU"),
  title = "Neuron + generic markers -- NEU subclusters vs reference")
ggsave(plot = p14, filename = file.path(out_plot_dir, "01_dotplot_zhang_markers_NEU_vsRef.pdf"),
       width = 25, height = 8)

# 7. STROMAL --------------------------------------------------------------
# explore a panel for the stomal for fibroblast and red plood cells
sobj_stromal <- readRDS("../../out/object/analysis_R44/27_STROMAL_subcluster_HarmonySample.rds")
DimPlot(sobj_stromal,group.by = "RNA_snn_res.0.4",label = T, repel = T) + NoLegend()

qc_stomal_panel <- list(
  Fibro     = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA"),
  RBC = c("HBB", "HBA1", "HBA2", "ALAS2", "GYPA"),
  # add also generic markers
  IMMUNE = c("CX3CR1","P2RY12","C3","CSF1R", "CD74","C1QB"),
  LYM = c("IGHG1", "CD38","SKAP1", "CD8A", "CD2"),
  OL = c("MOG","MBP","MAG","NLGN4X","OLIG1","OLIG2"),
  ASTRO = c("AQP4", "GFAP", "CD44", "AQP1"),
  NEURONS = c("CUX2","SYP", "NEFL","SYT1"),
  VAS = c("VWF","FLT1","CLDN5","PDGFRB"),
  SCHWANN = c("PMP22","MPZ","PRX"),
  EPENDYMA = c("CFAP299","DNAH7","DNAH9"),
  STROMAL = c("LAMA2","RBMS3","CEMIP","GPC6")
)

# add the subcluster clustering info from the neuron to the full object
sobj_test <- sobj
df_sobj_stromal <- data.frame(STROMAL_res0.4 = sobj_stromal@meta.data[,c("RNA_snn_res.0.4")],
                              row.names = rownames(sobj_stromal@meta.data))
sobj_test <- AddMetaData(sobj_test,metadata = df_sobj_stromal)

p15 <- plot_marker_dotplot_focal(
  sobj_test, qc_stomal_panel,
  pop_col = "cell_id",
  group_col = "STROMAL_res0.4",
  populations = c("STROMAL"),
  title = "Stromal markers -- STROMAL subclusters (focal only)")
ggsave(plot = p15, filename = file.path(out_plot_dir, "01_dotplot_stromal_markers_STROMAL_focalOnly.pdf"),
       width = 15, height = 8)

p16 <- plot_marker_dotplot_vs_ref(
  sobj_test, qc_stomal_panel,
  pop_col = "cell_id",
  group_col = "STROMAL_res0.4",
  populations = c("STROMAL"),
  title = "Stromal + generic markers -- STROMAL subclusters vs reference")
ggsave(plot = p16, filename = file.path(out_plot_dir, "01_dotplot_stromal_markers_STROMAL_vsRef.pdf"),
       width = 25, height = 8)

