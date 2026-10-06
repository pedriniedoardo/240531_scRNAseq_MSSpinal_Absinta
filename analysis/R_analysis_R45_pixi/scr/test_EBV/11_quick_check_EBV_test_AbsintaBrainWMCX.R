# AIM ---------------------------------------------------------------------
# this is the same approach used in 07_quick_check_EBV_test_Spinal.R but for the AbsintaBrainWMCX dataset

# libraries -----------------------------------------------------------------
library(Seurat)
library(tidyverse)

# specify the version of Seurat Assay -------------------------------------
options(Seurat.object.assay.version = "v5")

# parameters ------------------------------------------------------------------
dir_postQC <- "/beegfs/scratch/ric.absinta/ric.absinta/pipeline/scrnaseq_cosr_standard_workflow_AbsintaBrainWMCX_testEBV/results/Seurat_base/object"
postQC_files <- list.files(dir_postQC, pattern = "_default_obj_postQC\\.rds$", full.names = TRUE)
message("per-sample postQC objects found: ", length(postQC_files))

out_table_dir <- "../../out/table/analysis_R45_pixi"
# dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
out_prefix <- "11_quick_check_EBV_test_AbsintaBrainWMCX"

# generic cell-type/lineage marker panel -- same panel reused across analysis_R45_pixi
# (e.g. scr/01_EDA_markerPanels_fullObject.R) for this tissue/disease context
list_markers <- list(
  IMMUNE = c("CX3CR1", "P2RY12", "C3", "CSF1R", "CD74", "C1QB"),
  LYM = c("IGHG1", "CD38", "SKAP1", "CD8A", "CD2"),
  OL = c("MOG", "MBP", "MAG", "NLGN4X", "OLIG1", "OLIG2"),
  ASTRO = c("AQP4", "GFAP", "CD44", "AQP1"),
  NEURONS = c("CUX2", "SYP", "NEFL", "SYT1"),
  VAS = c("VWF", "FLT1", "CLDN5", "PDGFRB"),
  SCHWANN = c("PMP22", "MPZ", "PRX"),
  EPENDYMA = c("CFAP299", "DNAH7", "DNAH9"),
  STROMAL = c("LAMA2", "RBMS3", "CEMIP", "GPC6")
)

# EBV/viral marker panel -- the four functionally-curated groups are identical to test/test_cellranger/test_viral_genome/analysis/R45_pixi/scr/03_plot_markers.R (both references' EBV annotation trace back to the same source GTF/FASTA -- see .../reference/GRCh38-2024-A_plus_EBV/ref/viral.source.txt).
# EBV_OTHER adds every remaining annotated EBV gene in this reference (viral.gtf / this pipeline's own cellranger features.tsv.gz, 86 genes total) that isn't individually curated above, so the check below covers the complete EBV annotation, not just the curated panel -- confirmed to give the same

list_markers_viral <- list(
  EBV_LATENT = c("EBNA-1", "EBNA-2", "EBNA-3A", "EBNA-3B-EBNA-3C", "EBNA-LP", "LMP-1", "RPMS1"),
  EBV_LYTIC_IE = c("BZLF1", "BRLF1"),
  EBV_LYTIC_E = c("BMRF1", "BALF2", "BSLF2-BMLF1", "BALF5", "BHRF1"),
  EBV_LYTIC_L = c("BLLF1", "BcLF1", "BLRF2", "BALF4"),
  EBV_OTHER = c(
    "BNLF2b", "BNLF2a", "BARF1", "BALF1", "BARF0", "BALF3", "A73", "LF2", "LF1", "BILF1", "LF3",
    "BILF2", "BdRF1", "BVRF2", "BVLF1", "BVRF1", "BXRF1", "BXLF2", "BXLF1", "BTRF1", "BcRF1",
    "BDLF3", "BDLF2", "BDLF1", "BGLF2", "BGLF1", "BDLF4", "BDLF3.5", "BGRF1-BDRF1", "BBLF1",
    "BGLF5", "BGLF4", "BGLF3.5", "BGLF3", "BBRF3", "BBLF2-BBLF3", "BBRF2", "BBRF1", "BBLF4",
    "BKRF4", "BKRF3", "BKRF2", "BRRF2", "BRRF1", "BZLF2", "BLLF2", "BLRF1", "BLLF3", "BSRF1",
    "BSLF1", "BMRF2", "BaRF1", "BORF2", "BORF1", "BPLF1", "BOLF1", "BFRF3", "BFRF2", "BFRF1",
    "BFRF1A", "BFLF2", "BFLF1", "BHLF1", "BWRF1", "BCRF1",
    "rna-NC-007605.1:6956..7128", "rna-NC-007605.1:6629..6795", "BNRF1"
  )
)
viral_genes <- unique(unlist(list_markers_viral))
stopifnot(length(viral_genes) == 86)

# read + subset each sample to the marker panels only, then merge ---------------------------------
# RenameCells(add.cell.id=...) keeps every barcode traceable back to its sample after merging (10x
# barcodes repeat across GEM wells/samples, so this also avoids silent barcode collisions)
list_sobj <- map(postQC_files, function(f) {
  sample_id <- str_remove(basename(f), "_default_obj_postQC\\.rds$")
  obj <- readRDS(f)
  genes_keep <- intersect(c(unlist(list_markers), viral_genes), rownames(obj))
  obj <- obj[genes_keep, ]
  RenameCells(obj, add.cell.id = sample_id)
})

merged <- merge(list_sobj[[1]], list_sobj[-1]) %>%
  JoinLayers()
rm(list_sobj); gc()

dim(merged)

# aggregate check: total UMI per viral gene across the whole cohort -------------------------------
viral_genes_in_merged <- intersect(viral_genes, rownames(merged))
rowSums(GetAssayData(merged, assay = "RNA", layer = "data")[viral_genes_in_merged, , drop = FALSE])
rowSums(GetAssayData(merged, assay = "RNA", layer = "counts")[viral_genes_in_merged, , drop = FALSE])

rowSums(GetAssayData(merged, assay = "RNA", layer = "counts")[viral_genes_in_merged, , drop = FALSE])[rowSums(GetAssayData(merged, assay = "RNA", layer = "counts")[viral_genes_in_merged, , drop = FALSE])>0]

# which cell(s), if any, actually carry a nonzero viral count? -------------------------------------
# (generalized over the whole panel, not hardcoded to one gene)
viral_counts <- as.matrix(GetAssayData(merged, assay = "RNA", layer = "counts")[viral_genes_in_merged, , drop = FALSE])
hits <- which(viral_counts != 0, arr.ind = TRUE)

df_hits <- tibble(
  gene = rownames(viral_counts)[hits[, "row"]],
  barcode = colnames(viral_counts)[hits[, "col"]],
  count = viral_counts[hits]
) %>%
  mutate(sample_id = merged@meta.data[barcode, "orig.ident"],
         nCount_RNA = merged@meta.data[barcode, "nCount_RNA"],
         scDblFinder.class = merged@meta.data[barcode, "scDblFinder.class"]) %>%
  arrange(desc(count))

if (nrow(df_hits) == 0) {
  message("no cell in the entire cohort has a nonzero viral-panel count")
} else {
  message(nrow(df_hits), " viral-positive cell/gene hit(s) found across the cohort -- see df_hits")
}
df_hits

# cross-check against the existing manually-annotated full object (different pipeline run/reference, but same underlying FASTQs/barcodes) -- a hit's barcode only matches there if that same cell also passed THAT pipeline's own, independently-run QC
sobj_annotated <- readRDS("/beegfs/scratch/ric.cosr/pedrini.edoardo/project_edoardo/220501_scRNAseq_MSbrain_Absinta/out/object/revision/120_WMCX_ManualClean4_harmonySkipIntegration_AllSoupX_4000_AnnotationSCType_manualAnnotation.rds")
DefaultAssay(sobj_annotated) <- "RNA"
DimPlot(sobj_annotated,group.by = "expertAnno.l2")

# test
rownames(sobj_annotated@meta.data) %>% str_subset(pattern = "GCCATTCTCAACCGAT-1|CTCACTGTCGGTAGAG-1|TACCTGCCAAGCTGTT-1|AGACACTGTGCGACAA-1")
sobj_annotated@meta.data %>%
  rownames_to_column("barcode") %>%
  filter(str_detect(barcode, pattern = "GCCATTCTCAACCGAT-1|CTCACTGTCGGTAGAG-1|TACCTGCCAAGCTGTT-1|AGACACTGTGCGACAA-1"))
