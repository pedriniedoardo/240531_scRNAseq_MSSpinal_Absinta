# AIM ---------------------------------------------------------------------
# extends scr/07_quick_check_EBV_test.R's EBV/viral-signal check from the postQC-filtered called cells to EVERY raw barcode cellranger ever saw (empty droplets + ambient background + real cells), for the "new reference transcriptome" run (ric.absinta pipeline "scrnaseq_cosr_standard_workflow_AbsintaSpinal_snRNAseq_testEBV").
# Earlier in this investigation, summing the whole EBV gene block across the RAW matrix of one sample gave 211 UMI (vs. 1 UMI across ALL postQC-filtered real cells in the whole 25-sample cohort) -- this script finds out which raw barcode(s) that 211 actually sits in, and whether any of them look like real cells (called by cellranger and/or surviving this pipeline's own postQC) or just background.
#
# per-sample raw_feature_bc_matrix directories have ~1.2M barcodes each (vs ~10-15k called cells in the postQC objects), so this is done in two passes per sample to stay tractable:
# 1) collect nCount/nFeature/percent.mt from the FULL (un-subset) raw matrix for every barcode -- this is what actually characterizes a barcode as a real cell vs. background/ambient, and has to be computed before any gene subsetting;
# 2) only THEN trim the matrix down to the 86-gene EBV panel (same panel as scr/07_quick_check_EBV_test.R) and pull out whatever barcodes have a nonzero count. Only barcodes/samples with an actual hit get cross-checked against the officially called-cell list (filtered_feature_bc_matrix) and against the postQC Seurat object, so this stays a "quick check": nothing at raw-barcode scale (~30M barcodes total across the cohort) is ever merged or held in memory beyond the per-sample loop.

# libraries -----------------------------------------------------------------
library(Seurat)
library(tidyverse)

# specify the version of Seurat Assay -------------------------------------
options(Seurat.object.assay.version = "v5")

# parameters ------------------------------------------------------------------
dir_cellranger <- "/beegfs/scratch/ric.absinta/ric.absinta/pipeline/scrnaseq_cosr_standard_workflow_AbsintaSpinal_snRNAseq_testEBV/results/cellranger/merged"
dir_postQC <- "/beegfs/scratch/ric.absinta/ric.absinta/pipeline/scrnaseq_cosr_standard_workflow_AbsintaSpinal_snRNAseq_testEBV/results/Seurat_base/object"
sample_ids <- list.dirs(dir_cellranger, recursive = FALSE, full.names = FALSE)
message("samples found: ", length(sample_ids))

out_table_dir <- "../../out/table/analysis_R45_pixi"
# dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
out_prefix <- "07_quick_check_EBV_test_rawBarcodes"

# EBV/viral marker panel -- identical (all 86 annotated genes) to scr/07_quick_check_EBV_test.R
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

# processing --------------------------------------------------------------
# per sample: load the raw (full-barcode) matrix, characterize every barcode, then trim to viral genes and keep only nonzero hits
# sample_id <- "s2000_090"

# each sample returns a list: a one-row summary (always) + its viral hits (possibly empty)
list_res <- map(sample_ids, function(sample_id) {
  raw_dir <- file.path(dir_cellranger, sample_id, "outs", "raw_feature_bc_matrix")
  if (!dir.exists(raw_dir)) {
    message("sample: ", sample_id, "  -- no raw_feature_bc_matrix, skipping")
    return(list(summary = tibble(sample_id = sample_id, n_raw_barcodes = NA_integer_, n_genes = NA_integer_,
                                 n_viral_genes_present = NA_integer_, n_viral_hits = 0L, total_viral_UMI = 0, n_called_cell_hits = 0L),
                hits = tibble()))
  }
  message("sample: ", sample_id, "  loading raw matrix...")
  
  # genes x ALL barcodes, full transcriptome
  raw_mat <- Read10X(data.dir = raw_dir)  
  # CreateSeuratObject sanitizes "_" -> "-" in feature names internally (Assay5 constructor) 
  # Read10X's raw output still has underscores (BSLF2_BMLF1, EBNA-3B_EBNA-3C, BGRF1_BDRF1, BBLF2_BBLF3, both rna-NC_007605.1:... entries), which wouldn't match list_markers_viral's hyphenated names. Building the object also gets nCount_RNA nFeature_RNA computed for free (min.cells=0/min.features=0 defaults keep every raw barcode).
  # this is also in keeping with the processin in 07_quick_check_EBV_test.R
  obj <- CreateSeuratObject(counts = raw_mat, project = sample_id)
  rm(raw_mat); gc()
  message("  ", ncol(obj), " raw barcodes x ", nrow(obj), " genes")

  # percent.mt still needs its own call
  obj$pct_mt_full <- PercentageFeatureSet(obj, pattern = "^MT-")

  # officially called cells for this sample (cellranger's own cell calling)
  called_barcodes_file <- file.path(dir_cellranger, sample_id, "outs", "filtered_feature_bc_matrix", "barcodes.tsv.gz")
  called_barcodes <- if (file.exists(called_barcodes_file)) {
    readLines(gzfile(called_barcodes_file))
  } else {
    character(0)
  }

  # now trim to the viral panel and pull out nonzero entries (sparse-native, no densifying at ~1.2M-barcode scale)
  n_genes_total <- nrow(obj)
  genes_present <- intersect(viral_genes, rownames(obj))
  # check if there is any gene missing
  setdiff(viral_genes, rownames(obj))
  
  viral_counts <- GetAssayData(obj, assay = "RNA", layer = "counts")[genes_present, , drop = FALSE]
  meta <- obj@meta.data
  rm(obj); gc()

  n_barcodes <- nrow(meta)
  nz <- summary(viral_counts)  # i/j/x for every nonzero entry, sparse-native
  if (nrow(nz) == 0) {
    message("  no nonzero viral-panel entries in this sample's raw matrix")
    return(list(summary = tibble(sample_id = sample_id, n_raw_barcodes = n_barcodes, n_genes = n_genes_total,
                                 n_viral_genes_present = length(genes_present), n_viral_hits = 0L, total_viral_UMI = 0, n_called_cell_hits = 0L),
                hits = tibble()))
  }

  out <- tibble(
    sample_id = sample_id,
    gene = rownames(viral_counts)[nz$i],
    barcode = colnames(viral_counts)[nz$j],
    count = nz$x,
    nCount_RNA_full = meta$nCount_RNA[nz$j],
    nFeature_RNA_full = meta$nFeature_RNA[nz$j],
    pct_mt_full = meta$pct_mt_full[nz$j],
    is_called_cell = colnames(viral_counts)[nz$j] %in% called_barcodes
  )
  message("  ", nrow(out), " nonzero viral-panel (gene, barcode) hit(s), total UMI: ", sum(out$count))
  return(list(summary = tibble(sample_id = sample_id, n_raw_barcodes = n_barcodes, n_genes = n_genes_total,
                               n_viral_genes_present = length(genes_present), n_viral_hits = nrow(out),
                               total_viral_UMI = sum(out$count), n_called_cell_hits = sum(out$is_called_cell)),
              hits = out))
})

# per-sample summary (same numbers as printed in the log), written even when there are no hits
df_summary <- map_dfr(list_res, "summary")
df_hits <- map_dfr(list_res, "hits")
df_summary
write_tsv(df_summary, file.path(out_table_dir, paste0(out_prefix, "_perSampleSummary.tsv")))

# cross-check hits against this pipeline's own postQC Seurat object (doublet-called, QC-filtered) only opened for samples that actually had a raw hit, not all 25, to stay "quick"
if (nrow(df_hits) == 0) {
  message("no barcode in the entire raw dataset (all droplets, all samples) has a nonzero viral-panel count")
  df_hits_full <- df_hits
} else {
  hit_samples <- unique(df_hits$sample_id)
  # sid <- "Sample_2015_017_300"
  df_postQC_barcodes <- map_dfr(hit_samples, function(sid) {
    f <- file.path(dir_postQC, paste0(sid, "_default_obj_postQC.rds"))
    if (!file.exists(f)) return(tibble(sample_id = character(0), barcode = character(0)))
    obj <- readRDS(f)
    out <- obj@meta.data %>% rownames_to_column("barcode") %>% dplyr::select(barcode,
                                                                             sample_id=orig.ident,
                                                                             scDblFinder.class,
                                                                             nCount_RNA,
                                                                             nFeature_RNA,
                                                                             percent.mt,
                                                                             percent.ribo,
                                                                             percent.globin,
                                                                             discard_multi,
                                                                             discard_single)
    rm(obj); gc()
    return(out)
  })
  df_hits_full <- df_hits %>%
    mutate(in_postQC = paste(sample_id, barcode) %in% paste(df_postQC_barcodes$sample_id, df_postQC_barcodes$barcode)) %>%
    left_join(df_postQC_barcodes,by = c("sample_id","barcode")) %>%
    arrange(desc(count))
  message(nrow(df_hits_full), " viral-positive (gene, barcode) hit(s) found across all raw barcodes, ",
          "all samples -- ", sum(df_hits$is_called_cell), " in cellranger's called-cell list, ",
          sum(df_hits_full$in_postQC), " surviving this pipeline's postQC")
}
df_hits_full

write_tsv(df_hits_full, file.path(out_table_dir, paste0(out_prefix, "_viralPositiveBarcodes_raw.tsv")))
