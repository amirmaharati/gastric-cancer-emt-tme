############################################################
# EMT-Associated lncRNA Prioritization
#
# Project:
#   Gastric Cancer EMT–TME Analysis
#
# Dataset:
#   GSE183904
#
# Inputs:
#   objects/seurat_epithelial_EMT_GSE183904.rds
#
#   results/tables/NicheNet_top_TFs.csv
#
# External annotation:
#   GENCODE release 49
#   Long non-coding RNA annotation
#
# Workflow:
#   01. Load epithelial EMT Seurat object
#   02. Download/import GENCODE v49 lncRNA annotations
#   03. Identify lncRNAs present in scRNA-seq data
#   04. Load prioritized NicheNet TFs
#   05. Calculate TF-lncRNA Spearman correlations
#   06. Apply historical correlation criterion:
#         |rho| >= 0.5
#         p < 0.05
#   07. Calculate BH-adjusted correlation p-values
#   08. Cluster correlation-associated lncRNAs
#   09. Perform EMT_high vs EMT_low differential expression
#  10. Retain significant differentially expressed lncRNAs:
#         adjusted p < 0.05
#         |log2FC| > 0.25
#  11. Intersect DE and TF-correlation candidate sets
#  12. Generate high-confidence EMT-associated lncRNA list
#
# Main output for TCGA validation:
#   results/tables/
#   EMT_high_confidence_lncRNAs.csv
#
# IMPORTANT:
#   The TF-lncRNA correlation analysis reproduces the
#   historical cell-level analysis. Cells are therefore
#   treated as observations for this exploratory correlation
#   step. A sample-aware analysis can be added later as a
#   sensitivity/validation analysis.
############################################################


# ==========================================================
# 0. Required packages
# ==========================================================

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "dplyr",
  "tidyr",
  "tibble",
  "ggplot2",
  "Matrix",
  "matrixStats",
  "pheatmap",
  "rtracklayer",
  "R.utils"
)


missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]


if (length(missing_packages) > 0) {
  
  stop(
    paste0(
      "\nMissing required packages:\n",
      paste(
        missing_packages,
        collapse = ", "
      ),
      "\n\nInstall them before running this script.\n"
    ),
    call. = FALSE
  )
  
}


suppressPackageStartupMessages({
  
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(Matrix)
  library(matrixStats)
  library(pheatmap)
  library(rtracklayer)
  library(R.utils)
  
})


# ==========================================================
# 1. Reproducibility and analysis parameters
# ==========================================================

set.seed(123)


# ----------------------------------------------------------
# EMT populations
# ----------------------------------------------------------

emt_states_to_use <- c(
  "EMT_low",
  "EMT_intermediate",
  "EMT_high"
)


de_group_1 <- "EMT_high"

de_group_2 <- "EMT_low"


# ----------------------------------------------------------
# Historical TF-lncRNA correlation thresholds
# ----------------------------------------------------------

correlation_rho_cutoff <- 0.50

correlation_p_cutoff <- 0.05


# ----------------------------------------------------------
# Additional FDR threshold
#
# This does NOT replace the historical criterion.
# It is produced as a stricter sensitivity result.
# ----------------------------------------------------------

correlation_fdr_cutoff <- 0.05


# ----------------------------------------------------------
# EMT-high versus EMT-low differential expression
# ----------------------------------------------------------

de_min_pct <- 0.10

de_log2fc_cutoff <- 0.25

de_fdr_cutoff <- 0.05


# ----------------------------------------------------------
# lncRNA correlation clustering
#
# Historical workflow used k = 3.
# ----------------------------------------------------------

lncrna_k_clusters <- 3


# ----------------------------------------------------------
# Correlation calculation block size
#
# lncRNAs are processed in blocks to avoid creating one
# extremely large dense expression matrix.
# ----------------------------------------------------------

lncrna_block_size <- 250


# ----------------------------------------------------------
# Maximum number of high-confidence lncRNAs shown in
# violin plots.
# ----------------------------------------------------------

max_candidate_plot_genes <- 12


# ==========================================================
# 2. Project directories
# ==========================================================

base_dir <- "."


data_dir <- file.path(
  base_dir,
  "data"
)


gencode_dir <- file.path(
  data_dir,
  "gencode"
)


objects_dir <- file.path(
  base_dir,
  "objects"
)


results_dir <- file.path(
  base_dir,
  "results"
)


tables_dir <- file.path(
  results_dir,
  "tables"
)


plots_dir <- file.path(
  results_dir,
  "figures"
)


dir.create(
  gencode_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  objects_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  tables_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  plots_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==========================================================
# 3. Helper functions
# ==========================================================

log_message <- function(...) {
  
  cat(
    format(
      Sys.time(),
      "[%Y-%m-%d %H:%M:%S]"
    ),
    ...,
    "\n"
  )
  
}


checkpoint <- function(
    condition,
    message
) {
  
  if (!condition) {
    
    stop(
      message,
      call. = FALSE
    )
    
  }
  
}


# ----------------------------------------------------------
# Retrieve normalized RNA expression.
#
# Compatible with SeuratObject v4/v5.
# ----------------------------------------------------------

get_normalized_expression <- function(
    seurat_object
) {
  
  if (
    packageVersion("SeuratObject") >=
    package_version("5.0.0")
  ) {
    
    expression_matrix <- LayerData(
      object = seurat_object,
      assay = "RNA",
      layer = "data"
    )
    
  } else {
    
    expression_matrix <- GetAssayData(
      object = seurat_object,
      assay = "RNA",
      slot = "data"
    )
    
  }
  
  
  expression_matrix
  
}


# ----------------------------------------------------------
# Rank-transform each row and convert each ranked row to
# mean 0 / SD 1.
#
# Pearson correlation of these ranked vectors is equivalent
# to Spearman correlation.
# ----------------------------------------------------------

rank_and_scale_rows <- function(
    expression_matrix
) {
  
  dense_matrix <- as.matrix(
    expression_matrix
  )
  
  
  row_sd <- matrixStats::rowSds(
    dense_matrix
  )
  
  
  keep <- is.finite(
    row_sd
  ) &
    row_sd > 0
  
  
  dense_matrix <- dense_matrix[
    keep,
    ,
    drop = FALSE
  ]
  
  
  checkpoint(
    nrow(dense_matrix) > 0,
    "No non-constant genes remain for correlation analysis."
  )
  
  
  ranked_matrix <- t(
    apply(
      dense_matrix,
      1,
      rank,
      ties.method = "average"
    )
  )
  
  
  if (
    is.null(
      dim(
        ranked_matrix
      )
    )
  ) {
    
    ranked_matrix <- matrix(
      ranked_matrix,
      nrow = 1
    )
    
    
    rownames(ranked_matrix) <- rownames(
      dense_matrix
    )
    
  }
  
  
  ranked_mean <- rowMeans(
    ranked_matrix
  )
  
  
  ranked_sd <- matrixStats::rowSds(
    ranked_matrix
  )
  
  
  keep_rank <- is.finite(
    ranked_sd
  ) &
    ranked_sd > 0
  
  
  ranked_matrix <- ranked_matrix[
    keep_rank,
    ,
    drop = FALSE
  ]
  
  
  ranked_mean <- ranked_mean[
    keep_rank
  ]
  
  
  ranked_sd <- ranked_sd[
    keep_rank
  ]
  
  
  scaled_matrix <- sweep(
    ranked_matrix,
    1,
    ranked_mean,
    "-"
  )
  
  
  scaled_matrix <- sweep(
    scaled_matrix,
    1,
    ranked_sd,
    "/"
  )
  
  
  scaled_matrix
  
}


# ----------------------------------------------------------
# Approximate two-sided p-value for Spearman correlation
# using the t approximation employed in the original code.
# ----------------------------------------------------------

correlation_pvalue <- function(
    rho,
    n
) {
  
  denominator <- pmax(
    1e-12,
    1 - rho ^ 2
  )
  
  
  t_statistic <- rho *
    sqrt(
      (n - 2) /
        denominator
    )
  
  
  2 * pt(
    -abs(
      t_statistic
    ),
    df = n - 2
  )
  
}


# ==========================================================
# 4. Define input files
# ==========================================================

seurat_file <- file.path(
  objects_dir,
  "seurat_epithelial_EMT_GSE183904.rds"
)


tf_file <- file.path(
  tables_dir,
  "NicheNet_top_TFs.csv"
)


checkpoint(
  file.exists(seurat_file),
  paste0(
    "Epithelial EMT object not found:\n",
    seurat_file,
    "\n\nRun 04_emt_scoring.R first."
  )
)


checkpoint(
  file.exists(tf_file),
  paste0(
    "NicheNet TF table not found:\n",
    tf_file,
    "\n\nRun 06_nichenet_ligand_analysis.R first."
  )
)


# ==========================================================
# 5. Load epithelial EMT object
# ==========================================================

log_message(
  "Loading epithelial EMT Seurat object..."
)


seurat_epithelial <- readRDS(
  seurat_file
)


checkpoint(
  inherits(
    seurat_epithelial,
    "Seurat"
  ),
  "Input object is not a valid Seurat object."
)


checkpoint(
  "EMT_GMM_cluster" %in%
    colnames(
      seurat_epithelial@meta.data
    ),
  "EMT_GMM_cluster is missing."
)


cat(
  "\nEMT-state distribution:\n"
)


print(
  table(
    seurat_epithelial$EMT_GMM_cluster
  )
)


for (
  current_group
  in c(
    de_group_1,
    de_group_2
  )
) {
  
  checkpoint(
    current_group %in%
      unique(
        as.character(
          seurat_epithelial$EMT_GMM_cluster
        )
      ),
    paste0(
      "EMT group not found: ",
      current_group
    )
  )
  
}


DefaultAssay(
  seurat_epithelial
) <- "RNA"


# ==========================================================
# 6. Download GENCODE v49 lncRNA annotation
# ==========================================================

gencode_url <- paste0(
  "https://ftp.ebi.ac.uk/pub/databases/gencode/",
  "Gencode_human/release_49/",
  "gencode.v49.long_noncoding_RNAs.gtf.gz"
)


gencode_gtf_gz <- file.path(
  gencode_dir,
  "gencode.v49.long_noncoding_RNAs.gtf.gz"
)


gencode_gtf <- file.path(
  gencode_dir,
  "gencode.v49.long_noncoding_RNAs.gtf"
)


if (
  !file.exists(
    gencode_gtf_gz
  ) &&
  !file.exists(
    gencode_gtf
  )
) {
  
  log_message(
    "Downloading GENCODE v49 lncRNA annotation..."
  )
  
  
  download.file(
    
    url =
      gencode_url,
    
    destfile =
      gencode_gtf_gz,
    
    mode =
      "wb"
    
  )
  
}


checkpoint(
  file.exists(
    gencode_gtf_gz
  ) ||
    file.exists(
      gencode_gtf
    ),
  "GENCODE v49 lncRNA GTF could not be obtained."
)


# ==========================================================
# 7. Decompress GENCODE GTF when required
# ==========================================================

if (
  !file.exists(
    gencode_gtf
  )
) {
  
  log_message(
    "Decompressing GENCODE v49 GTF..."
  )
  
  
  R.utils::gunzip(
    
    filename =
      gencode_gtf_gz,
    
    destname =
      gencode_gtf,
    
    overwrite =
      FALSE,
    
    remove =
      FALSE
    
  )
  
}


checkpoint(
  file.exists(
    gencode_gtf
  ),
  "Decompressed GENCODE GTF not found."
)


# ==========================================================
# 8. Import GENCODE lncRNA annotations
# ==========================================================

log_message(
  "Importing GENCODE v49 lncRNA annotation..."
)


lnc_gtf <- rtracklayer::import(
  gencode_gtf
)


checkpoint(
  length(
    lnc_gtf
  ) > 0,
  "GENCODE GTF import returned no annotations."
)


# ==========================================================
# 9. Build gene-level lncRNA annotation table
# ==========================================================

gene_rows <- lnc_gtf[
  lnc_gtf$type ==
    "gene"
]


checkpoint(
  length(
    gene_rows
  ) > 0,
  "No gene-level entries detected in GENCODE lncRNA GTF."
)


gene_metadata_names <- names(
  S4Vectors::mcols(
    gene_rows
  )
)


checkpoint(
  "gene_name" %in%
    gene_metadata_names,
  "GENCODE gene_name field is missing."
)


checkpoint(
  "gene_id" %in%
    gene_metadata_names,
  "GENCODE gene_id field is missing."
)


if (
  "gene_type" %in%
  gene_metadata_names
) {
  
  gene_type_vector <- as.character(
    gene_rows$gene_type
  )
  
} else if (
  "gene_biotype" %in%
  gene_metadata_names
) {
  
  gene_type_vector <- as.character(
    gene_rows$gene_biotype
  )
  
} else {
  
  gene_type_vector <- rep(
    "long_noncoding_RNA",
    length(
      gene_rows
    )
  )
  
}


lnc_annotation <- data.frame(
  
  gene_name =
    as.character(
      gene_rows$gene_name
    ),
  
  gene_id_versioned =
    as.character(
      gene_rows$gene_id
    ),
  
  gene_id =
    sub(
      "\\..*$",
      "",
      as.character(
        gene_rows$gene_id
      )
    ),
  
  gene_type =
    gene_type_vector,
  
  stringsAsFactors =
    FALSE
  
)


lnc_annotation <- lnc_annotation[
  !is.na(
    lnc_annotation$gene_name
  ) &
    lnc_annotation$gene_name != "",
  ,
  drop = FALSE
]


lnc_annotation <- lnc_annotation %>%
  
  dplyr::distinct(
    gene_name,
    .keep_all = TRUE
  ) %>%
  
  dplyr::arrange(
    gene_name
  )


checkpoint(
  nrow(
    lnc_annotation
  ) > 0,
  "No valid lncRNA gene annotations were generated."
)


log_message(
  "GENCODE v49 lncRNA genes:",
  nrow(
    lnc_annotation
  )
)


write.csv(
  
  lnc_annotation,
  
  file.path(
    tables_dir,
    "GENCODE_v49_lncRNA_annotation_used.csv"
  ),
  
  row.names =
    FALSE
  
)


rm(
  lnc_gtf,
  gene_rows
)


gc()


# ==========================================================
# 10. Identify GENCODE lncRNAs present in GSE183904
# ==========================================================

dataset_genes <- rownames(
  seurat_epithelial
)


lnc_genes <- intersect(
  
  dataset_genes,
  
  lnc_annotation$gene_name
  
)


checkpoint(
  length(
    lnc_genes
  ) > 0,
  "No GENCODE v49 lncRNAs were detected in the dataset."
)


log_message(
  "GENCODE lncRNAs present in epithelial dataset:",
  length(
    lnc_genes
  )
)


lnc_present_table <- lnc_annotation %>%
  
  dplyr::filter(
    gene_name %in%
      lnc_genes
  )


write.csv(
  
  lnc_present_table,
  
  file.path(
    tables_dir,
    "GENCODE_v49_lncRNAs_present_GSE183904.csv"
  ),
  
  row.names =
    FALSE
  
)


############################################################
#
#                NICENET TF SELECTION
#
############################################################


# ==========================================================
# 11. Load NicheNet-prioritized transcription factors
# ==========================================================

tf_table <- read.csv(
  tf_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


checkpoint(
  "TF" %in%
    colnames(
      tf_table
    ),
  "NicheNet_top_TFs.csv does not contain a TF column."
)


tfs <- unique(
  as.character(
    tf_table$TF
  )
)


tfs <- tfs[
  !is.na(tfs) &
    tfs != ""
]


tfs_present <- intersect(
  tfs,
  dataset_genes
)


checkpoint(
  length(
    tfs_present
  ) > 0,
  "None of the NicheNet-prioritized TFs were found in the dataset."
)


log_message(
  "NicheNet-prioritized TFs:",
  length(
    tfs
  )
)


log_message(
  "NicheNet TFs present in expression matrix:",
  length(
    tfs_present
  )
)


write.csv(
  
  data.frame(
    TF =
      tfs_present
  ),
  
  file.path(
    tables_dir,
    "NicheNet_TFs_present_for_lncRNA_analysis.csv"
  ),
  
  row.names =
    FALSE
  
)


############################################################
#
#            TF-lncRNA CORRELATION ANALYSIS
#
############################################################


# ==========================================================
# 12. Restrict expression to EMT epithelial populations
# ==========================================================

cells_to_use <- colnames(
  seurat_epithelial
)[
  as.character(
    seurat_epithelial$EMT_GMM_cluster
  ) %in%
    emt_states_to_use
]


checkpoint(
  length(
    cells_to_use
  ) > 2,
  "Too few EMT epithelial cells for correlation analysis."
)


seurat_emt <- subset(
  seurat_epithelial,
  cells = cells_to_use
)


expression_matrix <- get_normalized_expression(
  seurat_emt
)


checkpoint(
  ncol(
    expression_matrix
  ) ==
    length(
      cells_to_use
    ),
  "Expression matrix/cell count mismatch."
)


n_cells <- ncol(
  expression_matrix
)


log_message(
  "Cells used for TF-lncRNA correlation:",
  n_cells
)


# ==========================================================
# 13. Prepare TF rank matrix
# ==========================================================

log_message(
  "Preparing TF expression ranks..."
)


tf_expression <- expression_matrix[
  tfs_present,
  ,
  drop = FALSE
]


tf_scaled <- rank_and_scale_rows(
  tf_expression
)


tfs_for_correlation <- rownames(
  tf_scaled
)


checkpoint(
  length(
    tfs_for_correlation
  ) > 0,
  "All TFs had zero expression variance."
)


log_message(
  "Non-constant TFs:",
  length(
    tfs_for_correlation
  )
)


# ==========================================================
# 14. Calculate TF-lncRNA Spearman correlations in blocks
# ==========================================================

log_message(
  "Calculating TF-lncRNA Spearman correlations..."
)


lncrna_blocks <- split(
  
  lnc_genes,
  
  ceiling(
    seq_along(
      lnc_genes
    ) /
      lncrna_block_size
  )
  
)


correlation_results <- vector(
  "list",
  length(
    lncrna_blocks
  )
)


for (
  block_index
  in seq_along(
    lncrna_blocks
  )
) {
  
  block_genes <- lncrna_blocks[
    [block_index]
  ]
  
  
  log_message(
    paste0(
      "Correlation block ",
      block_index,
      "/",
      length(
        lncrna_blocks
      ),
      " (",
      length(
        block_genes
      ),
      " lncRNAs)"
    )
  )
  
  
  lnc_expression_block <- expression_matrix[
    block_genes,
    ,
    drop = FALSE
  ]
  
  
  lnc_scaled <- try(
    
    rank_and_scale_rows(
      lnc_expression_block
    ),
    
    silent = TRUE
    
  )
  
  
  if (
    inherits(
      lnc_scaled,
      "try-error"
    )
  ) {
    
    log_message(
      "Skipping block with no variable lncRNAs:",
      block_index
    )
    
    next
    
  }
  
  
  if (
    nrow(
      lnc_scaled
    ) == 0
  ) {
    
    next
    
  }
  
  
  # --------------------------------------------------------
  # Pearson correlation of standardized ranks =
  # Spearman rank correlation
  # --------------------------------------------------------
  
  correlation_matrix <- (
    
    tf_scaled %*%
      t(
        lnc_scaled
      )
    
  ) / (
    n_cells - 1
  )
  
  
  # Clamp tiny floating-point excursions beyond [-1, 1]
  correlation_matrix[
    correlation_matrix > 1
  ] <- 1
  
  
  correlation_matrix[
    correlation_matrix < -1
  ] <- -1
  
  
  correlation_df <- data.frame(
    
    TF =
      rep(
        rownames(
          correlation_matrix
        ),
        times =
          ncol(
            correlation_matrix
          )
      ),
    
    lncRNA =
      rep(
        colnames(
          correlation_matrix
        ),
        each =
          nrow(
            correlation_matrix
          )
      ),
    
    rho =
      as.vector(
        correlation_matrix
      ),
    
    stringsAsFactors =
      FALSE
    
  )
  
  
  correlation_df$p_value <-
    
    correlation_pvalue(
      
      rho =
        correlation_df$rho,
      
      n =
        n_cells
      
    )
  
  
  correlation_results[
    [block_index]
  ] <- correlation_df
  
  
  rm(
    lnc_expression_block,
    lnc_scaled,
    correlation_matrix,
    correlation_df
  )
  
  
  gc()
  
}


correlation_results <- correlation_results[
  !vapply(
    correlation_results,
    is.null,
    FUN.VALUE = logical(1)
  )
]


checkpoint(
  length(
    correlation_results
  ) > 0,
  "No TF-lncRNA correlations were generated."
)


tf_lnc_correlations <- dplyr::bind_rows(
  correlation_results
)


rm(
  correlation_results
)


gc()


checkpoint(
  nrow(
    tf_lnc_correlations
  ) > 0,
  "TF-lncRNA correlation table is empty."
)


# ==========================================================
# 15. Multiple-testing correction
# ==========================================================

tf_lnc_correlations$p_adj_BH <- p.adjust(
  
  tf_lnc_correlations$p_value,
  
  method =
    "BH"
  
)


tf_lnc_correlations$abs_rho <- abs(
  tf_lnc_correlations$rho
)


tf_lnc_correlations <- tf_lnc_correlations %>%
  
  dplyr::arrange(
    dplyr::desc(
      abs_rho
    ),
    p_value
  )


log_message(
  "Total TF-lncRNA pairs tested:",
  nrow(
    tf_lnc_correlations
  )
)


# ==========================================================
# 16. Save complete correlation result as RDS
#
# This can be large, so we do not write the entire matrix as
# a CSV intended for GitHub.
# ==========================================================

all_correlations_file <- file.path(
  
  objects_dir,
  
  "TF_lncRNA_correlations_all_GSE183904.rds"
  
)


saveRDS(
  
  tf_lnc_correlations,
  
  all_correlations_file
  
)


# ==========================================================
# 17. Historical correlation criterion
#
# Reproduces original analysis:
#   |rho| >= 0.5
#   raw p < 0.05
# ==========================================================

significant_correlations_historical <-
  
  tf_lnc_correlations %>%
  
  dplyr::filter(
    
    abs_rho >=
      correlation_rho_cutoff,
    
    p_value <
      correlation_p_cutoff
    
  )


checkpoint(
  nrow(
    significant_correlations_historical
  ) > 0,
  paste0(
    "No TF-lncRNA pairs passed the historical criterion: ",
    "|rho| >= ",
    correlation_rho_cutoff,
    " and p < ",
    correlation_p_cutoff,
    "."
  )
)


log_message(
  "Pairs passing historical correlation criterion:",
  nrow(
    significant_correlations_historical
  )
)


write.csv(
  
  significant_correlations_historical,
  
  file.path(
    tables_dir,
    "TF_lncRNA_pairs_historical_criterion.csv"
  ),
  
  row.names =
    FALSE
  
)


# ==========================================================
# 18. FDR-controlled correlation result
#
# This is an additional stricter sensitivity analysis.
# ==========================================================

significant_correlations_fdr <-
  
  tf_lnc_correlations %>%
  
  dplyr::filter(
    
    abs_rho >=
      correlation_rho_cutoff,
    
    p_adj_BH <
      correlation_fdr_cutoff
    
  )


log_message(
  "Pairs passing |rho| and BH-FDR criterion:",
  nrow(
    significant_correlations_fdr
  )
)


write.csv(
  
  significant_correlations_fdr,
  
  file.path(
    tables_dir,
    "TF_lncRNA_pairs_BH_FDR05.csv"
  ),
  
  row.names =
    FALSE
  
)


# ==========================================================
# 19. Identify correlation-associated lncRNAs
# ==========================================================

correlation_lncRNAs <- unique(
  
  significant_correlations_historical$lncRNA
  
)


log_message(
  "lncRNAs associated with prioritized TFs:",
  length(
    correlation_lncRNAs
  )
)


write.csv(
  
  data.frame(
    lncRNA =
      correlation_lncRNAs
  ),
  
  file.path(
    tables_dir,
    "TF_associated_lncRNAs.csv"
  ),
  
  row.names =
    FALSE
  
)


############################################################
#
#         lncRNA CORRELATION MODULES / CLUSTERS
#
############################################################


# ==========================================================
# 20. lncRNA-lncRNA correlation matrix
# ==========================================================

lncrna_cluster_table <- NULL


if (
  length(
    correlation_lncRNAs
  ) >= 2
) {
  
  log_message(
    "Calculating lncRNA-lncRNA correlation matrix..."
  )
  
  
  lnc_candidate_expression <- as.matrix(
    
    expression_matrix[
      correlation_lncRNAs,
      ,
      drop = FALSE
    ]
    
  )
  
  
  candidate_sd <- matrixStats::rowSds(
    lnc_candidate_expression
  )
  
  
  keep_candidate <- is.finite(
    candidate_sd
  ) &
    candidate_sd > 0
  
  
  lnc_candidate_expression <-
    lnc_candidate_expression[
      keep_candidate,
      ,
      drop = FALSE
    ]
  
  
  checkpoint(
    nrow(
      lnc_candidate_expression
    ) >= 2,
    "Too few variable correlation-associated lncRNAs."
  )
  
  
  lnc_correlation_matrix <- suppressWarnings(
    
    cor(
      
      t(
        lnc_candidate_expression
      ),
      
      method =
        "spearman",
      
      use =
        "pairwise.complete.obs"
      
    )
    
  )
  
  
  lnc_correlation_matrix[
    !is.finite(
      lnc_correlation_matrix
    )
  ] <- 0
  
  
  diag(
    lnc_correlation_matrix
  ) <- 1
  
  
  # ========================================================
  # 21. K-means correlation-profile clustering
  #
  # Historical workflow used k = 3.
  # ========================================================
  
  unique_profiles <- nrow(
    
    unique(
      
      round(
        lnc_correlation_matrix,
        digits = 8
      )
      
    )
    
  )
  
  
  actual_k <- min(
    
    lncrna_k_clusters,
    
    nrow(
      lnc_correlation_matrix
    ),
    
    unique_profiles
    
  )
  
  
  if (
    actual_k >= 2
  ) {
    
    set.seed(123)
    
    
    km <- kmeans(
      
      lnc_correlation_matrix,
      
      centers =
        actual_k,
      
      nstart =
        50
      
    )
    
    
    lncrna_cluster_table <- data.frame(
      
      lncRNA =
        names(
          km$cluster
        ),
      
      cluster =
        as.integer(
          km$cluster
        ),
      
      stringsAsFactors =
        FALSE
      
    )
    
  } else {
    
    lncrna_cluster_table <- data.frame(
      
      lncRNA =
        rownames(
          lnc_correlation_matrix
        ),
      
      cluster =
        1L,
      
      stringsAsFactors =
        FALSE
      
    )
    
  }
  
  
  lncrna_cluster_table <- lncrna_cluster_table %>%
    
    dplyr::arrange(
      cluster,
      lncRNA
    )
  
  
  write.csv(
    
    lncrna_cluster_table,
    
    file.path(
      tables_dir,
      "lncRNA_correlation_clusters.csv"
    ),
    
    row.names =
      FALSE
    
  )
  
  
  # ========================================================
  # 22. lncRNA-lncRNA correlation heatmap
  # ========================================================
  
  heatmap_annotation <- data.frame(
    
    cluster =
      factor(
        lncrna_cluster_table$cluster
      )
    
  )
  
  
  rownames(
    heatmap_annotation
  ) <- lncrna_cluster_table$lncRNA
  
  
  heatmap_annotation <- heatmap_annotation[
    rownames(
      lnc_correlation_matrix
    ),
    ,
    drop = FALSE
  ]
  
  
  pheatmap::pheatmap(
    
    lnc_correlation_matrix,
    
    annotation_row =
      heatmap_annotation,
    
    annotation_col =
      heatmap_annotation,
    
    show_rownames =
      nrow(
        lnc_correlation_matrix
      ) <= 50,
    
    show_colnames =
      ncol(
        lnc_correlation_matrix
      ) <= 50,
    
    border_color =
      NA,
    
    main =
      "Correlation-Associated lncRNA Modules",
    
    filename = file.path(
      plots_dir,
      "lncRNA_01_correlation_modules_heatmap.pdf"
    ),
    
    width =
      10,
    
    height =
      9
    
  )
  
}


# ==========================================================
# 23. TF-lncRNA correlation heatmap
# ==========================================================

sig_heatmap_df <-
  
  significant_correlations_historical %>%
  
  dplyr::select(
    TF,
    lncRNA,
    rho
  )


tf_lnc_heatmap_df <- sig_heatmap_df %>%
  
  tidyr::pivot_wider(
    
    names_from =
      lncRNA,
    
    values_from =
      rho,
    
    values_fill =
      0
    
  )


if (
  nrow(
    tf_lnc_heatmap_df
  ) > 0 &&
  ncol(
    tf_lnc_heatmap_df
  ) > 1
) {
  
  tf_names <- tf_lnc_heatmap_df$TF
  
  
  tf_lnc_heatmap_matrix <- as.matrix(
    
    tf_lnc_heatmap_df[
      ,
      setdiff(
        colnames(
          tf_lnc_heatmap_df
        ),
        "TF"
      ),
      drop = FALSE
    ]
    
  )
  
  
  rownames(
    tf_lnc_heatmap_matrix
  ) <- tf_names
  
  
  pheatmap::pheatmap(
    
    tf_lnc_heatmap_matrix,
    
    cluster_rows =
      TRUE,
    
    cluster_cols =
      TRUE,
    
    border_color =
      NA,
    
    main =
      "NicheNet TF-lncRNA Spearman Correlations",
    
    filename = file.path(
      plots_dir,
      "lncRNA_02_TF_lncRNA_correlation_heatmap.pdf"
    ),
    
    width =
      12,
    
    height =
      9
    
  )
  
}


############################################################
#
#          EMT-HIGH vs EMT-LOW DIFFERENTIAL EXPRESSION
#
############################################################


# ==========================================================
# 24. Set EMT state as Seurat identity
# ==========================================================

Idents(
  seurat_epithelial
) <- seurat_epithelial$EMT_GMM_cluster


cat(
  "\nCells used for EMT-high vs EMT-low differential expression:\n"
)


print(
  table(
    Idents(
      seurat_epithelial
    )[
      Idents(
        seurat_epithelial
      ) %in%
        c(
          de_group_1,
          de_group_2
        )
    ]
  )
)


# ==========================================================
# 25. Run Wilcoxon differential expression
# ==========================================================

log_message(
  paste0(
    "Running differential expression: ",
    de_group_1,
    " vs ",
    de_group_2,
    "..."
  )
)


deg_emthigh_vs_low <- FindMarkers(
  
  object =
    seurat_epithelial,
  
  ident.1 =
    de_group_1,
  
  ident.2 =
    de_group_2,
  
  assay =
    "RNA",
  
  logfc.threshold =
    0,
  
  min.pct =
    de_min_pct,
  
  test.use =
    "wilcox",
  
  only.pos =
    FALSE,
  
  verbose =
    TRUE
  
)


checkpoint(
  nrow(
    deg_emthigh_vs_low
  ) > 0,
  "FindMarkers returned no genes."
)


deg_emthigh_vs_low$gene <- rownames(
  deg_emthigh_vs_low
)


# ==========================================================
# 26. Standardize fold-change column
# ==========================================================

fc_column_candidates <- c(
  "avg_log2FC",
  "avg_logFC"
)


fc_column <- intersect(
  
  fc_column_candidates,
  
  colnames(
    deg_emthigh_vs_low
  )
  
)[1]


checkpoint(
  !is.na(
    fc_column
  ),
  "Could not identify average log fold-change column."
)


deg_emthigh_vs_low$avg_log2FC_standardized <-
  
  deg_emthigh_vs_low[
    [fc_column]
  ]


deg_emthigh_vs_low <- deg_emthigh_vs_low[
  ,
  c(
    "gene",
    setdiff(
      colnames(
        deg_emthigh_vs_low
      ),
      "gene"
    )
  ),
  drop = FALSE
]


write.csv(
  
  deg_emthigh_vs_low,
  
  file.path(
    tables_dir,
    "EMT_high_vs_low_all_genes.csv"
  ),
  
  row.names =
    FALSE
  
)


# ==========================================================
# 27. Significant DE genes
#
# Historical code:
#   p_val_adj < 0.05
#   |avg_log2FC| > 0.25
# ==========================================================

deg_significant <- deg_emthigh_vs_low %>%
  
  dplyr::filter(
    
    !is.na(
      p_val_adj
    ),
    
    p_val_adj <
      de_fdr_cutoff,
    
    abs(
      avg_log2FC_standardized
    ) >
      de_log2fc_cutoff
    
  )


log_message(
  "Significant EMT-high vs EMT-low genes:",
  nrow(
    deg_significant
  )
)


# ==========================================================
# 28. Restrict DE results to GENCODE v49 lncRNAs
# ==========================================================

deg_lnc_all <- deg_emthigh_vs_low %>%
  
  dplyr::filter(
    gene %in%
      lnc_genes
  ) %>%
  
  dplyr::left_join(
    
    lnc_annotation,
    
    by = c(
      "gene" =
        "gene_name"
    )
    
  )


deg_lnc_significant <- deg_significant %>%
  
  dplyr::filter(
    gene %in%
      lnc_genes
  ) %>%
  
  dplyr::mutate(
    
    direction =
      dplyr::case_when(
        
        avg_log2FC_standardized > 0 ~
          "Higher_in_EMT_high",
        
        avg_log2FC_standardized < 0 ~
          "Lower_in_EMT_high",
        
        TRUE ~
          "No_change"
        
      )
    
  ) %>%
  
  dplyr::left_join(
    
    lnc_annotation,
    
    by = c(
      "gene" =
        "gene_name"
    )
    
  ) %>%
  
  dplyr::arrange(
    p_val_adj
  )


log_message(
  "lncRNAs tested:",
  nrow(
    deg_lnc_all
  )
)


log_message(
  "Significant DE lncRNAs:",
  nrow(
    deg_lnc_significant
  )
)


write.csv(
  
  deg_lnc_all,
  
  file.path(
    tables_dir,
    "EMT_high_vs_low_all_lncRNAs.csv"
  ),
  
  row.names =
    FALSE
  
)


write.csv(
  
  deg_lnc_significant,
  
  file.path(
    tables_dir,
    "EMT_high_vs_low_significant_lncRNAs.csv"
  ),
  
  row.names =
    FALSE
  
)


# ==========================================================
# 29. lncRNA differential-expression volcano plot
# ==========================================================

volcano_df <- deg_lnc_all


volcano_df$neg_log10_FDR <-
  
  -log10(
    
    pmax(
      volcano_df$p_val_adj,
      .Machine$double.xmin
    )
    
  )


volcano_df$significant <- with(
  
  volcano_df,
  
  !is.na(
    p_val_adj
  ) &
    
    p_val_adj <
    de_fdr_cutoff &
    
    abs(
      avg_log2FC_standardized
    ) >
    de_log2fc_cutoff
  
)


p_volcano <- ggplot(
  
  volcano_df,
  
  aes(
    x =
      avg_log2FC_standardized,
    
    y =
      neg_log10_FDR
  )
  
) +
  
  geom_point(
    
    aes(
      shape =
        significant
    ),
    
    alpha =
      0.6,
    
    size =
      1.3
    
  ) +
  
  geom_vline(
    
    xintercept =
      c(
        -de_log2fc_cutoff,
        de_log2fc_cutoff
      ),
    
    linetype =
      "dashed"
    
  ) +
  
  geom_hline(
    
    yintercept =
      -log10(
        de_fdr_cutoff
      ),
    
    linetype =
      "dashed"
    
  ) +
  
  labs(
    
    title =
      "EMT-high vs EMT-low lncRNA Differential Expression",
    
    x =
      "Average log2 fold change",
    
    y =
      "-log10 adjusted p-value",
    
    shape =
      "Significant"
    
  ) +
  
  theme_classic()


ggsave(
  
  filename = file.path(
    plots_dir,
    "lncRNA_03_EMT_high_vs_low_volcano.pdf"
  ),
  
  plot =
    p_volcano,
  
  width =
    8,
  
  height =
    6
  
)


############################################################
#
#         HIGH-CONFIDENCE EMT-ASSOCIATED lncRNAs
#
############################################################


# ==========================================================
# 30. Best TF correlation for each candidate lncRNA
# ==========================================================

best_tf_per_lnc <- significant_correlations_historical %>%
  
  dplyr::group_by(
    lncRNA
  ) %>%
  
  dplyr::slice_max(
    
    order_by =
      abs_rho,
    
    n =
      1,
    
    with_ties =
      FALSE
    
  ) %>%
  
  dplyr::ungroup() %>%
  
  dplyr::rename(
    
    best_TF =
      TF,
    
    best_TF_rho =
      rho,
    
    best_TF_p =
      p_value,
    
    best_TF_q_BH =
      p_adj_BH
    
  ) %>%
  
  dplyr::select(
    
    lncRNA,
    
    best_TF,
    
    best_TF_rho,
    
    best_TF_p,
    
    best_TF_q_BH
    
  )


# ==========================================================
# 31. Intersect significant DE lncRNAs with TF-associated
#     lncRNAs
# ==========================================================

high_confidence_lncRNAs <- deg_lnc_significant %>%
  
  dplyr::filter(
    
    gene %in%
      correlation_lncRNAs
    
  ) %>%
  
  dplyr::left_join(
    
    best_tf_per_lnc,
    
    by = c(
      "gene" =
        "lncRNA"
    )
    
  )


# Add correlation-cluster assignment where available
if (
  !is.null(
    lncrna_cluster_table
  )
) {
  
  high_confidence_lncRNAs <-
    
    high_confidence_lncRNAs %>%
    
    dplyr::left_join(
      
      lncrna_cluster_table,
      
      by = c(
        "gene" =
          "lncRNA"
      )
      
    )
  
} else {
  
  high_confidence_lncRNAs$cluster <- NA_integer_
  
}


high_confidence_lncRNAs <-
  
  high_confidence_lncRNAs %>%
  
  dplyr::mutate(
    
    selection_basis =
      paste0(
        "DE: BH<",
        de_fdr_cutoff,
        ", |log2FC|>",
        de_log2fc_cutoff,
        "; TF-correlation: |rho|>=",
        correlation_rho_cutoff,
        ", raw p<",
        correlation_p_cutoff
      )
    
  ) %>%
  
  dplyr::arrange(
    
    p_val_adj,
    
    dplyr::desc(
      abs(
        avg_log2FC_standardized
      )
    )
    
  )


checkpoint(
  nrow(
    high_confidence_lncRNAs
  ) > 0,
  paste0(
    "No lncRNAs were shared between significant EMT DE ",
    "and TF-correlation candidate sets."
  )
)


log_message(
  "High-confidence EMT-associated lncRNAs:",
  nrow(
    high_confidence_lncRNAs
  )
)


write.csv(
  
  high_confidence_lncRNAs,
  
  file.path(
    tables_dir,
    "EMT_high_confidence_lncRNAs.csv"
  ),
  
  row.names =
    FALSE
  
)


# ==========================================================
# 32. Simple candidate gene list for downstream TCGA
# ==========================================================

tcga_candidate_list <- data.frame(
  
  gene =
    unique(
      high_confidence_lncRNAs$gene
    )
  
)


write.csv(
  
  tcga_candidate_list,
  
  file.path(
    tables_dir,
    "EMT_lncRNA_candidates_for_TCGA.csv"
  ),
  
  row.names =
    FALSE
  
)


cat(
  "\nHigh-confidence EMT-associated lncRNAs:\n"
)


print(
  high_confidence_lncRNAs[
    ,
    intersect(
      c(
        "gene",
        "gene_type",
        "avg_log2FC_standardized",
        "p_val_adj",
        "direction",
        "best_TF",
        "best_TF_rho",
        "best_TF_p",
        "best_TF_q_BH",
        "cluster"
      ),
      colnames(
        high_confidence_lncRNAs
      )
    ),
    drop = FALSE
  ]
)


# ==========================================================
# 33. Stricter BH-FDR correlation sensitivity candidates
# ==========================================================

fdr_correlation_lncRNAs <- unique(
  significant_correlations_fdr$lncRNA
)


high_confidence_fdr_sensitivity <-
  
  deg_lnc_significant %>%
  
  dplyr::filter(
    
    gene %in%
      fdr_correlation_lncRNAs
    
  ) %>%
  
  dplyr::arrange(
    p_val_adj
  )


write.csv(
  
  high_confidence_fdr_sensitivity,
  
  file.path(
    tables_dir,
    "EMT_high_confidence_lncRNAs_BH_correlation_sensitivity.csv"
  ),
  
  row.names =
    FALSE
  
)


log_message(
  "High-confidence candidates also surviving TF-correlation BH-FDR:",
  nrow(
    high_confidence_fdr_sensitivity
  )
)


# ==========================================================
# 34. Candidate expression violin plots
# ==========================================================

candidate_genes <- unique(
  high_confidence_lncRNAs$gene
)


candidate_genes <- intersect(
  candidate_genes,
  rownames(
    seurat_epithelial
  )
)


if (
  length(
    candidate_genes
  ) > 0
) {
  
  genes_to_plot <- head(
    
    candidate_genes,
    
    max_candidate_plot_genes
    
  )
  
  
  p_candidate_violin <- VlnPlot(
    
    seurat_epithelial,
    
    features =
      genes_to_plot,
    
    group.by =
      "EMT_GMM_cluster",
    
    pt.size =
      0,
    
    ncol =
      3
    
  )
  
  
  ggsave(
    
    filename = file.path(
      plots_dir,
      "lncRNA_04_high_confidence_candidates_violin.pdf"
    ),
    
    plot =
      p_candidate_violin,
    
    width =
      12,
    
    height =
      max(
        6,
        ceiling(
          length(
            genes_to_plot
          ) /
            3
        ) * 3
      )
    
  )
  
}


# ==========================================================
# 35. Candidate FeaturePlots
# ==========================================================

if (
  length(
    candidate_genes
  ) > 0
) {
  
  genes_to_plot <- head(
    
    candidate_genes,
    
    min(
      9,
      length(
        candidate_genes
      )
    )
    
  )
  
  
  p_candidate_feature <- FeaturePlot(
    
    seurat_epithelial,
    
    features =
      genes_to_plot,
    
    reduction =
      "umap",
    
    raster =
      TRUE,
    
    ncol =
      3
    
  )
  
  
  ggsave(
    
    filename = file.path(
      plots_dir,
      "lncRNA_05_high_confidence_candidates_featureplot.pdf"
    ),
    
    plot =
      p_candidate_feature,
    
    width =
      12,
    
    height =
      max(
        5,
        ceiling(
          length(
            genes_to_plot
          ) /
            3
        ) * 4
      )
    
  )
  
}


############################################################
#
#                   SAVE ANALYSIS OBJECT
#
############################################################


# ==========================================================
# 36. Save reusable lncRNA-analysis object
# ==========================================================

lncrna_analysis <- list(
  
  gencode_release =
    49,
  
  lnc_annotation =
    lnc_annotation,
  
  lnc_genes_detected =
    lnc_genes,
  
  NicheNet_TFs =
    tfs_present,
  
  correlation_threshold =
    correlation_rho_cutoff,
  
  historical_p_threshold =
    correlation_p_cutoff,
  
  BH_FDR_threshold =
    correlation_fdr_cutoff,
  
  significant_correlations_historical =
    significant_correlations_historical,
  
  significant_correlations_FDR =
    significant_correlations_fdr,
  
  correlation_clusters =
    lncrna_cluster_table,
  
  differential_expression_all =
    deg_emthigh_vs_low,
  
  differential_expression_lncRNAs =
    deg_lnc_significant,
  
  high_confidence_lncRNAs =
    high_confidence_lncRNAs,
  
  high_confidence_FDR_sensitivity =
    high_confidence_fdr_sensitivity
  
)


lncrna_object_file <- file.path(
  
  objects_dir,
  
  "EMT_lncRNA_analysis_GSE183904.rds"
  
)


saveRDS(
  
  lncrna_analysis,
  
  lncrna_object_file
  
)


# ==========================================================
# 37. Save session information
# ==========================================================

capture.output(
  
  sessionInfo(),
  
  file = file.path(
    results_dir,
    "sessionInfo_07_lncRNA_analysis.txt"
  )
  
)


# ==========================================================
# 38. Completion summary
# ==========================================================

cat(
  
  "\n",
  "=====================================================\n",
  " EMT lncRNA ANALYSIS COMPLETED\n",
  "=====================================================\n",
  "\n",
  "GENCODE release:\n",
  "49",
  "\n\n",
  "GENCODE lncRNAs present in dataset:\n",
  length(
    lnc_genes
  ),
  "\n\n",
  "NicheNet TFs used:\n",
  length(
    tfs_for_correlation
  ),
  "\n\n",
  "TF-lncRNA pairs tested:\n",
  nrow(
    tf_lnc_correlations
  ),
  "\n\n",
  "Pairs passing historical |rho|/p criterion:\n",
  nrow(
    significant_correlations_historical
  ),
  "\n\n",
  "Pairs passing |rho|/BH-FDR criterion:\n",
  nrow(
    significant_correlations_fdr
  ),
  "\n\n",
  "Significant EMT-high vs EMT-low lncRNAs:\n",
  nrow(
    deg_lnc_significant
  ),
  "\n\n",
  "High-confidence EMT-associated lncRNAs:\n",
  nrow(
    high_confidence_lncRNAs
  ),
  "\n\n",
  "TCGA candidate file:\n",
  file.path(
    tables_dir,
    "EMT_lncRNA_candidates_for_TCGA.csv"
  ),
  "\n\n",
  "Analysis object:\n",
  lncrna_object_file,
  "\n",
  "=====================================================\n"
  
)


log_message(
  "07_emt_lncrna_analysis.R completed."
)