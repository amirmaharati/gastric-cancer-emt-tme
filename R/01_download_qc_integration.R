############################################################
# GSE183904 scRNA-seq Preprocessing, Doublet Removal,
# Harmony Integration, and Clustering
#
# Project:
#   Gastric Cancer EMT–TME Analysis
#
# Dataset:
#   GEO accession: GSE183904
#
# Workflow:
#   01. Download GSE183904 raw count matrices
#   02. Construct combined Seurat object
#   03. Calculate quality-control metrics
#   04. Apply QC filtering
#   05. Split cells by biological sample
#   06. Run DoubletFinder independently per sample
#   07. Remove predicted doublets
#   08. Merge singlet-only samples
#   09. Normalize merged object
#   10. Identify variable features
#   11. Scale data
#   12. PCA
#   13. Harmony batch correction
#   14. UMAP
#   15. Graph-based clustering
#
# Main output:
#   objects/seurat_integrated_GSE183904.rds
#
# Additional outputs:
#   objects/seurat_raw_GSE183904.rds
#   objects/seurat_QC_GSE183904.rds
#   objects/seurat_singlets_GSE183904.rds
#
#   results/figures/
#   results/tables/
#
# IMPORTANT:
#   The historical project code confirmed that DoubletFinder
#   was part of the workflow, but the actual DoubletFinder
#   function calls were not preserved in the supplied script.
#   The DoubletFinder implementation below reconstructs this
#   missing step using the current package workflow.
############################################################


# ==========================================================
# 0. Required packages
# ==========================================================

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "data.table",
  "Matrix",
  "ggplot2",
  "harmony",
  "DoubletFinder"
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
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(harmony)
  library(DoubletFinder)
  
})


# ==========================================================
# 1. Global parameters
# ==========================================================

set.seed(123)


# Number of variable genes
n_variable_features <- 2000


# Global clustering resolution
global_cluster_resolution <- 1.0


# Preliminary clustering resolution used only for
# estimating homotypic doublets within each sample
doublet_cluster_resolution <- 0.5


# DoubletFinder pN
doublet_pN <- 0.25


# ----------------------------------------------------------
# IMPORTANT
# ----------------------------------------------------------
#
# The exact expected doublet rate used in the historical
# gastric-cancer analysis was not present in the supplied
# source code.
#
# 0.075 is therefore currently an explicit reconstruction
# parameter, NOT a recovered historical parameter.
#
# Change this later if the original loading-specific rate
# is recovered.
#
# ----------------------------------------------------------

doublet_rate <- 0.075


# Minimum number of cells required to attempt
# sample-level DoubletFinder analysis
minimum_cells_for_doubletfinder <- 100


# ==========================================================
# 2. Define project directories
# ==========================================================

base_dir <- "."


data_dir <- file.path(
  base_dir,
  "data"
)


raw_dir <- file.path(
  data_dir,
  "GSE183904_RAW"
)


objects_dir <- file.path(
  base_dir,
  "objects"
)


results_dir <- file.path(
  base_dir,
  "results"
)


plots_dir <- file.path(
  results_dir,
  "figures"
)


tables_dir <- file.path(
  results_dir,
  "tables"
)


doublet_tables_dir <- file.path(
  tables_dir,
  "doubletfinder"
)


dir.create(
  data_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  raw_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  objects_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  plots_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  tables_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  doublet_tables_dir,
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
# Function to choose PCs explaining >=90% of the variance
# represented by the computed PCA dimensions
# ----------------------------------------------------------

choose_pcs <- function(
    seurat_object,
    reduction = "pca",
    threshold = 0.90
) {
  
  pca_sd <- Stdev(
    seurat_object,
    reduction = reduction
  )
  
  
  variance_explained <- (
    pca_sd ^ 2
  ) / sum(
    pca_sd ^ 2
  )
  
  
  cumulative_variance <- cumsum(
    variance_explained
  )
  
  
  selected_pc <- which(
    cumulative_variance >= threshold
  )[1]
  
  
  if (is.na(selected_pc)) {
    
    selected_pc <- length(
      pca_sd
    )
    
  }
  
  
  selected_pc
}


# ==========================================================
# 4. Download GSE183904 raw data
# ==========================================================

geo_url <- paste0(
  "https://www.ncbi.nlm.nih.gov/geo/download/",
  "?acc=GSE183904&format=file"
)


tar_file <- file.path(
  data_dir,
  "GSE183904_RAW.tar"
)


if (!file.exists(tar_file)) {
  
  log_message(
    "Downloading GSE183904 raw archive..."
  )
  
  
  download.file(
    url = geo_url,
    destfile = tar_file,
    mode = "wb"
  )
  
} else {
  
  log_message(
    "GSE183904_RAW.tar already exists.",
    "Skipping download."
  )
  
}


checkpoint(
  file.exists(tar_file),
  paste0(
    "Download failed. File not found:\n",
    tar_file
  )
)


# ==========================================================
# 5. Extract GEO archive
# ==========================================================

raw_files <- list.files(
  raw_dir,
  pattern = "\\.csv\\.gz$",
  full.names = TRUE
)


if (length(raw_files) == 0) {
  
  log_message(
    "Extracting GSE183904_RAW.tar..."
  )
  
  
  untar(
    tar_file,
    exdir = raw_dir
  )
  
  
  raw_files <- list.files(
    raw_dir,
    pattern = "\\.csv\\.gz$",
    full.names = TRUE
  )
  
}


checkpoint(
  length(raw_files) > 0,
  paste0(
    "No .csv.gz count matrices found in:\n",
    raw_dir
  )
)


raw_files <- sort(
  raw_files
)


log_message(
  "Number of raw sample files:",
  length(raw_files)
)


# ==========================================================
# 6. Read raw count matrices
# ==========================================================

log_message(
  "Reading individual sample matrices..."
)


matrix_list <- vector(
  mode = "list",
  length = length(raw_files)
)


names(matrix_list) <- basename(
  raw_files
)


for (i in seq_along(raw_files)) {
  
  current_file <- raw_files[i]
  
  
  sample_id <- sub(
    "_.*",
    "",
    basename(current_file)
  )
  
  
  log_message(
    paste0(
      "Reading sample ",
      i,
      "/",
      length(raw_files),
      ": ",
      sample_id
    )
  )
  
  
  expression_df <- fread(
    current_file
  )
  
  
  checkpoint(
    ncol(expression_df) > 1,
    paste(
      "No expression columns detected in:",
      basename(current_file)
    )
  )
  
  
  gene_names <- expression_df[[1]]
  
  
  expression_df[[1]] <- NULL
  
  
  # --------------------------------------------------------
  # Prefix every cell barcode with its sample identifier.
  #
  # This preserves globally unique cell names after samples
  # are combined.
  # --------------------------------------------------------
  
  colnames(expression_df) <- paste0(
    sample_id,
    "_",
    colnames(expression_df)
  )
  
  
  expression_matrix <- as(
    as.matrix(expression_df),
    "dgCMatrix"
  )
  
  
  rownames(expression_matrix) <- gene_names
  
  
  matrix_list[[i]] <- expression_matrix
  
  
  rm(
    expression_df,
    expression_matrix
  )
  
  
  gc()
  
}


# ==========================================================
# 7. Combine count matrices
# ==========================================================

log_message(
  "Combining raw matrices..."
)


count_matrix <- do.call(
  cbind,
  matrix_list
)


rm(
  matrix_list
)


gc()


checkpoint(
  nrow(count_matrix) > 0 &&
    ncol(count_matrix) > 0,
  "Combined expression matrix is empty."
)


log_message(
  "Combined matrix dimensions:",
  nrow(count_matrix),
  "genes x",
  ncol(count_matrix),
  "cells"
)


# ==========================================================
# 8. Create raw Seurat object
# ==========================================================

log_message(
  "Creating raw Seurat object..."
)


seurat_raw <- CreateSeuratObject(
  counts = count_matrix,
  project = "GSE183904",
  min.cells = 0,
  min.features = 0
)


rm(
  count_matrix
)


gc()


# ----------------------------------------------------------
# Derive biological sample identity from prefixed barcode
# ----------------------------------------------------------

seurat_raw$sample <- sub(
  "_.*",
  "",
  colnames(seurat_raw)
)


# Keep orig.ident aligned with sample identity as well
seurat_raw$orig.ident <- seurat_raw$sample


checkpoint(
  length(unique(seurat_raw$sample)) > 1,
  paste0(
    "Only one sample identifier detected. ",
    "Check cell barcode naming."
  )
)


log_message(
  "Detected samples:",
  length(unique(seurat_raw$sample))
)


# ==========================================================
# 9. Save raw Seurat object
# ==========================================================

raw_object_file <- file.path(
  objects_dir,
  "seurat_raw_GSE183904.rds"
)


saveRDS(
  seurat_raw,
  raw_object_file
)


log_message(
  "Saved raw object:",
  raw_object_file
)


# ==========================================================
# 10. Calculate QC metrics
# ==========================================================

log_message(
  "Calculating QC metrics..."
)


seurat_raw[["percent.mt"]] <- PercentageFeatureSet(
  seurat_raw,
  pattern = "^MT-"
)


seurat_raw[["percent.ribo"]] <- PercentageFeatureSet(
  seurat_raw,
  pattern = "^RPL|^RPS"
)


seurat_raw[["percent.hb"]] <- PercentageFeatureSet(
  seurat_raw,
  pattern = "^HB[AB]"
)


# ==========================================================
# 11. Save QC plots
# ==========================================================

qc_plot_file <- file.path(
  plots_dir,
  "QC_vlnplots_GSE183904.pdf"
)


pdf(
  qc_plot_file,
  width = 18,
  height = 10
)


print(
  
  VlnPlot(
    seurat_raw,
    features = c(
      "nFeature_RNA",
      "nCount_RNA",
      "percent.mt"
    ),
    group.by = "sample",
    ncol = 3,
    pt.size = 0.1
  )
  
)


print(
  
  VlnPlot(
    seurat_raw,
    features = "percent.hb",
    group.by = "sample",
    pt.size = 0.1
  ) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed"
    )
  
)


print(
  
  VlnPlot(
    seurat_raw,
    features = "percent.ribo",
    group.by = "sample",
    pt.size = 0.1
  ) +
    geom_hline(
      yintercept = 45,
      linetype = "dashed"
    )
  
)


print(
  
  VlnPlot(
    seurat_raw,
    features = "percent.mt",
    group.by = "sample",
    pt.size = 0.1
  ) +
    geom_hline(
      yintercept = 30,
      linetype = "dashed"
    )
  
)


dev.off()


log_message(
  "Saved QC plots:",
  qc_plot_file
)


# ==========================================================
# 12. Cell counts before QC
# ==========================================================

before_qc <- as.data.frame(
  table(
    seurat_raw$sample
  )
)


colnames(before_qc) <- c(
  "sample",
  "cell_count"
)


before_qc_file <- file.path(
  tables_dir,
  "cell_numbers_before_QC_GSE183904.csv"
)


write.csv(
  before_qc,
  before_qc_file,
  row.names = FALSE
)


# ==========================================================
# 13. Apply QC filters
# ==========================================================

log_message(
  "Applying quality-control filters..."
)


seurat_qc <- subset(
  seurat_raw,
  subset =
    nFeature_RNA > 300 &
    nCount_RNA > 1000 &
    nCount_RNA < 50000 &
    percent.mt < 30 &
    percent.ribo < 45 &
    percent.hb < 1
)


checkpoint(
  ncol(seurat_qc) > 0,
  "QC filtering removed all cells."
)


log_message(
  "Cells before QC:",
  ncol(seurat_raw)
)


log_message(
  "Cells after QC:",
  ncol(seurat_qc)
)


log_message(
  "Cells removed by QC:",
  ncol(seurat_raw) -
    ncol(seurat_qc)
)


# ==========================================================
# 14. Cell counts after QC
# ==========================================================

after_qc <- as.data.frame(
  table(
    seurat_qc$sample
  )
)


colnames(after_qc) <- c(
  "sample",
  "cell_count"
)


after_qc_file <- file.path(
  tables_dir,
  "cell_numbers_after_QC_GSE183904.csv"
)


write.csv(
  after_qc,
  after_qc_file,
  row.names = FALSE
)


# ==========================================================
# 15. Save QC-filtered object
# ==========================================================

qc_object_file <- file.path(
  objects_dir,
  "seurat_QC_GSE183904.rds"
)


saveRDS(
  seurat_qc,
  qc_object_file
)


log_message(
  "Saved QC-filtered object:",
  qc_object_file
)


rm(
  seurat_raw
)


gc()


############################################################
#
#                DOUBLET DETECTION
#
############################################################


# ==========================================================
# 16. Split QC object by biological sample
# ==========================================================

log_message(
  "Splitting QC-filtered dataset by sample..."
)


sample_list <- SplitObject(
  seurat_qc,
  split.by = "sample"
)


checkpoint(
  length(sample_list) > 0,
  "SplitObject returned no samples."
)


log_message(
  "Samples prepared for DoubletFinder:",
  length(sample_list)
)


# ==========================================================
# 17. DoubletFinder function
# ==========================================================

run_doubletfinder <- function(
    sample_obj,
    sample_name,
    expected_doublet_rate = doublet_rate,
    pN = doublet_pN,
    cluster_resolution = doublet_cluster_resolution
) {
  
  set.seed(123)
  
  
  log_message(
    paste0(
      "[DoubletFinder] Starting sample: ",
      sample_name
    )
  )
  
  
  initial_cells <- ncol(
    sample_obj
  )
  
  
  checkpoint(
    initial_cells >= minimum_cells_for_doubletfinder,
    paste0(
      "Sample ",
      sample_name,
      " contains only ",
      initial_cells,
      " cells after QC. ",
      "DoubletFinder was not run."
    )
  )
  
  
  # --------------------------------------------------------
  # 17A. Normalize sample independently
  # --------------------------------------------------------
  
  sample_obj <- NormalizeData(
    sample_obj,
    verbose = FALSE
  )
  
  
  sample_obj <- FindVariableFeatures(
    sample_obj,
    selection.method = "vst",
    nfeatures = n_variable_features,
    verbose = FALSE
  )
  
  
  sample_obj <- ScaleData(
    sample_obj,
    verbose = FALSE
  )
  
  
  # --------------------------------------------------------
  # 17B. PCA
  # --------------------------------------------------------
  
  max_pcs <- min(
    30,
    initial_cells - 1,
    length(
      VariableFeatures(sample_obj)
    )
  )
  
  
  checkpoint(
    max_pcs >= 5,
    paste0(
      "Insufficient PCA dimensions for sample ",
      sample_name
    )
  )
  
  
  sample_obj <- RunPCA(
    sample_obj,
    npcs = max_pcs,
    seed.use = 123,
    verbose = FALSE
  )
  
  
  num_pcs <- choose_pcs(
    sample_obj,
    threshold = 0.90
  )
  
  
  num_pcs <- min(
    num_pcs,
    max_pcs
  )
  
  
  num_pcs <- max(
    num_pcs,
    2
  )
  
  
  pcs_to_use <- seq_len(
    num_pcs
  )
  
  
  log_message(
    paste0(
      "[DoubletFinder] ",
      sample_name,
      ": using ",
      num_pcs,
      " PCs"
    )
  )
  
  
  # --------------------------------------------------------
  # 17C. Preliminary clustering
  #
  # Used to estimate the proportion of homotypic doublets.
  # These clusters are NOT the final study clusters.
  # --------------------------------------------------------
  
  sample_obj <- FindNeighbors(
    sample_obj,
    reduction = "pca",
    dims = pcs_to_use,
    verbose = FALSE
  )
  
  
  sample_obj <- FindClusters(
    sample_obj,
    resolution = cluster_resolution,
    random.seed = 123,
    verbose = FALSE
  )
  
  
  preliminary_annotations <- as.character(
    Idents(sample_obj)
  )
  
  
  homotypic_proportion <- modelHomotypic(
    preliminary_annotations
  )
  
  
  # --------------------------------------------------------
  # 17D. Determine optimal pK
  # --------------------------------------------------------
  
  log_message(
    paste0(
      "[DoubletFinder] ",
      sample_name,
      ": starting pK sweep"
    )
  )
  
  
  sweep_results <- paramSweep(
    sample_obj,
    PCs = pcs_to_use,
    sct = FALSE
  )
  
  
  sweep_stats <- summarizeSweep(
    sweep_results,
    GT = FALSE
  )
  
  
  bcmvn <- find.pK(
    sweep_stats
  )
  
  
  checkpoint(
    nrow(bcmvn) > 0,
    paste0(
      "No pK results produced for sample ",
      sample_name
    )
  )
  
  
  bcmvn$pK_numeric <- as.numeric(
    as.character(
      bcmvn$pK
    )
  )
  
  
  valid_pk <- which(
    is.finite(bcmvn$pK_numeric) &
      is.finite(bcmvn$BCmetric)
  )
  
  
  checkpoint(
    length(valid_pk) > 0,
    paste0(
      "No valid pK candidates for sample ",
      sample_name
    )
  )
  
  
  best_index <- valid_pk[
    which.max(
      bcmvn$BCmetric[
        valid_pk
      ]
    )
  ]
  
  
  optimal_pK <- bcmvn$pK_numeric[
    best_index
  ]
  
  
  checkpoint(
    is.finite(optimal_pK),
    paste0(
      "Failed to determine optimal pK for ",
      sample_name
    )
  )
  
  
  # Save pK sweep
  pk_file <- file.path(
    doublet_tables_dir,
    paste0(
      "DoubletFinder_pK_",
      sample_name,
      ".csv"
    )
  )
  
  
  write.csv(
    bcmvn,
    pk_file,
    row.names = FALSE
  )
  
  
  log_message(
    paste0(
      "[DoubletFinder] ",
      sample_name,
      ": optimal pK = ",
      optimal_pK
    )
  )
  
  
  # --------------------------------------------------------
  # 17E. Estimate expected doublets
  # --------------------------------------------------------
  
  nExp_poisson <- round(
    expected_doublet_rate *
      initial_cells
  )
  
  
  nExp_poisson <- max(
    1,
    nExp_poisson
  )
  
  
  nExp_adjusted <- round(
    nExp_poisson *
      (
        1 -
          homotypic_proportion
      )
  )
  
  
  nExp_adjusted <- max(
    1,
    nExp_adjusted
  )
  
  
  log_message(
    paste0(
      "[DoubletFinder] ",
      sample_name,
      ": expected doublets = ",
      nExp_poisson,
      "; homotypic-adjusted = ",
      nExp_adjusted
    )
  )
  
  
  # --------------------------------------------------------
  # 17F. First DoubletFinder run
  #
  # Generates the pANN score using unadjusted Poisson
  # estimate.
  # --------------------------------------------------------
  
  metadata_before_df <- colnames(
    sample_obj@meta.data
  )
  
  
  sample_obj <- doubletFinder(
    sample_obj,
    PCs = pcs_to_use,
    pN = pN,
    pK = optimal_pK,
    nExp = nExp_poisson,
    reuse.pANN = NULL,
    sct = FALSE
  )
  
  
  new_metadata <- setdiff(
    colnames(sample_obj@meta.data),
    metadata_before_df
  )
  
  
  new_pann_columns <- grep(
    "^pANN",
    new_metadata,
    value = TRUE
  )
  
  
  checkpoint(
    length(new_pann_columns) > 0,
    paste0(
      "DoubletFinder failed to create pANN for ",
      sample_name
    )
  )
  
  
  pANN_column <- tail(
    new_pann_columns,
    1
  )
  
  
  # --------------------------------------------------------
  # 17G. Second DoubletFinder classification
  #
  # Reuses pANN and applies the homotypic-adjusted estimate.
  # --------------------------------------------------------
  
  classification_before <- grep(
    "^DF.classifications",
    colnames(sample_obj@meta.data),
    value = TRUE
  )
  
  
  sample_obj <- doubletFinder(
    sample_obj,
    PCs = pcs_to_use,
    pN = pN,
    pK = optimal_pK,
    nExp = nExp_adjusted,
    reuse.pANN = pANN_column,
    sct = FALSE
  )
  
  
  classification_after <- grep(
    "^DF.classifications",
    colnames(sample_obj@meta.data),
    value = TRUE
  )
  
  
  final_classification_candidates <- setdiff(
    classification_after,
    classification_before
  )
  
  
  if (
    length(final_classification_candidates) == 0
  ) {
    
    final_classification_column <- tail(
      classification_after,
      1
    )
    
  } else {
    
    final_classification_column <- tail(
      final_classification_candidates,
      1
    )
    
  }
  
  
  checkpoint(
    length(final_classification_column) == 1,
    paste0(
      "Could not identify final DoubletFinder ",
      "classification column for ",
      sample_name
    )
  )
  
  
  # --------------------------------------------------------
  # 17H. Save standardized doublet metadata
  # --------------------------------------------------------
  
  sample_obj$doublet_status <- sample_obj@meta.data[
    [
      final_classification_column
    ]
  ]
  
  
  sample_obj$doublet_pANN <- sample_obj@meta.data[
    [
      pANN_column
    ]
  ]
  
  
  doublet_table <- table(
    sample_obj$doublet_status
  )
  
  
  n_detected_doublets <- sum(
    sample_obj$doublet_status ==
      "Doublet",
    na.rm = TRUE
  )
  
  
  n_detected_singlets <- sum(
    sample_obj$doublet_status ==
      "Singlet",
    na.rm = TRUE
  )
  
  
  log_message(
    paste0(
      "[DoubletFinder] ",
      sample_name,
      ": predicted doublets = ",
      n_detected_doublets
    )
  )
  
  
  # --------------------------------------------------------
  # 17I. Remove predicted doublets
  # --------------------------------------------------------
  
  singlet_obj <- subset(
    sample_obj,
    subset =
      doublet_status ==
      "Singlet"
  )
  
  
  checkpoint(
    ncol(singlet_obj) > 0,
    paste0(
      "No singlets remained for sample ",
      sample_name
    )
  )
  
  
  # --------------------------------------------------------
  # 17J. Create summary row
  # --------------------------------------------------------
  
  summary_row <- data.frame(
    
    sample = sample_name,
    
    cells_before_doubletfinder =
      initial_cells,
    
    expected_doublet_rate =
      expected_doublet_rate,
    
    expected_doublets_poisson =
      nExp_poisson,
    
    homotypic_proportion =
      homotypic_proportion,
    
    expected_doublets_adjusted =
      nExp_adjusted,
    
    optimal_pK =
      optimal_pK,
    
    predicted_doublets =
      n_detected_doublets,
    
    retained_singlets =
      n_detected_singlets,
    
    cells_after_doubletfinder =
      ncol(singlet_obj),
    
    stringsAsFactors = FALSE
    
  )
  
  
  return(
    list(
      singlets = singlet_obj,
      summary = summary_row
    )
  )
  
}


# ==========================================================
# 18. Run DoubletFinder sample-by-sample
# ==========================================================

log_message(
  "Starting sample-level DoubletFinder analysis..."
)


doublet_results <- vector(
  mode = "list",
  length = length(sample_list)
)


names(doublet_results) <- names(
  sample_list
)


for (sample_name in names(sample_list)) {
  
  log_message(
    paste0(
      "Processing ",
      sample_name
    )
  )
  
  
  doublet_results[[sample_name]] <-
    run_doubletfinder(
      
      sample_obj =
        sample_list[[sample_name]],
      
      sample_name =
        sample_name,
      
      expected_doublet_rate =
        doublet_rate,
      
      pN =
        doublet_pN,
      
      cluster_resolution =
        doublet_cluster_resolution
      
    )
  
}


# ==========================================================
# 19. Extract singlet objects
# ==========================================================

singlet_list <- lapply(
  doublet_results,
  function(x) {
    x$singlets
  }
)


names(singlet_list) <- names(
  doublet_results
)


checkpoint(
  length(singlet_list) > 0,
  "No singlet objects were generated."
)


# ==========================================================
# 20. DoubletFinder summary table
# ==========================================================

doublet_summary <- do.call(
  rbind,
  lapply(
    doublet_results,
    function(x) {
      x$summary
    }
  )
)


rownames(doublet_summary) <- NULL


doublet_summary_file <- file.path(
  tables_dir,
  "DoubletFinder_summary_GSE183904.csv"
)


write.csv(
  doublet_summary,
  doublet_summary_file,
  row.names = FALSE
)


log_message(
  "DoubletFinder summary saved:",
  doublet_summary_file
)


# ==========================================================
# 21. Remove temporary objects before merge
# ==========================================================

rm(
  sample_list,
  doublet_results,
  seurat_qc
)


gc()


# ==========================================================
# 22. Merge singlet-only samples
# ==========================================================

log_message(
  "Merging singlet-only sample objects..."
)


sample_names <- names(
  singlet_list
)


checkpoint(
  length(sample_names) > 0,
  "No sample names found in singlet list."
)


if (length(sample_names) == 1) {
  
  seurat_singlets <- singlet_list[
    [1]
  ]
  
} else {
  
  seurat_singlets <- merge(
    
    x =
      singlet_list[[1]],
    
    y =
      singlet_list[-1],
    
    project =
      "GSE183904_singlets",
    
    merge.data =
      FALSE
    
  )
  
}


rm(
  singlet_list
)


gc()


# ==========================================================
# 23. Join Seurat v5 assay layers when required
# ==========================================================

if (
  packageVersion("SeuratObject") >=
  package_version("5.0.0")
) {
  
  log_message(
    "Joining Seurat v5 assay layers..."
  )
  
  
  seurat_singlets <- JoinLayers(
    seurat_singlets
  )
  
}


# ==========================================================
# 24. Reset sample/orig.ident consistency
# ==========================================================

checkpoint(
  "sample" %in%
    colnames(
      seurat_singlets@meta.data
    ),
  "Sample metadata was lost during merge."
)


seurat_singlets$orig.ident <-
  seurat_singlets$sample


# ==========================================================
# 25. Save singlet-only object
# ==========================================================

singlet_object_file <- file.path(
  objects_dir,
  "seurat_singlets_GSE183904.rds"
)


saveRDS(
  seurat_singlets,
  singlet_object_file
)


log_message(
  "Saved singlet-only object:",
  singlet_object_file
)


# ==========================================================
# 26. Global normalization
# ==========================================================

log_message(
  "Normalizing merged singlet dataset..."
)


seurat_singlets <- NormalizeData(
  seurat_singlets
)


# ==========================================================
# 27. Variable-feature selection
# ==========================================================

log_message(
  "Identifying variable features..."
)


seurat_singlets <- FindVariableFeatures(
  seurat_singlets,
  selection.method = "vst",
  nfeatures = n_variable_features
)


# ==========================================================
# 28. Scale expression data
# ==========================================================

log_message(
  "Scaling expression data..."
)


seurat_singlets <- ScaleData(
  seurat_singlets
)


# ==========================================================
# 29. PCA
# ==========================================================

log_message(
  "Running PCA..."
)


seurat_singlets <- RunPCA(
  seurat_singlets,
  npcs = 50,
  seed.use = 123
)


# ==========================================================
# 30. Save elbow plot
# ==========================================================

elbow_file <- file.path(
  plots_dir,
  "ElbowPlot_GSE183904.pdf"
)


pdf(
  elbow_file,
  width = 6,
  height = 5
)


print(
  ElbowPlot(
    seurat_singlets,
    ndims = 50
  )
)


dev.off()


# ==========================================================
# 31. Select PCs
# ==========================================================

num_pcs <- choose_pcs(
  seurat_singlets,
  threshold = 0.90
)


checkpoint(
  !is.na(num_pcs) &&
    num_pcs > 0,
  "Unable to determine number of PCs."
)


pcs_to_use <- seq_len(
  num_pcs
)


log_message(
  "PCs retained for Harmony:",
  num_pcs
)


# ==========================================================
# 32. Pre-Harmony UMAP
# ==========================================================

log_message(
  "Running UMAP before Harmony integration..."
)


seurat_singlets <- RunUMAP(
  seurat_singlets,
  reduction = "pca",
  dims = pcs_to_use,
  seed.use = 123,
  reduction.name = "umap.preharmony",
  reduction.key = "UMAPPRE_"
)


pre_harmony_plot <- DimPlot(
  seurat_singlets,
  reduction = "umap.preharmony",
  group.by = "sample"
) +
  ggtitle(
    "GSE183904 — Before Harmony"
  )


ggsave(
  filename = file.path(
    plots_dir,
    "UMAP_preHarmony_GSE183904.pdf"
  ),
  plot = pre_harmony_plot,
  width = 10,
  height = 8
)


# ==========================================================
# 33. Harmony batch correction
# ==========================================================

log_message(
  "Running Harmony integration..."
)


seurat_integrated <- RunHarmony(
  
  object =
    seurat_singlets,
  
  group.by.vars =
    "sample",
  
  reduction.use =
    "pca",
  
  dims.use =
    pcs_to_use,
  
  reduction.save =
    "harmony",
  
  plot_convergence =
    FALSE
  
)


checkpoint(
  "harmony" %in%
    names(
      seurat_integrated@reductions
    ),
  "Harmony reduction was not generated."
)


# ==========================================================
# 34. Harmony-based UMAP
# ==========================================================

log_message(
  "Running Harmony-based UMAP..."
)


seurat_integrated <- RunUMAP(
  
  seurat_integrated,
  
  reduction =
    "harmony",
  
  dims =
    pcs_to_use,
  
  seed.use =
    123
  
)


# ==========================================================
# 35. Nearest-neighbor graph
# ==========================================================

log_message(
  "Constructing nearest-neighbor graph..."
)


seurat_integrated <- FindNeighbors(
  
  seurat_integrated,
  
  reduction =
    "harmony",
  
  dims =
    pcs_to_use
  
)


# ==========================================================
# 36. Graph-based clustering
# ==========================================================

log_message(
  "Running global clustering..."
)


seurat_integrated <- FindClusters(
  
  seurat_integrated,
  
  resolution =
    global_cluster_resolution,
  
  random.seed =
    123
  
)


# ==========================================================
# 37. Sanity checks
# ==========================================================

checkpoint(
  "umap" %in%
    names(
      seurat_integrated@reductions
    ),
  "Final UMAP reduction missing."
)


checkpoint(
  "seurat_clusters" %in%
    colnames(
      seurat_integrated@meta.data
    ),
  "seurat_clusters missing from metadata."
)


log_message(
  "Final number of cells:",
  ncol(seurat_integrated)
)


log_message(
  "Final number of clusters:",
  length(
    unique(
      seurat_integrated$seurat_clusters
    )
  )
)


# ==========================================================
# 38. Final cluster UMAP
# ==========================================================

cluster_plot <- DimPlot(
  
  seurat_integrated,
  
  reduction =
    "umap",
  
  group.by =
    "seurat_clusters",
  
  label =
    TRUE,
  
  repel =
    TRUE,
  
  raster =
    FALSE
  
) +
  ggtitle(
    "GSE183904 — Seurat Clusters"
  )


ggsave(
  
  filename = file.path(
    plots_dir,
    "UMAP_clusters_GSE183904.pdf"
  ),
  
  plot =
    cluster_plot,
  
  width =
    10,
  
  height =
    8
  
)


# ==========================================================
# 39. Sample UMAP
# ==========================================================

sample_plot <- DimPlot(
  
  seurat_integrated,
  
  reduction =
    "umap",
  
  group.by =
    "sample"
  
) +
  ggtitle(
    "GSE183904 — Samples"
  )


ggsave(
  
  filename = file.path(
    plots_dir,
    "UMAP_samples_GSE183904.pdf"
  ),
  
  plot =
    sample_plot,
  
  width =
    10,
  
  height =
    8
  
)


# ==========================================================
# 40. Cluster-size summary
# ==========================================================

cluster_counts <- as.data.frame(
  table(
    seurat_integrated$seurat_clusters
  )
)


colnames(
  cluster_counts
) <- c(
  "seurat_cluster",
  "cell_count"
)


write.csv(
  
  cluster_counts,
  
  file.path(
    tables_dir,
    "cluster_cell_counts_GSE183904.csv"
  ),
  
  row.names =
    FALSE
  
)


# ==========================================================
# 41. Final per-sample cell counts
# ==========================================================

final_sample_counts <- as.data.frame(
  table(
    seurat_integrated$sample
  )
)


colnames(
  final_sample_counts
) <- c(
  "sample",
  "final_cell_count"
)


write.csv(
  
  final_sample_counts,
  
  file.path(
    tables_dir,
    "final_cell_numbers_GSE183904.csv"
  ),
  
  row.names =
    FALSE
  
)


# ==========================================================
# 42. Save final integrated object
# ==========================================================

integrated_object_file <- file.path(
  objects_dir,
  "seurat_integrated_GSE183904.rds"
)


saveRDS(
  
  seurat_integrated,
  
  integrated_object_file
  
)


log_message(
  "Final integrated object saved:",
  integrated_object_file
)


# ==========================================================
# 43. Save session information
# ==========================================================

session_info_file <- file.path(
  results_dir,
  "sessionInfo_01_preprocessing.txt"
)


capture.output(
  
  sessionInfo(),
  
  file =
    session_info_file
  
)


# ==========================================================
# 44. Final pipeline summary
# ==========================================================

cat(
  "\n",
  "=====================================================\n",
  " GSE183904 PREPROCESSING COMPLETED SUCCESSFULLY\n",
  "=====================================================\n",
  "\n",
  "QC-filtered object:\n",
  qc_object_file,
  "\n\n",
  "Doublet-filtered object:\n",
  singlet_object_file,
  "\n\n",
  "Final integrated object:\n",
  integrated_object_file,
  "\n\n",
  "Figures:\n",
  plots_dir,
  "\n\n",
  "Tables:\n",
  tables_dir,
  "\n",
  "=====================================================\n"
)


log_message(
  "01_download_qc_integration.R completed."
)