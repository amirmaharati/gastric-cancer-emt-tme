############################################################
# Epithelial Re-Clustering and EMT State Classification
#
# Project:
#   Gastric Cancer EMT–TME Analysis
#
# Dataset:
#   GSE183904
#
# Input:
#   objects/seurat_annotated_GSE183904.rds
#
# Outputs:
#   objects/seurat_epithelial_EMT_GSE183904.rds
#   objects/seurat_final_GSE183904.rds
#   objects/EMT_GMM_model_GSE183904.rds
#
# Workflow:
#   01. Load annotated global Seurat object
#   02. Subset epithelial cells
#   03. Re-normalize epithelial cells
#   04. Select 2,000 variable genes
#   05. Scale expression
#   06. PCA
#   07. Harmony integration
#   08. UMAP and epithelial re-clustering (resolution = 2.0)
#   09. Calculate epithelial and mesenchymal module scores
#   10. Calculate mean normalized-expression scores
#   11. Calculate EMT C-score (MES - EPI)
#   12. Z-transform EMT scores
#   13. Fit 5-component Gaussian mixture model
#   14. Order GMM components by EMT C-score
#   15. Define EMT_rank1 ... EMT_rank5
#   16. Collapse ranks into:
#          EMT_low
#          EMT_intermediate
#          EMT_high
#   17. Transfer EMT annotations to global object
#
# NOTE:
#   Trajectory/pseudotime analysis is intentionally NOT
#   performed here. It belongs to:
#
#   05_emt_trajectory_monocle3.R
############################################################


# ==========================================================
# 0. Required packages
# ==========================================================

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "ggplot2",
  "harmony",
  "mclust"
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
  library(Matrix)
  library(ggplot2)
  library(harmony)
  library(mclust)
  
})


# ==========================================================
# 1. Reproducibility and analysis parameters
# ==========================================================

set.seed(123)


n_variable_features <- 2000


epithelial_cluster_resolution <- 2.0


pca_variance_threshold <- 0.90


max_pcs <- 50


# ----------------------------------------------------------
# EMT GMM
# ----------------------------------------------------------

gmm_components <- 5


# Original analysis used a 2-dimensional representation
# of mesenchymal and epithelial scores for GMM fitting.

gmm_mode <- "2D"


# Minimum fraction of each supplied EMT gene set that must
# be found in the expression matrix.

min_fraction_genes_present <- 0.40


# ==========================================================
# 2. EMT gene sets
# ==========================================================

emt_mesenchymal_genes <- c(
  
  "VIM",
  "FN1",
  "COL1A1",
  "COL1A2",
  "COL3A1",
  "ITGA5",
  "ITGB1",
  "SNAI1",
  "SNAI2",
  "ZEB1",
  "ZEB2",
  "TWIST1",
  "TAGLN",
  "ACTA2",
  "MMP2",
  "MMP9",
  "CDH2",
  "S100A4"
  
)


emt_epithelial_genes <- c(
  
  "CDH1",
  "EPCAM",
  "KRT8",
  "KRT18",
  "KRT19",
  "MUC1",
  "CLDN3",
  "CLDN4",
  "CLDN7",
  "OCLN",
  "TJP1"
  
)


# ==========================================================
# 3. Project directories
# ==========================================================

base_dir <- "."


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


# ==========================================================
# 4. Helper functions
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
# Determine number of PCs explaining the requested
# proportion of variance.
# ----------------------------------------------------------

choose_pcs <- function(
    seurat_object,
    threshold = 0.90
) {
  
  pca_sd <- Stdev(
    seurat_object,
    reduction = "pca"
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


# ----------------------------------------------------------
# Retrieve normalized RNA expression.
#
# Compatible with SeuratObject v4/v5 workflows.
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


# ==========================================================
# 5. Load annotated global Seurat object
# ==========================================================

input_file <- file.path(
  objects_dir,
  "seurat_annotated_GSE183904.rds"
)


checkpoint(
  file.exists(input_file),
  paste0(
    "Annotated Seurat object not found:\n",
    input_file,
    "\n\nRun scripts 01, 02, and 03 first."
  )
)


log_message(
  "Loading annotated GSE183904 object..."
)


seurat_global <- readRDS(
  input_file
)


checkpoint(
  inherits(seurat_global, "Seurat"),
  "Input file is not a valid Seurat object."
)


checkpoint(
  "cell_type" %in%
    colnames(
      seurat_global@meta.data
    ),
  "cell_type metadata is missing."
)


checkpoint(
  "sample" %in%
    colnames(
      seurat_global@meta.data
    ),
  "sample metadata is missing."
)


log_message(
  "Total cells:",
  ncol(seurat_global)
)


# ==========================================================
# 6. Subset epithelial cells
# ==========================================================

log_message(
  "Subsetting epithelial cells..."
)


seurat_epithelial <- subset(
  seurat_global,
  subset = cell_type == "Epithelial"
)


checkpoint(
  ncol(seurat_epithelial) > 0,
  "No cells annotated as Epithelial were found."
)


log_message(
  "Epithelial cells:",
  ncol(seurat_epithelial)
)


# ==========================================================
# 7. Re-normalize epithelial subset
# ==========================================================

DefaultAssay(
  seurat_epithelial
) <- "RNA"


log_message(
  "Normalizing epithelial-cell expression..."
)


seurat_epithelial <- NormalizeData(
  seurat_epithelial,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)


# ==========================================================
# 8. Select variable features
# ==========================================================

log_message(
  "Selecting variable features..."
)


seurat_epithelial <- FindVariableFeatures(
  
  seurat_epithelial,
  
  selection.method =
    "vst",
  
  nfeatures =
    n_variable_features
  
)


# ==========================================================
# 9. Scale expression
# ==========================================================

log_message(
  "Scaling epithelial-cell expression..."
)


seurat_epithelial <- ScaleData(
  seurat_epithelial
)


# ==========================================================
# 10. PCA
# ==========================================================

log_message(
  "Running epithelial-cell PCA..."
)


seurat_epithelial <- RunPCA(
  
  seurat_epithelial,
  
  npcs =
    max_pcs,
  
  seed.use =
    123
  
)


# ==========================================================
# 11. Epithelial PCA elbow plot
# ==========================================================

elbow_plot <- ElbowPlot(
  seurat_epithelial,
  ndims = max_pcs
)


ggsave(
  
  filename = file.path(
    plots_dir,
    "EMT_01_epithelial_elbow_plot.pdf"
  ),
  
  plot =
    elbow_plot,
  
  width =
    6,
  
  height =
    5
  
)


# ==========================================================
# 12. Select PCs explaining >=90% variance
# ==========================================================

num_pcs <- choose_pcs(
  
  seurat_epithelial,
  
  threshold =
    pca_variance_threshold
  
)


num_pcs <- min(
  num_pcs,
  max_pcs
)


checkpoint(
  num_pcs >= 2,
  "Fewer than two usable PCs were identified."
)


pcs_to_use <- seq_len(
  num_pcs
)


log_message(
  "PCs retained:",
  num_pcs
)


# ==========================================================
# 13. Pre-Harmony epithelial UMAP
# ==========================================================

log_message(
  "Running pre-Harmony epithelial UMAP..."
)


seurat_epithelial <- RunUMAP(
  
  seurat_epithelial,
  
  reduction =
    "pca",
  
  dims =
    pcs_to_use,
  
  seed.use =
    123,
  
  reduction.name =
    "umap.preharmony",
  
  reduction.key =
    "UMAPPRE_"
  
)


p_pre_harmony <- DimPlot(
  
  seurat_epithelial,
  
  reduction =
    "umap.preharmony",
  
  group.by =
    "sample"
  
) +
  ggtitle(
    "Epithelial Cells — Before Harmony"
  )


ggsave(
  
  filename = file.path(
    plots_dir,
    "EMT_02_epithelial_preHarmony_UMAP.pdf"
  ),
  
  plot =
    p_pre_harmony,
  
  width =
    10,
  
  height =
    8
  
)


# ==========================================================
# 14. Harmony integration
# ==========================================================

log_message(
  "Running Harmony on epithelial cells..."
)


seurat_epithelial <- RunHarmony(
  
  object =
    seurat_epithelial,
  
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
      seurat_epithelial@reductions
    ),
  "Harmony reduction was not created."
)


# ==========================================================
# 15. Harmony-based UMAP
# ==========================================================

seurat_epithelial <- RunUMAP(
  
  seurat_epithelial,
  
  reduction =
    "harmony",
  
  dims =
    pcs_to_use,
  
  seed.use =
    123
  
)


# ==========================================================
# 16. Epithelial nearest-neighbor graph
# ==========================================================

seurat_epithelial <- FindNeighbors(
  
  seurat_epithelial,
  
  reduction =
    "harmony",
  
  dims =
    pcs_to_use
  
)


# ==========================================================
# 17. Epithelial clustering
#
# Final selected resolution = 2.0
# ==========================================================

log_message(
  "Clustering epithelial cells at resolution 2.0..."
)


seurat_epithelial <- FindClusters(
  
  seurat_epithelial,
  
  resolution =
    epithelial_cluster_resolution,
  
  random.seed =
    123
  
)


checkpoint(
  "seurat_clusters" %in%
    colnames(
      seurat_epithelial@meta.data
    ),
  "Epithelial clustering failed."
)


log_message(
  "Number of epithelial clusters:",
  length(
    unique(
      seurat_epithelial$seurat_clusters
    )
  )
)


# ==========================================================
# 18. Save epithelial cluster UMAP
# ==========================================================

p_epi_clusters <- DimPlot(
  
  seurat_epithelial,
  
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
    "Epithelial Subclusters — Resolution 2.0"
  )


ggsave(
  
  filename = file.path(
    plots_dir,
    "EMT_03_epithelial_clusters_resolution_2.pdf"
  ),
  
  plot =
    p_epi_clusters,
  
  width =
    10,
  
  height =
    8
  
)


############################################################
#
#                    EMT SCORING
#
############################################################


# ==========================================================
# 19. Retrieve normalized expression
# ==========================================================

expression_matrix <- get_normalized_expression(
  seurat_epithelial
)


gene_names <- rownames(
  expression_matrix
)


# ==========================================================
# 20. Identify available EMT genes
# ==========================================================

emt_mes <- intersect(
  emt_mesenchymal_genes,
  gene_names
)


emt_epi <- intersect(
  emt_epithelial_genes,
  gene_names
)


mes_fraction <- length(
  emt_mes
) / length(
  emt_mesenchymal_genes
)


epi_fraction <- length(
  emt_epi
) / length(
  emt_epithelial_genes
)


log_message(
  "Mesenchymal genes found:",
  length(emt_mes),
  "/",
  length(emt_mesenchymal_genes)
)


log_message(
  "Epithelial genes found:",
  length(emt_epi),
  "/",
  length(emt_epithelial_genes)
)


checkpoint(
  mes_fraction >=
    min_fraction_genes_present,
  "Too few mesenchymal EMT genes were detected."
)


checkpoint(
  epi_fraction >=
    min_fraction_genes_present,
  "Too few epithelial EMT genes were detected."
)


# ==========================================================
# 21. Save EMT gene-set availability
# ==========================================================

emt_gene_table <- rbind(
  
  data.frame(
    
    gene =
      emt_mesenchymal_genes,
    
    gene_set =
      "Mesenchymal",
    
    detected =
      emt_mesenchymal_genes %in%
      gene_names
    
  ),
  
  data.frame(
    
    gene =
      emt_epithelial_genes,
    
    gene_set =
      "Epithelial",
    
    detected =
      emt_epithelial_genes %in%
      gene_names
    
  )
  
)


write.csv(
  
  emt_gene_table,
  
  file.path(
    tables_dir,
    "EMT_gene_sets_GSE183904.csv"
  ),
  
  row.names =
    FALSE
  
)


# ==========================================================
# 22. Mesenchymal AddModuleScore
# ==========================================================

log_message(
  "Calculating mesenchymal module score..."
)


seurat_epithelial <- AddModuleScore(
  
  object =
    seurat_epithelial,
  
  features =
    list(
      emt_mes
    ),
  
  assay =
    "RNA",
  
  name =
    "EMT_MES_AMS",
  
  seed =
    123
  
)


checkpoint(
  "EMT_MES_AMS1" %in%
    colnames(
      seurat_epithelial@meta.data
    ),
  "Mesenchymal AddModuleScore was not generated."
)


seurat_epithelial$EMT_MES_AddModule <-
  seurat_epithelial$EMT_MES_AMS1


# ==========================================================
# 23. Epithelial AddModuleScore
# ==========================================================

log_message(
  "Calculating epithelial module score..."
)


seurat_epithelial <- AddModuleScore(
  
  object =
    seurat_epithelial,
  
  features =
    list(
      emt_epi
    ),
  
  assay =
    "RNA",
  
  name =
    "EMT_EPI_AMS",
  
  seed =
    123
  
)


checkpoint(
  "EMT_EPI_AMS1" %in%
    colnames(
      seurat_epithelial@meta.data
    ),
  "Epithelial AddModuleScore was not generated."
)


seurat_epithelial$EMT_EPI_AddModule <-
  seurat_epithelial$EMT_EPI_AMS1


# ==========================================================
# 24. Mean normalized-expression scores
# ==========================================================

log_message(
  "Calculating mean-expression EMT scores..."
)


seurat_epithelial$EMT_MES_mean <-
  
  Matrix::colMeans(
    
    expression_matrix[
      emt_mes,
      ,
      drop = FALSE
    ]
    
  )


seurat_epithelial$EMT_EPI_mean <-
  
  Matrix::colMeans(
    
    expression_matrix[
      emt_epi,
      ,
      drop = FALSE
    ]
    
  )


# ==========================================================
# 25. Define primary EMT MES and EPI scores
#
# AddModuleScore is used as the primary score.
# Mean expression is retained as an independent
# concordance check.
# ==========================================================

seurat_epithelial$EMT_MES <-
  seurat_epithelial$EMT_MES_AddModule


seurat_epithelial$EMT_EPI <-
  seurat_epithelial$EMT_EPI_AddModule


# ==========================================================
# 26. EMT composite C-score
#
# C-score = Mesenchymal score - Epithelial score
# ==========================================================

seurat_epithelial$EMT_Cscore <-
  
  seurat_epithelial$EMT_MES -
  seurat_epithelial$EMT_EPI


# ==========================================================
# 27. Z-transform EMT scores
# ==========================================================

seurat_epithelial$EMT_MES_z <-
  
  as.numeric(
    scale(
      seurat_epithelial$EMT_MES
    )
  )


seurat_epithelial$EMT_EPI_z <-
  
  as.numeric(
    scale(
      seurat_epithelial$EMT_EPI
    )
  )


seurat_epithelial$EMT_Cscore_z <-
  
  as.numeric(
    scale(
      seurat_epithelial$EMT_Cscore
    )
  )


checkpoint(
  all(
    is.finite(
      seurat_epithelial$EMT_Cscore_z
    )
  ),
  "Non-finite EMT C-scores were generated."
)


# ==========================================================
# 28. Assess concordance with mean-expression scores
# ==========================================================

score_concordance <- data.frame(
  
  score = c(
    "Mesenchymal",
    "Epithelial"
  ),
  
  spearman_rho = c(
    
    cor(
      seurat_epithelial$EMT_MES_AddModule,
      seurat_epithelial$EMT_MES_mean,
      method = "spearman",
      use = "complete.obs"
    ),
    
    cor(
      seurat_epithelial$EMT_EPI_AddModule,
      seurat_epithelial$EMT_EPI_mean,
      method = "spearman",
      use = "complete.obs"
    )
    
  )
  
)


write.csv(
  
  score_concordance,
  
  file.path(
    tables_dir,
    "EMT_score_concordance_GSE183904.csv"
  ),
  
  row.names =
    FALSE
  
)


print(
  score_concordance
)


############################################################
#
#               GAUSSIAN MIXTURE MODEL
#
############################################################


# ==========================================================
# 29. Construct GMM input
# ==========================================================

if (gmm_mode == "2D") {
  
  gmm_input <- cbind(
    
    MES =
      seurat_epithelial$EMT_MES_z,
    
    EPI =
      seurat_epithelial$EMT_EPI_z
    
  )
  
} else if (gmm_mode == "1D") {
  
  gmm_input <-
    seurat_epithelial$EMT_Cscore_z
  
} else {
  
  stop(
    "gmm_mode must be '1D' or '2D'.",
    call. = FALSE
  )
  
}


checkpoint(
  all(
    is.finite(
      gmm_input
    )
  ),
  "GMM input contains non-finite values."
)


# ==========================================================
# 30. Fit exactly five GMM components
# ==========================================================

log_message(
  "Fitting 5-component EMT Gaussian mixture model..."
)


set.seed(123)


gmm_fit <- mclust::Mclust(
  
  data =
    gmm_input,
  
  G =
    gmm_components
  
)


checkpoint(
  !is.null(
    gmm_fit$classification
  ),
  "GMM classification failed."
)


checkpoint(
  gmm_fit$G ==
    gmm_components,
  paste0(
    "Expected ",
    gmm_components,
    " GMM components but model returned ",
    gmm_fit$G,
    "."
  )
)


log_message(
  "GMM model:",
  gmm_fit$modelName
)


# ==========================================================
# 31. Extract component assignments
# ==========================================================

gmm_component <-
  as.integer(
    gmm_fit$classification
  )


seurat_epithelial$EMT_GMM_component <-
  gmm_component


# ==========================================================
# 32. Order GMM components by EMT C-score
#
# Lowest mean C-score = EMT_rank1
# Highest mean C-score = EMT_rank5
# ==========================================================

component_mean_cscore <- tapply(
  
  seurat_epithelial$EMT_Cscore_z,
  
  gmm_component,
  
  mean,
  
  na.rm = TRUE
  
)


ordered_components <- as.integer(
  
  names(
    
    sort(
      component_mean_cscore,
      decreasing = FALSE
    )
    
  )
  
)


checkpoint(
  length(
    ordered_components
  ) == 5,
  "Five ordered GMM components were not identified."
)


rank_names <- paste0(
  "EMT_rank",
  seq_len(
    gmm_components
  )
)


component_to_rank <- setNames(
  
  rank_names,
  
  as.character(
    ordered_components
  )
  
)


emt_rank <- component_to_rank[
  as.character(
    gmm_component
  )
]


seurat_epithelial$EMT_GMM_rank <- factor(
  
  unname(
    emt_rank
  ),
  
  levels =
    rank_names
  
)


# ==========================================================
# 33. Collapse five ranks into three biological EMT groups
#
# rank1 + rank2 = EMT_low
# rank3         = EMT_intermediate
# rank4 + rank5 = EMT_high
# ==========================================================

rank_to_group <- c(
  
  EMT_rank1 =
    "EMT_low",
  
  EMT_rank2 =
    "EMT_low",
  
  EMT_rank3 =
    "EMT_intermediate",
  
  EMT_rank4 =
    "EMT_high",
  
  EMT_rank5 =
    "EMT_high"
  
)


emt_group <- rank_to_group[
  as.character(
    seurat_epithelial$EMT_GMM_rank
  )
]


seurat_epithelial$EMT_GMM_cluster <- factor(
  
  unname(
    emt_group
  ),
  
  levels = c(
    "EMT_low",
    "EMT_intermediate",
    "EMT_high"
  )
  
)


checkpoint(
  all(
    !is.na(
      seurat_epithelial$EMT_GMM_cluster
    )
  ),
  "Some epithelial cells were not assigned to an EMT group."
)


# ==========================================================
# 34. Save GMM posterior probabilities
# ==========================================================

gmm_posterior <- gmm_fit$z


checkpoint(
  nrow(gmm_posterior) ==
    ncol(seurat_epithelial),
  "GMM posterior probability dimensions do not match cells."
)


for (
  component_id
  in seq_len(
    gmm_components
  )
) {
  
  column_name <- paste0(
    "EMT_GMM_probability_component_",
    component_id
  )
  
  
  seurat_epithelial[
    [column_name]
  ] <- gmm_posterior[
    ,
    component_id
  ]
  
}


seurat_epithelial$EMT_GMM_confidence <-
  
  apply(
    gmm_posterior,
    1,
    max
  )


# ==========================================================
# 35. Create GMM component summary
# ==========================================================

component_counts <- table(
  seurat_epithelial$EMT_GMM_component
)


component_summary <- data.frame(
  
  GMM_component =
    ordered_components,
  
  EMT_rank =
    component_to_rank[
      as.character(
        ordered_components
      )
    ],
  
  EMT_group =
    rank_to_group[
      component_to_rank[
        as.character(
          ordered_components
        )
      ]
    ],
  
  mean_Cscore_z =
    as.numeric(
      component_mean_cscore[
        as.character(
          ordered_components
        )
      ]
    ),
  
  cell_count =
    as.numeric(
      component_counts[
        as.character(
          ordered_components
        )
      ]
    ),
  
  stringsAsFactors =
    FALSE
  
)


write.csv(
  
  component_summary,
  
  file.path(
    tables_dir,
    "EMT_GMM_component_summary_GSE183904.csv"
  ),
  
  row.names =
    FALSE
  
)


cat(
  "\nGMM component ordering:\n"
)


print(
  component_summary
)


# ==========================================================
# 36. EMT group counts
# ==========================================================

emt_group_counts <- as.data.frame(
  table(
    seurat_epithelial$EMT_GMM_cluster
  )
)


colnames(
  emt_group_counts
) <- c(
  "EMT_group",
  "cell_count"
)


emt_group_counts$fraction <-
  
  emt_group_counts$cell_count /
  sum(
    emt_group_counts$cell_count
  )


write.csv(
  
  emt_group_counts,
  
  file.path(
    tables_dir,
    "EMT_group_counts_GSE183904.csv"
  ),
  
  row.names =
    FALSE
  
)


cat(
  "\nEMT group counts:\n"
)


print(
  emt_group_counts
)


# ==========================================================
# 37. EMT rank counts
# ==========================================================

emt_rank_counts <- as.data.frame(
  table(
    seurat_epithelial$EMT_GMM_rank
  )
)


colnames(
  emt_rank_counts
) <- c(
  "EMT_rank",
  "cell_count"
)


write.csv(
  
  emt_rank_counts,
  
  file.path(
    tables_dir,
    "EMT_rank_counts_GSE183904.csv"
  ),
  
  row.names =
    FALSE
  
)


############################################################
#
#                    EMT FIGURES
#
############################################################


# ==========================================================
# 38. EMT five-rank UMAP
# ==========================================================

p_rank <- DimPlot(
  
  seurat_epithelial,
  
  reduction =
    "umap",
  
  group.by =
    "EMT_GMM_rank",
  
  raster =
    FALSE
  
) +
  ggtitle(
    "Epithelial Cells — Five EMT States"
  )


ggsave(
  
  filename = file.path(
    plots_dir,
    "EMT_04_five_state_GMM_UMAP.pdf"
  ),
  
  plot =
    p_rank,
  
  width =
    10,
  
  height =
    8
  
)


# ==========================================================
# 39. EMT low/intermediate/high UMAP
# ==========================================================

p_group <- DimPlot(
  
  seurat_epithelial,
  
  reduction =
    "umap",
  
  group.by =
    "EMT_GMM_cluster",
  
  raster =
    FALSE
  
) +
  ggtitle(
    "Epithelial Cells — EMT State"
  )


ggsave(
  
  filename = file.path(
    plots_dir,
    "EMT_05_low_intermediate_high_UMAP.pdf"
  ),
  
  plot =
    p_group,
  
  width =
    10,
  
  height =
    8
  
)


# ==========================================================
# 40. EMT score FeaturePlots
# ==========================================================

p_scores <- FeaturePlot(
  
  seurat_epithelial,
  
  features = c(
    "EMT_EPI_z",
    "EMT_MES_z",
    "EMT_Cscore_z"
  ),
  
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
    "EMT_06_score_featureplots.pdf"
  ),
  
  plot =
    p_scores,
  
  width =
    15,
  
  height =
    5
  
)


# ==========================================================
# 41. Canonical EMT marker FeaturePlots
# ==========================================================

marker_features <- intersect(
  
  c(
    "EPCAM",
    "CDH1",
    "VIM",
    "FN1"
  ),
  
  rownames(
    seurat_epithelial
  )
  
)


if (
  length(
    marker_features
  ) > 0
) {
  
  p_markers <- FeaturePlot(
    
    seurat_epithelial,
    
    features =
      marker_features,
    
    reduction =
      "umap",
    
    raster =
      TRUE,
    
    ncol =
      2
    
  )
  
  
  ggsave(
    
    filename = file.path(
      plots_dir,
      "EMT_07_canonical_marker_featureplots.pdf"
    ),
    
    plot =
      p_markers,
    
    width =
      10,
    
    height =
      9
    
  )
  
}


# ==========================================================
# 42. Canonical EMT marker violin plots
# ==========================================================

if (
  length(
    marker_features
  ) > 0
) {
  
  p_marker_violin <- VlnPlot(
    
    seurat_epithelial,
    
    features =
      marker_features,
    
    group.by =
      "EMT_GMM_cluster",
    
    pt.size =
      0,
    
    ncol =
      2
    
  )
  
  
  ggsave(
    
    filename = file.path(
      plots_dir,
      "EMT_08_canonical_marker_violin.pdf"
    ),
    
    plot =
      p_marker_violin,
    
    width =
      10,
    
    height =
      8
    
  )
  
}


# ==========================================================
# 43. EMT C-score violin plot
# ==========================================================

p_cscore <- VlnPlot(
  
  seurat_epithelial,
  
  features =
    "EMT_Cscore_z",
  
  group.by =
    "EMT_GMM_cluster",
  
  pt.size =
    0
  
) +
  ggtitle(
    "EMT C-score by EMT State"
  )


ggsave(
  
  filename = file.path(
    plots_dir,
    "EMT_09_Cscore_by_state.pdf"
  ),
  
  plot =
    p_cscore,
  
  width =
    8,
  
  height =
    6
  
)


# ==========================================================
# 44. MES-versus-EPI score plot
# ==========================================================

score_dataframe <- data.frame(
  
  EMT_MES_z =
    seurat_epithelial$EMT_MES_z,
  
  EMT_EPI_z =
    seurat_epithelial$EMT_EPI_z,
  
  EMT_Cscore_z =
    seurat_epithelial$EMT_Cscore_z,
  
  EMT_rank =
    seurat_epithelial$EMT_GMM_rank,
  
  EMT_group =
    seurat_epithelial$EMT_GMM_cluster
  
)


p_score_scatter <- ggplot(
  
  score_dataframe,
  
  aes(
    x = EMT_EPI_z,
    y = EMT_MES_z,
    group = EMT_group
  )
  
) +
  
  geom_point(
    alpha = 0.25,
    size = 0.4
  ) +
  
  facet_wrap(
    ~ EMT_group
  ) +
  
  labs(
    x = "Epithelial score (z)",
    y = "Mesenchymal score (z)",
    title = "EMT Score Landscape"
  ) +
  
  theme_classic()


ggsave(
  
  filename = file.path(
    plots_dir,
    "EMT_10_MES_vs_EPI_scores.pdf"
  ),
  
  plot =
    p_score_scatter,
  
  width =
    11,
  
  height =
    6
  
)


############################################################
#
#              SAVE EPITHELIAL OBJECT
#
############################################################


# ==========================================================
# 45. Save five-component GMM model
# ==========================================================

gmm_model_file <- file.path(
  objects_dir,
  "EMT_GMM_model_GSE183904.rds"
)


saveRDS(
  gmm_fit,
  gmm_model_file
)


# ==========================================================
# 46. Save epithelial EMT object
# ==========================================================

epithelial_output_file <- file.path(
  
  objects_dir,
  
  "seurat_epithelial_EMT_GSE183904.rds"
  
)


saveRDS(
  
  seurat_epithelial,
  
  epithelial_output_file
  
)


log_message(
  "Saved epithelial EMT object:",
  epithelial_output_file
)


############################################################
#
#       TRANSFER EMT RESULTS TO GLOBAL OBJECT
#
############################################################


# ==========================================================
# 47. Confirm epithelial cell correspondence
# ==========================================================

epithelial_cell_names <- colnames(
  seurat_epithelial
)


global_cell_names <- colnames(
  seurat_global
)


cell_match <- match(
  epithelial_cell_names,
  global_cell_names
)


checkpoint(
  all(
    !is.na(
      cell_match
    )
  ),
  "Some epithelial cells could not be matched to global object."
)


# ==========================================================
# 48. Metadata columns to transfer
# ==========================================================

metadata_to_transfer <- c(
  
  "EMT_MES_AddModule",
  "EMT_EPI_AddModule",
  
  "EMT_MES_mean",
  "EMT_EPI_mean",
  
  "EMT_MES",
  "EMT_EPI",
  
  "EMT_Cscore",
  
  "EMT_MES_z",
  "EMT_EPI_z",
  "EMT_Cscore_z",
  
  "EMT_GMM_component",
  "EMT_GMM_rank",
  "EMT_GMM_cluster",
  "EMT_GMM_confidence"
  
)


# ==========================================================
# 49. Transfer EMT metadata
# ==========================================================

for (
  metadata_column
  in metadata_to_transfer
) {
  
  source_vector <-
    seurat_epithelial@meta.data[
      [metadata_column]
    ]
  
  
  if (
    is.numeric(
      source_vector
    ) ||
    is.integer(
      source_vector
    )
  ) {
    
    seurat_global@meta.data[
      [metadata_column]
    ] <- NA_real_
    
  } else {
    
    seurat_global@meta.data[
      [metadata_column]
    ] <- NA_character_
    
  }
  
  
  seurat_global@meta.data[
    cell_match,
    metadata_column
  ] <- as.vector(
    source_vector
  )
  
}


# ==========================================================
# 50. Restore ordered EMT factors
# ==========================================================

seurat_global$EMT_GMM_rank <- factor(
  
  seurat_global$EMT_GMM_rank,
  
  levels = c(
    "EMT_rank1",
    "EMT_rank2",
    "EMT_rank3",
    "EMT_rank4",
    "EMT_rank5"
  )
  
)


seurat_global$EMT_GMM_cluster <- factor(
  
  seurat_global$EMT_GMM_cluster,
  
  levels = c(
    "EMT_low",
    "EMT_intermediate",
    "EMT_high"
  )
  
)


# ==========================================================
# 51. Create cell_type_final
#
# Non-epithelial cells retain their original cell type.
# Epithelial cells become:
#
#   EMT_low
#   EMT_intermediate
#   EMT_high
# ==========================================================

seurat_global$cell_type_final <-
  
  as.character(
    seurat_global$cell_type
  )


global_epithelial_cells <- which(
  seurat_global$cell_type ==
    "Epithelial"
)


checkpoint(
  length(
    global_epithelial_cells
  ) ==
    ncol(
      seurat_epithelial
    ),
  paste0(
    "Mismatch between epithelial-cell counts in ",
    "global and epithelial objects."
  )
)


seurat_global$cell_type_final[
  global_epithelial_cells
] <- as.character(
  
  seurat_global$EMT_GMM_cluster[
    global_epithelial_cells
  ]
  
)


seurat_global$cell_type_final <- factor(
  seurat_global$cell_type_final
)


# ==========================================================
# 52. Final global summary
# ==========================================================

cat(
  "\nFinal EMT status in global object:\n"
)


print(
  table(
    seurat_global$EMT_GMM_cluster,
    useNA = "ifany"
  )
)


cat(
  "\nFinal cell types:\n"
)


print(
  table(
    seurat_global$cell_type_final,
    useNA = "ifany"
  )
)


# ==========================================================
# 53. Save final global object
# ==========================================================

final_output_file <- file.path(
  
  objects_dir,
  
  "seurat_final_GSE183904.rds"
  
)


saveRDS(
  
  seurat_global,
  
  final_output_file
  
)


log_message(
  "Saved final global object:",
  final_output_file
)


# ==========================================================
# 54. Save session information
# ==========================================================

capture.output(
  
  sessionInfo(),
  
  file = file.path(
    results_dir,
    "sessionInfo_04_EMT_scoring.txt"
  )
  
)


# ==========================================================
# 55. Completion message
# ==========================================================

cat(
  
  "\n",
  "=====================================================\n",
  " EMT SCORING COMPLETED SUCCESSFULLY\n",
  "=====================================================\n",
  "\n",
  "Epithelial reclustering resolution:\n",
  epithelial_cluster_resolution,
  "\n\n",
  "GMM components:\n",
  gmm_components,
  "\n\n",
  "Epithelial EMT object:\n",
  epithelial_output_file,
  "\n\n",
  "Final global object:\n",
  final_output_file,
  "\n\n",
  "GMM model:\n",
  gmm_model_file,
  "\n",
  "=====================================================\n"
  
)


log_message(
  "04_emt_scoring.R completed."
)