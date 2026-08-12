############################################################
# Monocle3 EMT Trajectory Analysis
#
# Project:
#   Gastric Cancer EMT–TME Analysis
#
# Dataset:
#   GSE183904
#
# Input:
#   objects/seurat_epithelial_EMT_GSE183904.rds
#
# Outputs:
#   objects/monocle3_EMT_trajectory_GSE183904.rds
#   objects/seurat_epithelial_EMT_trajectory_GSE183904.rds
#
#   results/tables/EMT_trajectory_graph_test.csv
#   results/tables/EMT_trajectory_significant_genes.csv
#   results/tables/EMT_trajectory_gene_modules.csv
#   results/tables/EMT_trajectory_module_expression_by_EMT_status.csv
#   results/tables/EMT_trajectory_module_summary.csv
#   results/tables/EMT_trajectory_root_nodes.csv
#
# Workflow:
#   01. Load epithelial EMT Seurat object
#   02. Build Monocle3 cell_data_set
#   03. Preprocess using 50 PCs
#   04. UMAP dimensionality reduction
#   05. Monocle3 clustering
#   06. Learn principal trajectory graph
#   07. Identify trajectory roots from EMT-low cells
#   08. Order cells in pseudotime
#   09. Identify trajectory-associated genes
#  10. Filter q < 0.05 and Moran's I > 0.1
#  11. Identify trajectory-dependent gene modules
#  12. Aggregate modules across EMT states
#  13. Generate trajectory/module figures
############################################################


# ==========================================================
# 0. Required packages
# ==========================================================

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "monocle3",
  "Matrix",
  "dplyr",
  "tibble",
  "ggplot2",
  "pheatmap",
  "igraph",
  "SummarizedExperiment"
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
  library(monocle3)
  library(Matrix)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(pheatmap)
  library(igraph)
  library(SummarizedExperiment)
  
})


# ==========================================================
# 1. Reproducibility and analysis parameters
# ==========================================================

set.seed(123)


# Methods-defined Monocle dimensionality
monocle_num_dim <- 50


# Trajectory-gene thresholds
trajectory_q_cutoff <- 0.05

trajectory_morans_cutoff <- 0.10


# Gene-module resolution from final Methods
gene_module_resolution <- 0.001


# Number of CPU cores used by graph_test()
#
# Keep this conservative for cross-platform reproducibility.
# Increase locally if desired.

graph_test_cores <- 4


# EMT state used as biological trajectory origin
root_emt_state <- "EMT_low"


# ==========================================================
# 2. Project directories
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
# Retrieve RNA counts from Seurat
#
# Handles SeuratObject v4/v5.
# ----------------------------------------------------------

get_counts_matrix <- function(
    seurat_object
) {
  
  if (
    packageVersion("SeuratObject") >=
    package_version("5.0.0")
  ) {
    
    counts_matrix <- LayerData(
      object = seurat_object,
      assay = "RNA",
      layer = "counts"
    )
    
  } else {
    
    counts_matrix <- GetAssayData(
      object = seurat_object,
      assay = "RNA",
      slot = "counts"
    )
    
  }
  
  
  counts_matrix
  
}


# ==========================================================
# 4. Load epithelial EMT object
# ==========================================================

input_file <- file.path(
  objects_dir,
  "seurat_epithelial_EMT_GSE183904.rds"
)


checkpoint(
  file.exists(input_file),
  paste0(
    "Epithelial EMT object not found:\n",
    input_file,
    "\n\nRun 04_emt_scoring.R first."
  )
)


log_message(
  "Loading epithelial EMT object..."
)


seurat_epithelial <- readRDS(
  input_file
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


checkpoint(
  any(
    as.character(
      seurat_epithelial$EMT_GMM_cluster
    ) == root_emt_state
  ),
  paste0(
    "No ",
    root_emt_state,
    " cells were found."
  )
)


log_message(
  "Epithelial cells:",
  ncol(seurat_epithelial)
)


cat(
  "\nEMT-state distribution:\n"
)


print(
  table(
    seurat_epithelial$EMT_GMM_cluster
  )
)


# ==========================================================
# 5. Extract counts and metadata
# ==========================================================

log_message(
  "Extracting count matrix..."
)


counts_matrix <- get_counts_matrix(
  seurat_epithelial
)


checkpoint(
  nrow(counts_matrix) > 0 &&
    ncol(counts_matrix) > 0,
  "RNA count matrix is empty."
)


# Ensure matrix remains sparse
if (
  !inherits(
    counts_matrix,
    "sparseMatrix"
  )
) {
  
  counts_matrix <- as(
    counts_matrix,
    "dgCMatrix"
  )
  
}


cell_metadata <- seurat_epithelial@meta.data


checkpoint(
  identical(
    rownames(cell_metadata),
    colnames(counts_matrix)
  ),
  "Cell metadata and expression matrix are not aligned."
)


gene_metadata <- data.frame(
  
  gene_short_name =
    rownames(counts_matrix),
  
  row.names =
    rownames(counts_matrix),
  
  stringsAsFactors =
    FALSE
  
)


# ==========================================================
# 6. Construct Monocle3 cell_data_set
# ==========================================================

log_message(
  "Creating Monocle3 cell_data_set..."
)


cds <- new_cell_data_set(
  
  expression_data =
    counts_matrix,
  
  cell_metadata =
    cell_metadata,
  
  gene_metadata =
    gene_metadata
  
)


rm(
  counts_matrix,
  cell_metadata,
  gene_metadata
)


gc()


checkpoint(
  ncol(cds) ==
    ncol(seurat_epithelial),
  "Cell count changed during Monocle conversion."
)


# ==========================================================
# 7. Monocle3 preprocessing
#
# Final Methods specify 50 PCs.
# ==========================================================

log_message(
  "Preprocessing Monocle3 CDS using 50 dimensions..."
)


cds <- preprocess_cds(
  
  cds,
  
  method =
    "PCA",
  
  num_dim =
    monocle_num_dim
  
)


# ==========================================================
# 8. Monocle PC variance plot
# ==========================================================

p_pc_variance <- plot_pc_variance_explained(
  cds
)


ggsave(
  
  filename = file.path(
    plots_dir,
    "EMT_trajectory_01_pc_variance.pdf"
  ),
  
  plot =
    p_pc_variance,
  
  width =
    7,
  
  height =
    5
  
)


# ==========================================================
# 9. UMAP dimensionality reduction
# ==========================================================

log_message(
  "Running Monocle3 UMAP..."
)


set.seed(123)


cds <- reduce_dimension(
  
  cds,
  
  reduction_method =
    "UMAP"
  
)


# ==========================================================
# 10. Cluster cells in Monocle3
# ==========================================================

log_message(
  "Clustering epithelial cells in Monocle3..."
)


set.seed(123)


cds <- cluster_cells(
  
  cds,
  
  reduction_method =
    "UMAP"
  
)


# Store Monocle clustering information explicitly
SummarizedExperiment::colData(cds)$monocle_cluster <-
  
  as.character(
    clusters(cds)
  )


SummarizedExperiment::colData(cds)$monocle_partition <-
  
  as.character(
    partitions(cds)
  )


cat(
  "\nMonocle partitions:\n"
)


print(
  table(
    SummarizedExperiment::colData(cds)$monocle_partition
  )
)


# ==========================================================
# 11. Plot Monocle partitions before trajectory learning
# ==========================================================

p_partitions <- plot_cells(
  
  cds,
  
  color_cells_by =
    "monocle_partition",
  
  label_cell_groups =
    FALSE,
  
  label_leaves =
    FALSE,
  
  label_branch_points =
    FALSE,
  
  show_trajectory_graph =
    FALSE
  
)


ggsave(
  
  filename = file.path(
    plots_dir,
    "EMT_trajectory_02_monocle_partitions.pdf"
  ),
  
  plot =
    p_partitions,
  
  width =
    8,
  
  height =
    7
  
)


# ==========================================================
# 12. Learn principal trajectory graph
# ==========================================================

log_message(
  "Learning principal trajectory graph..."
)


set.seed(123)


cds <- learn_graph(
  cds
)


checkpoint(
  !is.null(
    principal_graph(cds)[["UMAP"]]
  ),
  "Monocle principal graph was not generated."
)


# ==========================================================
# 13. Select root principal nodes from EMT-low cells
#
# Instead of interactive order_cells(), root nodes are
# selected programmatically from graph vertices most highly
# occupied by EMT-low cells.
#
# If multiple Monocle partitions contain EMT-low cells,
# one EMT-low root is selected for each such partition.
# ==========================================================

get_emt_root_nodes <- function(
    cds,
    emt_column = "EMT_GMM_cluster",
    root_state = "EMT_low"
) {
  
  checkpoint(
    emt_column %in%
      colnames(
        SummarizedExperiment::colData(cds)
      ),
    paste0(
      emt_column,
      " is missing from Monocle metadata."
    )
  )
  
  
  emt_status <- as.character(
    SummarizedExperiment::colData(cds)[
      [emt_column]
    ]
  )
  
  
  names(emt_status) <- colnames(
    cds
  )
  
  
  root_cells <- names(
    emt_status
  )[
    emt_status ==
      root_state
  ]
  
  
  checkpoint(
    length(root_cells) > 0,
    paste0(
      "No cells belong to ",
      root_state,
      "."
    )
  )
  
  
  partition_vector <- partitions(
    cds
  )
  
  
  partition_vector <- partition_vector[
    colnames(cds)
  ]
  
  
  closest_vertex <-
    cds@principal_graph_aux[
      ["UMAP"]
    ]$pr_graph_cell_proj_closest_vertex
  
  
  closest_vertex <- as.matrix(
    closest_vertex[
      colnames(cds),
      ,
      drop = FALSE
    ]
  )
  
  
  checkpoint(
    nrow(closest_vertex) ==
      ncol(cds),
    "Could not align cells to principal graph vertices."
  )
  
  
  root_partitions <- unique(
    as.character(
      partition_vector[
        root_cells
      ]
    )
  )
  
  
  selected_root_nodes <- character(
    0
  )
  
  
  root_summary <- list()
  
  
  for (
    current_partition
    in root_partitions
  ) {
    
    partition_root_cells <- root_cells[
      as.character(
        partition_vector[
          root_cells
        ]
      ) ==
        current_partition
    ]
    
    
    vertex_table <- table(
      
      closest_vertex[
        partition_root_cells,
        1
      ]
      
    )
    
    
    checkpoint(
      length(vertex_table) > 0,
      paste0(
        "No graph vertices were found for EMT-low cells ",
        "in partition ",
        current_partition,
        "."
      )
    )
    
    
    selected_vertex_index <-
      as.numeric(
        names(
          which.max(
            vertex_table
          )
        )
      )
    
    
    selected_node <-
      igraph::V(
        principal_graph(cds)[
          ["UMAP"]
        ]
      )$name[
        selected_vertex_index
      ]
    
    
    selected_root_nodes <- c(
      selected_root_nodes,
      selected_node
    )
    
    
    root_summary[
      [length(root_summary) + 1]
    ] <- data.frame(
      
      partition =
        current_partition,
      
      root_principal_node =
        selected_node,
      
      EMT_low_cells_in_partition =
        length(
          partition_root_cells
        ),
      
      EMT_low_cells_at_selected_node =
        max(
          vertex_table
        ),
      
      stringsAsFactors =
        FALSE
      
    )
    
  }
  
  
  list(
    
    root_nodes =
      unique(
        selected_root_nodes
      ),
    
    summary =
      do.call(
        rbind,
        root_summary
      )
    
  )
  
}


root_result <- get_emt_root_nodes(
  
  cds,
  
  emt_column =
    "EMT_GMM_cluster",
  
  root_state =
    root_emt_state
  
)


root_nodes <- root_result$root_nodes


root_summary <- root_result$summary


checkpoint(
  length(root_nodes) > 0,
  "No trajectory root node was selected."
)


cat(
  "\nSelected EMT-low root principal nodes:\n"
)


print(
  root_summary
)


write.csv(
  
  root_summary,
  
  file.path(
    tables_dir,
    "EMT_trajectory_root_nodes.csv"
  ),
  
  row.names =
    FALSE
  
)


# ==========================================================
# 14. Order cells in pseudotime
# ==========================================================

log_message(
  "Ordering cells from EMT-low root nodes..."
)


cds <- order_cells(
  
  cds,
  
  reduction_method =
    "UMAP",
  
  root_pr_nodes =
    root_nodes
  
)


# ==========================================================
# 15. Extract pseudotime
# ==========================================================

cell_pseudotime <- pseudotime(
  cds
)


checkpoint(
  length(cell_pseudotime) ==
    ncol(cds),
  "Pseudotime vector has incorrect length."
)


SummarizedExperiment::colData(cds)$pseudotime <-
  
  as.numeric(
    cell_pseudotime[
      colnames(cds)
    ]
  )


finite_pseudotime <- is.finite(
  SummarizedExperiment::colData(cds)$pseudotime
)


n_finite <- sum(
  finite_pseudotime
)


n_infinite <- sum(
  !finite_pseudotime
)


log_message(
  "Cells with finite pseudotime:",
  n_finite
)


log_message(
  "Cells with infinite pseudotime:",
  n_infinite
)


if (n_infinite > 0) {
  
  warning(
    paste0(
      n_infinite,
      " cells have infinite pseudotime. ",
      "These cells occur in graph regions/partitions ",
      "not reachable from an EMT-low root."
    ),
    call. = FALSE
  )
  
}


pseudotime_summary <- data.frame(
  
  total_cells =
    ncol(cds),
  
  finite_pseudotime_cells =
    n_finite,
  
  infinite_pseudotime_cells =
    n_infinite,
  
  fraction_finite =
    n_finite /
    ncol(cds)
  
)


write.csv(
  
  pseudotime_summary,
  
  file.path(
    tables_dir,
    "EMT_trajectory_pseudotime_summary.csv"
  ),
  
  row.names =
    FALSE
  
)


# ==========================================================
# 16. Plot trajectory by EMT state
# ==========================================================

p_emt_trajectory <- plot_cells(
  
  cds,
  
  color_cells_by =
    "EMT_GMM_cluster",
  
  label_cell_groups =
    FALSE,
  
  label_leaves =
    FALSE,
  
  label_branch_points =
    FALSE,
  
  show_trajectory_graph =
    TRUE,
  
  graph_label_size =
    1.5
  
)


ggsave(
  
  filename = file.path(
    plots_dir,
    "EMT_trajectory_03_graph_by_EMT_state.pdf"
  ),
  
  plot =
    p_emt_trajectory,
  
  width =
    8,
  
  height =
    7
  
)


# ==========================================================
# 17. Plot trajectory by five EMT ranks
# ==========================================================

if (
  "EMT_GMM_rank" %in%
  colnames(
    SummarizedExperiment::colData(cds)
  )
) {
  
  p_rank_trajectory <- plot_cells(
    
    cds,
    
    color_cells_by =
      "EMT_GMM_rank",
    
    label_cell_groups =
      FALSE,
    
    label_leaves =
      FALSE,
    
    label_branch_points =
      FALSE,
    
    show_trajectory_graph =
      TRUE
    
  )
  
  
  ggsave(
    
    filename = file.path(
      plots_dir,
      "EMT_trajectory_04_graph_by_five_EMT_ranks.pdf"
    ),
    
    plot =
      p_rank_trajectory,
    
    width =
      8,
    
    height =
      7
    
  )
  
}


# ==========================================================
# 18. Plot pseudotime
# ==========================================================

p_pseudotime <- plot_cells(
  
  cds,
  
  color_cells_by =
    "pseudotime",
  
  label_cell_groups =
    FALSE,
  
  label_leaves =
    FALSE,
  
  label_branch_points =
    FALSE,
  
  show_trajectory_graph =
    TRUE
  
)


ggsave(
  
  filename = file.path(
    plots_dir,
    "EMT_trajectory_05_pseudotime.pdf"
  ),
  
  plot =
    p_pseudotime,
  
  width =
    8,
  
  height =
    7
  
)


############################################################
#
#             TRAJECTORY-ASSOCIATED GENES
#
############################################################


# ==========================================================
# 19. Graph-autocorrelation test
#
# principal_graph tests expression association with
# positions along the learned trajectory.
# ==========================================================

log_message(
  "Running graph_test on principal trajectory graph..."
)


trajectory_graph_test <- graph_test(
  
  cds,
  
  neighbor_graph =
    "principal_graph",
  
  cores =
    graph_test_cores
  
)


trajectory_graph_test_df <- as.data.frame(
  trajectory_graph_test
)


trajectory_graph_test_df$gene <- rownames(
  trajectory_graph_test_df
)


trajectory_graph_test_df <- trajectory_graph_test_df[
  ,
  c(
    "gene",
    setdiff(
      colnames(
        trajectory_graph_test_df
      ),
      "gene"
    )
  ),
  drop = FALSE
]


write.csv(
  
  trajectory_graph_test_df,
  
  file.path(
    tables_dir,
    "EMT_trajectory_graph_test.csv"
  ),
  
  row.names =
    FALSE
  
)


# ==========================================================
# 20. Select significant trajectory genes
#
# Final Methods:
#   q-value < 0.05
#   Moran's I > 0.1
# ==========================================================

trajectory_significant <- trajectory_graph_test_df %>%
  
  dplyr::filter(
    
    !is.na(q_value),
    
    !is.na(morans_I),
    
    q_value <
      trajectory_q_cutoff,
    
    morans_I >
      trajectory_morans_cutoff
    
  ) %>%
  
  dplyr::arrange(
    q_value,
    dplyr::desc(
      morans_I
    )
  )


checkpoint(
  nrow(
    trajectory_significant
  ) > 0,
  paste0(
    "No genes passed q < ",
    trajectory_q_cutoff,
    " and Moran's I > ",
    trajectory_morans_cutoff,
    "."
  )
)


log_message(
  "Significant trajectory-associated genes:",
  nrow(
    trajectory_significant
  )
)


write.csv(
  
  trajectory_significant,
  
  file.path(
    tables_dir,
    "EMT_trajectory_significant_genes.csv"
  ),
  
  row.names =
    FALSE
  
)


trajectory_gene_ids <- trajectory_significant$gene


trajectory_gene_ids <- intersect(
  
  trajectory_gene_ids,
  
  rownames(
    cds
  )
  
)


checkpoint(
  length(
    trajectory_gene_ids
  ) > 1,
  "Too few trajectory genes remain for module discovery."
)


############################################################
#
#                  GENE MODULES
#
############################################################


# ==========================================================
# 21. Ensure gene names are available to Monocle
# ==========================================================

SummarizedExperiment::rowData(
  cds
)$gene_short_name <- rownames(
  cds
)


# ==========================================================
# 22. Identify trajectory-dependent gene modules
#
# Final Methods:
#   resolution = 0.001
# ==========================================================

log_message(
  "Identifying trajectory-dependent gene modules..."
)


set.seed(123)


gene_module_df <- find_gene_modules(
  
  cds[
    trajectory_gene_ids,
    ,
    drop = FALSE
  ],
  
  reduction_method =
    "UMAP",
  
  resolution =
    gene_module_resolution
  
)


checkpoint(
  nrow(gene_module_df) > 0,
  "find_gene_modules() returned no genes."
)


checkpoint(
  "module" %in%
    colnames(
      gene_module_df
    ),
  "Module column missing from find_gene_modules output."
)


# ----------------------------------------------------------
# Standardize gene column for easier downstream use
# ----------------------------------------------------------

if (
  "id" %in%
  colnames(
    gene_module_df
  )
) {
  
  gene_module_df$gene <-
    gene_module_df$id
  
} else if (
  "gene_id" %in%
  colnames(
    gene_module_df
  )
) {
  
  gene_module_df$gene <-
    gene_module_df$gene_id
  
} else {
  
  gene_module_df$gene <-
    rownames(
      gene_module_df
    )
  
}


gene_module_df <- gene_module_df[
  ,
  c(
    "gene",
    setdiff(
      colnames(
        gene_module_df
      ),
      "gene"
    )
  ),
  drop = FALSE
]


write.csv(
  
  gene_module_df,
  
  file.path(
    tables_dir,
    "EMT_trajectory_gene_modules.csv"
  ),
  
  row.names =
    FALSE
  
)


# ==========================================================
# 23. Module-size summary
# ==========================================================

module_summary <- gene_module_df %>%
  
  dplyr::count(
    module,
    name = "n_genes"
  ) %>%
  
  dplyr::arrange(
    module
  )


write.csv(
  
  module_summary,
  
  file.path(
    tables_dir,
    "EMT_trajectory_module_summary.csv"
  ),
  
  row.names =
    FALSE
  
)


cat(
  "\nTrajectory module sizes:\n"
)


print(
  module_summary
)


# ==========================================================
# 24. Aggregate module expression by EMT state
# ==========================================================

checkpoint(
  "EMT_GMM_cluster" %in%
    colnames(
      SummarizedExperiment::colData(cds)
    ),
  "EMT_GMM_cluster missing from CDS."
)


emt_group <- factor(
  
  as.character(
    SummarizedExperiment::colData(cds)$EMT_GMM_cluster
  ),
  
  levels = c(
    "EMT_low",
    "EMT_intermediate",
    "EMT_high"
  )
  
)


cell_group_df <- tibble::tibble(
  
  cell =
    rownames(
      SummarizedExperiment::colData(cds)
    ),
  
  cell_group =
    emt_group
  
)


module_expression <- aggregate_gene_expression(
  
  cds,
  
  gene_module_df,
  
  cell_group_df
  
)


rownames(
  module_expression
) <- paste0(
  "Module ",
  rownames(
    module_expression
  )
)


# Force biologically meaningful EMT-state column order
desired_emt_order <- c(
  "EMT_low",
  "EMT_intermediate",
  "EMT_high"
)


available_emt_columns <- intersect(
  desired_emt_order,
  colnames(
    module_expression
  )
)


module_expression <- module_expression[
  ,
  available_emt_columns,
  drop = FALSE
]


checkpoint(
  ncol(
    module_expression
  ) >= 2,
  "Insufficient EMT groups for module aggregation."
)


write.csv(
  
  module_expression,
  
  file.path(
    tables_dir,
    "EMT_trajectory_module_expression_by_EMT_status.csv"
  ),
  
  row.names =
    TRUE
  
)


# ==========================================================
# 25. Module heatmap
# ==========================================================

module_row_sd <- apply(
  
  module_expression,
  
  1,
  
  sd,
  
  na.rm = TRUE
  
)


heatmap_matrix <- module_expression[
  is.finite(
    module_row_sd
  ) &
    module_row_sd > 0,
  ,
  drop = FALSE
]


checkpoint(
  nrow(
    heatmap_matrix
  ) > 0,
  "No variable modules remain for heatmap."
)


heatmap_file <- file.path(
  plots_dir,
  "EMT_trajectory_06_modules_by_EMT_state_heatmap.pdf"
)


pdf(
  
  heatmap_file,
  
  width =
    7,
  
  height =
    max(
      6,
      0.25 *
        nrow(
          heatmap_matrix
        )
    )
  
)


pheatmap::pheatmap(
  
  heatmap_matrix,
  
  scale =
    "row",
  
  cluster_rows =
    TRUE,
  
  cluster_cols =
    FALSE,
  
  clustering_method =
    "ward.D2",
  
  main =
    "Trajectory Gene Modules by EMT State",
  
  border_color =
    NA
  
)


dev.off()


log_message(
  "Module heatmap saved:",
  heatmap_file
)


# ==========================================================
# 26. Identify modules increasing across EMT states
#
# This is saved as a descriptive result only.
#
# We will review this table before deciding which modules
# should be passed to NicheNet in script 06.
# ==========================================================

if (
  all(
    c(
      "EMT_low",
      "EMT_intermediate",
      "EMT_high"
    ) %in%
    colnames(
      module_expression
    )
  )
) {
  
  module_trend <- data.frame(
    
    module =
      rownames(
        module_expression
      ),
    
    EMT_low =
      module_expression[
        ,
        "EMT_low"
      ],
    
    EMT_intermediate =
      module_expression[
        ,
        "EMT_intermediate"
      ],
    
    EMT_high =
      module_expression[
        ,
        "EMT_high"
      ],
    
    stringsAsFactors =
      FALSE
    
  )
  
  
  module_trend$high_minus_low <-
    
    module_trend$EMT_high -
    module_trend$EMT_low
  
  
  module_trend$monotonic_increase <-
    
    module_trend$EMT_low <=
    module_trend$EMT_intermediate &
    module_trend$EMT_intermediate <=
    module_trend$EMT_high
  
  
  module_trend <- module_trend[
    order(
      -module_trend$high_minus_low
    ),
  ]
  
  
  write.csv(
    
    module_trend,
    
    file.path(
      tables_dir,
      "EMT_trajectory_module_trends.csv"
    ),
    
    row.names =
      FALSE
    
  )
  
}


############################################################
#
#                 SAVE OBJECTS
#
############################################################


# ==========================================================
# 27. Save Monocle3 CDS
# ==========================================================

cds_output_file <- file.path(
  
  objects_dir,
  
  "monocle3_EMT_trajectory_GSE183904.rds"
  
)


saveRDS(
  
  cds,
  
  cds_output_file
  
)


log_message(
  "Saved Monocle3 trajectory object:",
  cds_output_file
)


# ==========================================================
# 28. Transfer pseudotime back to epithelial Seurat object
# ==========================================================

seurat_epithelial$Monocle3_pseudotime <- NA_real_


matched_cells <- match(
  
  colnames(
    seurat_epithelial
  ),
  
  colnames(
    cds
  )
  
)


checkpoint(
  all(
    !is.na(
      matched_cells
    )
  ),
  "Some epithelial cells could not be matched back to CDS."
)


seurat_epithelial$Monocle3_pseudotime <-
  
  SummarizedExperiment::colData(cds)$pseudotime[
    matched_cells
  ]


seurat_epithelial$Monocle3_cluster <-
  
  SummarizedExperiment::colData(cds)$monocle_cluster[
    matched_cells
  ]


seurat_epithelial$Monocle3_partition <-
  
  SummarizedExperiment::colData(cds)$monocle_partition[
    matched_cells
  ]


seurat_trajectory_file <- file.path(
  
  objects_dir,
  
  "seurat_epithelial_EMT_trajectory_GSE183904.rds"
  
)


saveRDS(
  
  seurat_epithelial,
  
  seurat_trajectory_file
  
)


log_message(
  "Saved epithelial Seurat object with pseudotime:",
  seurat_trajectory_file
)


# ==========================================================
# 29. Save session information
# ==========================================================

capture.output(
  
  sessionInfo(),
  
  file = file.path(
    results_dir,
    "sessionInfo_05_monocle3_trajectory.txt"
  )
  
)


# ==========================================================
# 30. Completion summary
# ==========================================================

cat(
  
  "\n",
  "=====================================================\n",
  " MONOCLE3 EMT TRAJECTORY COMPLETED\n",
  "=====================================================\n",
  "\n",
  "Root EMT state:\n",
  root_emt_state,
  "\n\n",
  "Root principal nodes:\n",
  paste(
    root_nodes,
    collapse = ", "
  ),
  "\n\n",
  "Cells with finite pseudotime:\n",
  n_finite,
  " / ",
  ncol(cds),
  "\n\n",
  "Significant trajectory genes:\n",
  nrow(
    trajectory_significant
  ),
  "\n\n",
  "Number of gene modules:\n",
  length(
    unique(
      gene_module_df$module
    )
  ),
  "\n\n",
  "Monocle CDS:\n",
  cds_output_file,
  "\n\n",
  "Seurat trajectory object:\n",
  seurat_trajectory_file,
  "\n",
  "=====================================================\n"
  
)


log_message(
  "05_emt_trajectory_monocle3.R completed."
)