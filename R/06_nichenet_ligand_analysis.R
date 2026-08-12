############################################################
# NicheNet T/NK-to-EMT Ligand Activity Analysis
#
# Project:
#   Gastric Cancer EMT–TME Analysis
#
# Dataset:
#   GSE183904
#
# Biological question:
#   Which ligands expressed by T/NK cells could contribute
#   to EMT-associated transcriptional programs in
#   EMT-low epithelial cells?
#
# Sender:
#   T/NK
#
# Receiver:
#   EMT_low
#
# Inputs:
#   objects/seurat_final_GSE183904.rds
#
#   results/tables/
#     EMT_trajectory_gene_modules.csv
#     EMT_trajectory_module_trends.csv
#
# NicheNet resources:
#   data/nichenet/
#
# Outputs:
#   results/tables/
#     NicheNet_expressed_genes_receiver.csv
#     NicheNet_expressed_genes_sender.csv
#     NicheNet_expressed_receptors.csv
#     NicheNet_potential_TNK_ligands.csv
#     NicheNet_selected_trajectory_modules.csv
#     NicheNet_geneset_of_interest.csv
#     NicheNet_ligand_activity_all.csv
#     NicheNet_top30_ligands.csv
#     NicheNet_active_ligand_target_links.csv
#     NicheNet_ligand_receptor_links.csv
#     NicheNet_signaling_GRN_edges.csv
#     NicheNet_top_TFs.csv
#
# Main figures:
#   NicheNet_01_ligand_activity_distribution.pdf
#   NicheNet_02_top_ligand_activity_heatmap.pdf
#   NicheNet_03_ligand_target_heatmap.pdf
#   NicheNet_04_ligand_receptor_heatmap.pdf
#
# Downstream:
#   NicheNet_top_TFs.csv is used by:
#   07_emt_lncrna_analysis.R
############################################################


# ==========================================================
# 0. Required packages
# ==========================================================

required_packages <- c(
  "nichenetr",
  "Seurat",
  "SeuratObject",
  "dplyr",
  "tibble",
  "tidyr",
  "ggplot2"
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
  
  library(nichenetr)
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  
})


# ==========================================================
# 1. Reproducibility and analysis parameters
# ==========================================================

set.seed(123)


# ----------------------------------------------------------
# Sender / receiver definition
# ----------------------------------------------------------

sender_celltype <- "T/NK"

receiver_celltype <- "EMT_low"


# ----------------------------------------------------------
# Expression threshold
#
# A gene is considered expressed when detected in at least
# 10% of cells in the relevant population.
# ----------------------------------------------------------

expression_pct <- 0.10


# ----------------------------------------------------------
# Ligand prioritization
# ----------------------------------------------------------

number_top_ligands <- 30


# Number of highest-priority targets extracted per ligand
number_targets_per_ligand <- 50


# Visualization cutoff for ligand-target regulatory
# potential
ligand_target_visualization_cutoff <- 0.40


# ----------------------------------------------------------
# Signaling-path analysis
# ----------------------------------------------------------

number_ligands_for_path <- 10

number_targets_for_path <- 100

top_n_regulators <- 4

number_top_tfs <- 50


# ----------------------------------------------------------
# Trajectory-module selection
#
# NULL:
#   automatically identify modules that increase from
#   EMT_low -> EMT_intermediate -> EMT_high.
#
# If an exact historical set needs to be reproduced later,
# it can instead be specified manually, e.g.:
#
# manual_emt_modules <- c("1", "2", "8")
#
# We deliberately do NOT hard-code historical module IDs
# because module numbering can change when Monocle is rerun.
# ----------------------------------------------------------

manual_emt_modules <- NULL


# ==========================================================
# 2. Project directories
# ==========================================================

base_dir <- "."


data_dir <- file.path(
  base_dir,
  "data"
)


nichenet_data_dir <- file.path(
  data_dir,
  "nichenet"
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
  nichenet_data_dir,
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
# Download external resource only if missing
# ----------------------------------------------------------

download_resource <- function(
    url,
    destination
) {
  
  if (!file.exists(destination)) {
    
    log_message(
      "Downloading:",
      basename(destination)
    )
    
    
    old_timeout <- getOption(
      "timeout"
    )
    
    
    options(
      timeout = max(
        1200,
        old_timeout
      )
    )
    
    
    on.exit(
      options(
        timeout = old_timeout
      ),
      add = TRUE
    )
    
    
    download.file(
      url = url,
      destfile = destination,
      mode = "wb"
    )
    
  } else {
    
    log_message(
      "Resource already present:",
      basename(destination)
    )
    
  }
  
  
  checkpoint(
    file.exists(destination),
    paste0(
      "Required resource could not be obtained:\n",
      destination
    )
  )
  
}


# ----------------------------------------------------------
# Convert old HGNC aliases when possible while preserving
# original symbols if no replacement is returned.
# ----------------------------------------------------------

standardize_symbols <- function(
    genes
) {
  
  genes <- unique(
    as.character(
      genes
    )
  )
  
  
  genes <- genes[
    !is.na(genes) &
      genes != ""
  ]
  
  
  converted <- try(
    
    nichenetr::convert_alias_to_symbols(
      genes,
      organism = "human",
      verbose = FALSE
    ),
    
    silent = TRUE
    
  )
  
  
  if (
    inherits(
      converted,
      "try-error"
    )
  ) {
    
    return(
      genes
    )
    
  }
  
  
  converted <- as.character(
    converted
  )
  
  
  invalid <- is.na(converted) |
    converted == ""
  
  
  converted[invalid] <-
    genes[invalid]
  
  
  unique(
    converted
  )
  
}


# ==========================================================
# 4. Define NicheNet v2 resources
# ==========================================================

resource_base_url <-
  "https://zenodo.org/record/7074291/files"


lr_network_file <- file.path(
  nichenet_data_dir,
  "lr_network_human_21122021.rds"
)


ligand_target_file <- file.path(
  nichenet_data_dir,
  "ligand_target_matrix_nsga2r_final.rds"
)


weighted_networks_file <- file.path(
  nichenet_data_dir,
  "weighted_networks_nsga2r_final.rds"
)


ligand_tf_file <- file.path(
  nichenet_data_dir,
  "ligand_tf_matrix_nsga2r_final.rds"
)


# ==========================================================
# 5. Download NicheNet resources if required
# ==========================================================

download_resource(
  
  paste0(
    resource_base_url,
    "/lr_network_human_21122021.rds"
  ),
  
  lr_network_file
  
)


download_resource(
  
  paste0(
    resource_base_url,
    "/ligand_target_matrix_nsga2r_final.rds"
  ),
  
  ligand_target_file
  
)


download_resource(
  
  paste0(
    resource_base_url,
    "/weighted_networks_nsga2r_final.rds"
  ),
  
  weighted_networks_file
  
)


download_resource(
  
  paste0(
    resource_base_url,
    "/ligand_tf_matrix_nsga2r_final.rds"
  ),
  
  ligand_tf_file
  
)


# ==========================================================
# 6. Load NicheNet resources
# ==========================================================

log_message(
  "Loading NicheNet prior models..."
)


lr_network <- readRDS(
  lr_network_file
)


ligand_target_matrix <- readRDS(
  ligand_target_file
)


weighted_networks <- readRDS(
  weighted_networks_file
)


ligand_tf_matrix <- readRDS(
  ligand_tf_file
)


checkpoint(
  all(
    c(
      "from",
      "to"
    ) %in%
      colnames(
        lr_network
      )
  ),
  "Unexpected ligand-receptor network structure."
)


checkpoint(
  is.matrix(
    ligand_target_matrix
  ),
  "ligand_target_matrix is not a matrix."
)


checkpoint(
  is.list(
    weighted_networks
  ),
  "weighted_networks is not a list."
)


checkpoint(
  is.matrix(
    ligand_tf_matrix
  ),
  "ligand_tf_matrix is not a matrix."
)


lr_network <- lr_network %>%
  
  dplyr::select(
    from,
    to
  ) %>%
  
  dplyr::distinct()


log_message(
  "Ligand-receptor interactions:",
  nrow(lr_network)
)


log_message(
  "Ligands in ligand-target model:",
  ncol(ligand_target_matrix)
)


log_message(
  "Targets in ligand-target model:",
  nrow(ligand_target_matrix)
)


# ==========================================================
# 7. Load final GSE183904 Seurat object
# ==========================================================

seurat_file <- file.path(
  objects_dir,
  "seurat_final_GSE183904.rds"
)


checkpoint(
  file.exists(seurat_file),
  paste0(
    "Final Seurat object not found:\n",
    seurat_file,
    "\n\nRun 04_emt_scoring.R first."
  )
)


log_message(
  "Loading final GSE183904 Seurat object..."
)


seurat_obj <- readRDS(
  seurat_file
)


checkpoint(
  inherits(
    seurat_obj,
    "Seurat"
  ),
  "Input object is not a valid Seurat object."
)


checkpoint(
  "cell_type_final" %in%
    colnames(
      seurat_obj@meta.data
    ),
  "cell_type_final metadata is missing."
)


# ==========================================================
# 8. Validate sender and receiver populations
# ==========================================================

available_celltypes <- unique(
  as.character(
    seurat_obj$cell_type_final
  )
)


checkpoint(
  sender_celltype %in%
    available_celltypes,
  paste0(
    "Sender population not found: ",
    sender_celltype
  )
)


checkpoint(
  receiver_celltype %in%
    available_celltypes,
  paste0(
    "Receiver population not found: ",
    receiver_celltype
  )
)


cat(
  "\nSender / receiver cell numbers:\n"
)


print(
  table(
    seurat_obj$cell_type_final[
      seurat_obj$cell_type_final %in%
        c(
          sender_celltype,
          receiver_celltype
        )
    ]
  )
)


# ==========================================================
# 9. Construct NicheNet-specific Seurat object
#
# Only sender and receiver cells are required for this
# analysis.
# ==========================================================

nichenet_seurat <- subset(
  
  seurat_obj,
  
  subset =
    cell_type_final %in%
    c(
      sender_celltype,
      receiver_celltype
    )
  
)


checkpoint(
  ncol(
    nichenet_seurat
  ) > 0,
  "NicheNet Seurat subset is empty."
)


DefaultAssay(
  nichenet_seurat
) <- "RNA"


Idents(
  nichenet_seurat
) <- nichenet_seurat$cell_type_final


# ==========================================================
# 10. Identify receiver-expressed genes
#
# Receiver = EMT_low
# ==========================================================

log_message(
  "Identifying genes expressed in EMT-low receiver cells..."
)


expressed_genes_receiver <-
  
  nichenetr::get_expressed_genes(
    
    receiver_celltype,
    
    nichenet_seurat,
    
    pct =
      expression_pct
    
  )


expressed_genes_receiver <-
  standardize_symbols(
    expressed_genes_receiver
  )


checkpoint(
  length(
    expressed_genes_receiver
  ) > 0,
  "No receiver-expressed genes were detected."
)


log_message(
  "Receiver-expressed genes:",
  length(
    expressed_genes_receiver
  )
)


write.csv(
  
  data.frame(
    gene =
      expressed_genes_receiver
  ),
  
  file.path(
    tables_dir,
    "NicheNet_expressed_genes_receiver.csv"
  ),
  
  row.names =
    FALSE
  
)


# ==========================================================
# 11. Identify T/NK-expressed genes
# ==========================================================

log_message(
  "Identifying genes expressed in T/NK sender cells..."
)


expressed_genes_sender <-
  
  nichenetr::get_expressed_genes(
    
    sender_celltype,
    
    nichenet_seurat,
    
    pct =
      expression_pct
    
  )


expressed_genes_sender <-
  standardize_symbols(
    expressed_genes_sender
  )


checkpoint(
  length(
    expressed_genes_sender
  ) > 0,
  "No sender-expressed genes were detected."
)


log_message(
  "T/NK-expressed genes:",
  length(
    expressed_genes_sender
  )
)


write.csv(
  
  data.frame(
    gene =
      expressed_genes_sender
  ),
  
  file.path(
    tables_dir,
    "NicheNet_expressed_genes_sender.csv"
  ),
  
  row.names =
    FALSE
  
)


############################################################
#
#         SENDER-FOCUSED POTENTIAL LIGANDS
#
############################################################


# ==========================================================
# 12. Receiver-expressed receptors
# ==========================================================

all_receptors <- unique(
  lr_network$to
)


expressed_receptors <- intersect(
  all_receptors,
  expressed_genes_receiver
)


checkpoint(
  length(
    expressed_receptors
  ) > 0,
  "No NicheNet receptors are expressed in EMT-low cells."
)


log_message(
  "Receiver-expressed receptors:",
  length(
    expressed_receptors
  )
)


write.csv(
  
  data.frame(
    receptor =
      expressed_receptors
  ),
  
  file.path(
    tables_dir,
    "NicheNet_expressed_receptors.csv"
  ),
  
  row.names =
    FALSE
  
)


# ==========================================================
# 13. Ligands with an expressed receiver receptor
# ==========================================================

receptor_supported_ligands <- lr_network %>%
  
  dplyr::filter(
    to %in%
      expressed_receptors
  ) %>%
  
  dplyr::pull(
    from
  ) %>%
  
  unique()


# ==========================================================
# 14. Require ligand expression in T/NK sender cells
#
# This converts the historical sender-agnostic logic into
# the intended sender-focused T/NK -> EMT-low analysis.
# ==========================================================

potential_ligands <- intersect(
  
  receptor_supported_ligands,
  
  expressed_genes_sender
  
)


# The ligand must also exist in the ligand-target model.
potential_ligands <- intersect(
  
  potential_ligands,
  
  colnames(
    ligand_target_matrix
  )
  
)


checkpoint(
  length(
    potential_ligands
  ) > 0,
  paste0(
    "No sender-focused potential ligands remain after ",
    "requiring T/NK ligand expression and EMT-low ",
    "receptor expression."
  )
)


log_message(
  "Potential T/NK ligands:",
  length(
    potential_ligands
  )
)


potential_ligand_lr <- lr_network %>%
  
  dplyr::filter(
    
    from %in%
      potential_ligands,
    
    to %in%
      expressed_receptors
    
  ) %>%
  
  dplyr::arrange(
    from,
    to
  )


write.csv(
  
  potential_ligand_lr,
  
  file.path(
    tables_dir,
    "NicheNet_potential_TNK_ligands.csv"
  ),
  
  row.names =
    FALSE
  
)


############################################################
#
#        EMT TRAJECTORY GENE SET OF INTEREST
#
############################################################


# ==========================================================
# 15. Load Monocle3 trajectory modules
# ==========================================================

module_file <- file.path(
  tables_dir,
  "EMT_trajectory_gene_modules.csv"
)


module_trend_file <- file.path(
  tables_dir,
  "EMT_trajectory_module_trends.csv"
)


checkpoint(
  file.exists(module_file),
  paste0(
    "Trajectory module file not found:\n",
    module_file,
    "\nRun 05_emt_trajectory_monocle3.R first."
  )
)


checkpoint(
  file.exists(module_trend_file),
  paste0(
    "Trajectory module-trend file not found:\n",
    module_trend_file,
    "\nRun 05_emt_trajectory_monocle3.R first."
  )
)


gene_module_df <- read.csv(
  module_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


module_trends <- read.csv(
  module_trend_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


checkpoint(
  all(
    c(
      "gene",
      "module"
    ) %in%
      colnames(
        gene_module_df
      )
  ),
  "Trajectory gene-module table has unexpected columns."
)


checkpoint(
  all(
    c(
      "module",
      "high_minus_low"
    ) %in%
      colnames(
        module_trends
      )
  ),
  "Trajectory module-trend table has unexpected columns."
)


# ==========================================================
# 16. Standardize module identifiers
# ==========================================================

gene_module_df$module_id <- as.character(
  gene_module_df$module
)


module_trends$module_id <- sub(
  
  "^Module[[:space:]]+",
  
  "",
  
  as.character(
    module_trends$module
  )
  
)


# ==========================================================
# 17. Select EMT-increasing trajectory modules
# ==========================================================

if (
  is.null(
    manual_emt_modules
  )
) {
  
  if (
    "monotonic_increase" %in%
    colnames(
      module_trends
    )
  ) {
    
    monotonic_flag <- as.logical(
      module_trends$monotonic_increase
    )
    
    
    selected_module_table <- module_trends[
      
      !is.na(
        monotonic_flag
      ) &
        
        monotonic_flag &
        
        module_trends$high_minus_low > 0,
      
      ,
      
      drop = FALSE
      
    ]
    
  } else {
    
    selected_module_table <- module_trends[
      
      module_trends$high_minus_low > 0,
      
      ,
      
      drop = FALSE
      
    ]
    
  }
  
  
  checkpoint(
    nrow(
      selected_module_table
    ) > 0,
    paste0(
      "No EMT-increasing trajectory modules were ",
      "identified automatically."
    )
  )
  
  
  selected_module_table <- selected_module_table[
    order(
      -selected_module_table$high_minus_low
    ),
    ,
    drop = FALSE
  ]
  
  
  selected_modules <- unique(
    selected_module_table$module_id
  )
  
} else {
  
  selected_modules <- as.character(
    manual_emt_modules
  )
  
  
  selected_module_table <- module_trends[
    module_trends$module_id %in%
      selected_modules,
    ,
    drop = FALSE
  ]
  
}


checkpoint(
  length(
    selected_modules
  ) > 0,
  "No trajectory modules were selected."
)


log_message(
  "Selected EMT-increasing modules:",
  paste(
    selected_modules,
    collapse = ", "
  )
)


write.csv(
  
  selected_module_table,
  
  file.path(
    tables_dir,
    "NicheNet_selected_trajectory_modules.csv"
  ),
  
  row.names =
    FALSE
  
)


# ==========================================================
# 18. Build trajectory-derived EMT gene set
# ==========================================================

geneset_oi_raw <- unique(
  
  gene_module_df$gene[
    gene_module_df$module_id %in%
      selected_modules
  ]
  
)


geneset_oi_raw <- standardize_symbols(
  geneset_oi_raw
)


log_message(
  "Trajectory genes before NicheNet filtering:",
  length(
    geneset_oi_raw
  )
)


# ==========================================================
# 19. Define NicheNet receiver background
# ==========================================================

background_expressed_genes <- intersect(
  
  expressed_genes_receiver,
  
  rownames(
    ligand_target_matrix
  )
  
)


checkpoint(
  length(
    background_expressed_genes
  ) > 0,
  "No receiver genes overlap the NicheNet ligand-target matrix."
)


# ==========================================================
# 20. Final gene set of interest
#
# Genes must:
#   1. come from EMT-increasing trajectory modules,
#   2. exist in the NicheNet target model,
#   3. belong to the receiver-expression background.
# ==========================================================

geneset_oi <- Reduce(
  
  intersect,
  
  list(
    
    geneset_oi_raw,
    
    rownames(
      ligand_target_matrix
    ),
    
    background_expressed_genes
    
  )
  
)


checkpoint(
  length(
    geneset_oi
  ) >= 10,
  paste0(
    "Only ",
    length(geneset_oi),
    " trajectory genes remain in the NicheNet gene set. ",
    "Inspect module selection and receiver-expression ",
    "criteria before interpreting ligand activity."
  )
)


log_message(
  "NicheNet EMT gene set:",
  length(
    geneset_oi
  )
)


log_message(
  "NicheNet background:",
  length(
    background_expressed_genes
  )
)


write.csv(
  
  data.frame(
    gene =
      geneset_oi
  ),
  
  file.path(
    tables_dir,
    "NicheNet_geneset_of_interest.csv"
  ),
  
  row.names =
    FALSE
  
)


write.csv(
  
  data.frame(
    gene =
      background_expressed_genes
  ),
  
  file.path(
    tables_dir,
    "NicheNet_background_genes.csv"
  ),
  
  row.names =
    FALSE
  
)


############################################################
#
#              LIGAND ACTIVITY ANALYSIS
#
############################################################


# ==========================================================
# 21. Predict ligand activities
# ==========================================================

log_message(
  "Running sender-focused NicheNet ligand activity analysis..."
)


ligand_activity <- predict_ligand_activities(
  
  geneset =
    geneset_oi,
  
  background_expressed_genes =
    background_expressed_genes,
  
  ligand_target_matrix =
    ligand_target_matrix,
  
  potential_ligands =
    potential_ligands
  
)


checkpoint(
  nrow(
    ligand_activity
  ) > 0,
  "NicheNet returned no ligand activity results."
)


checkpoint(
  "aupr_corrected" %in%
    colnames(
      ligand_activity
    ),
  "aupr_corrected is missing from NicheNet output."
)


ligand_activity <- ligand_activity %>%
  
  dplyr::arrange(
    dplyr::desc(
      aupr_corrected
    )
  ) %>%
  
  dplyr::mutate(
    rank =
      dplyr::row_number()
  )


write.csv(
  
  ligand_activity,
  
  file.path(
    tables_dir,
    "NicheNet_ligand_activity_all.csv"
  ),
  
  row.names =
    FALSE
  
)


# ==========================================================
# 22. Select top ligands
# ==========================================================

n_top <- min(
  number_top_ligands,
  nrow(
    ligand_activity
  )
)


top_ligand_activity <- ligand_activity %>%
  
  dplyr::slice_head(
    n =
      n_top
  )


best_upstream_ligands <-
  top_ligand_activity$test_ligand


checkpoint(
  length(
    best_upstream_ligands
  ) > 0,
  "No prioritized ligands were identified."
)


write.csv(
  
  top_ligand_activity,
  
  file.path(
    tables_dir,
    "NicheNet_top30_ligands.csv"
  ),
  
  row.names =
    FALSE
  
)


cat(
  "\nTop NicheNet ligands:\n"
)


print(
  top_ligand_activity[
    ,
    intersect(
      c(
        "rank",
        "test_ligand",
        "aupr_corrected",
        "aupr",
        "auroc",
        "pearson"
      ),
      colnames(
        top_ligand_activity
      )
    ),
    drop = FALSE
  ]
)


# ==========================================================
# 23. Ligand activity distribution
# ==========================================================

activity_cutoff <- min(
  top_ligand_activity$aupr_corrected,
  na.rm = TRUE
)


p_activity_distribution <- ggplot(
  
  ligand_activity,
  
  aes(
    x =
      aupr_corrected
  )
  
) +
  
  geom_histogram(
    bins = 30
  ) +
  
  geom_vline(
    
    xintercept =
      activity_cutoff,
    
    linetype =
      "dashed"
    
  ) +
  
  labs(
    
    title =
      "NicheNet Ligand Activity Distribution",
    
    subtitle =
      paste0(
        "Dashed line = top ",
        n_top,
        " ligand cutoff"
      ),
    
    x =
      "Corrected AUPR",
    
    y =
      "Number of candidate ligands"
    
  ) +
  
  theme_classic()


ggsave(
  
  filename = file.path(
    plots_dir,
    "NicheNet_01_ligand_activity_distribution.pdf"
  ),
  
  plot =
    p_activity_distribution,
  
  width =
    8,
  
  height =
    6
  
)


# ==========================================================
# 24. Top-ligand activity heatmap
# ==========================================================

vis_ligand_aupr <-
  
  top_ligand_activity %>%
  
  dplyr::select(
    test_ligand,
    aupr_corrected
  ) %>%
  
  tibble::column_to_rownames(
    "test_ligand"
  ) %>%
  
  dplyr::arrange(
    aupr_corrected
  ) %>%
  
  as.matrix()


p_ligand_activity <- nichenetr::make_heatmap_ggplot(
  
  vis_ligand_aupr,
  
  y_name =
    "Prioritized T/NK ligands",
  
  x_name =
    "Ligand activity",
  
  legend_title =
    "Corrected AUPR",
  
  color =
    "darkorange"
  
) +
  
  theme(
    axis.text.x.top =
      element_blank()
  )


ggsave(
  
  filename = file.path(
    plots_dir,
    "NicheNet_02_top_ligand_activity_heatmap.pdf"
  ),
  
  plot =
    p_ligand_activity,
  
  width =
    8,
  
  height =
    8
  
)


############################################################
#
#               LIGAND -> TARGET LINKS
#
############################################################


# ==========================================================
# 25. Infer active ligand-target links
# ==========================================================

log_message(
  "Inferring prioritized ligand-target relationships..."
)


active_ligand_target_links_df <-
  
  best_upstream_ligands %>%
  
  lapply(
    
    nichenetr::get_weighted_ligand_target_links,
    
    geneset =
      geneset_oi,
    
    ligand_target_matrix =
      ligand_target_matrix,
    
    n =
      number_targets_per_ligand
    
  ) %>%
  
  dplyr::bind_rows() %>%
  
  tidyr::drop_na()


checkpoint(
  nrow(
    active_ligand_target_links_df
  ) > 0,
  "No active ligand-target links were identified."
)


active_ligand_target_links_df <-
  
  active_ligand_target_links_df %>%
  
  dplyr::arrange(
    ligand,
    dplyr::desc(
      weight
    )
  )


write.csv(
  
  active_ligand_target_links_df,
  
  file.path(
    tables_dir,
    "NicheNet_active_ligand_target_links.csv"
  ),
  
  row.names =
    FALSE
  
)


# ==========================================================
# 26. Prepare ligand-target visualization matrix
# ==========================================================

active_ligand_target_links <-
  
  nichenetr::prepare_ligand_target_visualization(
    
    ligand_target_df =
      active_ligand_target_links_df,
    
    ligand_target_matrix =
      ligand_target_matrix,
    
    cutoff =
      ligand_target_visualization_cutoff
    
  )


order_ligands <- intersect(
  
  best_upstream_ligands,
  
  colnames(
    active_ligand_target_links
  )
  
)


order_ligands <- rev(
  order_ligands
)


order_targets <- unique(
  active_ligand_target_links_df$target
)


order_targets <- intersect(
  
  order_targets,
  
  rownames(
    active_ligand_target_links
  )
  
)


checkpoint(
  length(
    order_ligands
  ) > 0 &&
    length(
      order_targets
    ) > 0,
  "Ligand-target visualization matrix is empty."
)


vis_ligand_target <- t(
  
  active_ligand_target_links[
    order_targets,
    order_ligands,
    drop = FALSE
  ]
  
)


p_ligand_target <-
  
  nichenetr::make_heatmap_ggplot(
    
    vis_ligand_target,
    
    y_name =
      "Prioritized T/NK ligands",
    
    x_name =
      "EMT trajectory target genes",
    
    color =
      "purple",
    
    legend_title =
      "Regulatory potential"
    
  ) +
  
  scale_fill_gradient2(
    low = "whitesmoke",
    high = "purple"
  )


ggsave(
  
  filename = file.path(
    plots_dir,
    "NicheNet_03_ligand_target_heatmap.pdf"
  ),
  
  plot =
    p_ligand_target,
  
  width =
    14,
  
  height =
    10
  
)


############################################################
#
#              LIGAND -> RECEPTOR LINKS
#
############################################################


# ==========================================================
# 27. Infer prioritized ligand-receptor relationships
# ==========================================================

log_message(
  "Inferring prioritized ligand-receptor relationships..."
)


ligand_receptor_links_df <-
  
  nichenetr::get_weighted_ligand_receptor_links(
    
    best_upstream_ligands,
    
    expressed_receptors,
    
    lr_network,
    
    weighted_networks$lr_sig
    
  )


checkpoint(
  nrow(
    ligand_receptor_links_df
  ) > 0,
  "No prioritized ligand-receptor links were identified."
)


write.csv(
  
  ligand_receptor_links_df,
  
  file.path(
    tables_dir,
    "NicheNet_ligand_receptor_links.csv"
  ),
  
  row.names =
    FALSE
  
)


# ==========================================================
# 28. Ligand-receptor heatmap
# ==========================================================

vis_ligand_receptor_network <-
  
  nichenetr::prepare_ligand_receptor_visualization(
    
    ligand_receptor_links_df,
    
    best_upstream_ligands,
    
    order_hclust =
      "both"
    
  )


p_ligand_receptor <-
  
  nichenetr::make_heatmap_ggplot(
    
    t(
      vis_ligand_receptor_network
    ),
    
    y_name =
      "Prioritized T/NK ligands",
    
    x_name =
      "Receptors expressed in EMT-low cells",
    
    color =
      "mediumvioletred",
    
    legend_title =
      "Interaction potential"
    
  )


ggsave(
  
  filename = file.path(
    plots_dir,
    "NicheNet_04_ligand_receptor_heatmap.pdf"
  ),
  
  plot =
    p_ligand_receptor,
  
  width =
    12,
  
  height =
    9
  
)


############################################################
#
#         LIGAND -> SIGNALING -> TF -> TARGET
#
############################################################


# ==========================================================
# 29. Select ligands compatible with ligand-TF model
# ==========================================================

ligands_for_path <- intersect(
  
  best_upstream_ligands,
  
  colnames(
    ligand_tf_matrix
  )
  
)


ligands_for_path <- head(
  
  ligands_for_path,
  
  number_ligands_for_path
  
)


# ==========================================================
# 30. Select target genes compatible with GRN
# ==========================================================

grn_targets <- unique(
  weighted_networks$gr$to
)


active_targets <- unique(
  active_ligand_target_links_df$target
)


targets_for_path <- intersect(
  
  active_targets,
  
  grn_targets
  
)


targets_for_path <- head(
  
  targets_for_path,
  
  number_targets_for_path
  
)


checkpoint(
  length(
    ligands_for_path
  ) > 0,
  "No prioritized ligands are present in ligand_tf_matrix."
)


checkpoint(
  length(
    targets_for_path
  ) > 0,
  "No active targets are present in NicheNet GRN."
)


log_message(
  "Ligands used for signaling paths:",
  length(
    ligands_for_path
  )
)


log_message(
  "Targets used for signaling paths:",
  length(
    targets_for_path
  )
)


# ==========================================================
# 31. Infer signaling paths
# ==========================================================

log_message(
  "Inferring ligand-to-target signaling paths..."
)


active_signaling_network <-
  
  nichenetr::get_ligand_signaling_path(
    
    ligands_all =
      ligands_for_path,
    
    targets_all =
      targets_for_path,
    
    weighted_networks =
      weighted_networks,
    
    ligand_tf_matrix =
      ligand_tf_matrix,
    
    top_n_regulators =
      top_n_regulators,
    
    minmax_scaling =
      TRUE
    
  )


checkpoint(
  is.list(
    active_signaling_network
  ),
  "Signaling-path inference did not return a valid network."
)


# ==========================================================
# 32. Save signaling-network edges
# ==========================================================

if (
  nrow(
    active_signaling_network$sig
  ) > 0
) {
  
  write.csv(
    
    active_signaling_network$sig,
    
    file.path(
      tables_dir,
      "NicheNet_signaling_edges.csv"
    ),
    
    row.names =
      FALSE
    
  )
  
}


if (
  nrow(
    active_signaling_network$gr
  ) > 0
) {
  
  write.csv(
    
    active_signaling_network$gr,
    
    file.path(
      tables_dir,
      "NicheNet_signaling_GRN_edges.csv"
    ),
    
    row.names =
      FALSE
    
  )
  
}


# ==========================================================
# 33. Prioritize transcription factors
#
# GRN edges have:
#   from = transcription factor
#   to   = predicted target
# ==========================================================

tf_edges <- active_signaling_network$gr


checkpoint(
  nrow(
    tf_edges
  ) > 0,
  "No TF-target edges were returned."
)


checkpoint(
  all(
    c(
      "from",
      "to",
      "weight"
    ) %in%
      colnames(
        tf_edges
      )
  ),
  "Unexpected GRN edge structure."
)


tf_summary <- tf_edges %>%
  
  dplyr::group_by(
    TF = from
  ) %>%
  
  dplyr::summarise(
    
    n_targets =
      dplyr::n_distinct(
        to
      ),
    
    max_weight =
      max(
        weight,
        na.rm = TRUE
      ),
    
    mean_weight =
      mean(
        weight,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
    
  ) %>%
  
  dplyr::arrange(
    
    dplyr::desc(
      max_weight
    ),
    
    dplyr::desc(
      mean_weight
    )
    
  ) %>%
  
  dplyr::mutate(
    rank =
      dplyr::row_number()
  )


top_tfs <- tf_summary %>%
  
  dplyr::slice_head(
    n =
      min(
        number_top_tfs,
        nrow(
          tf_summary
        )
      )
  )


checkpoint(
  nrow(
    top_tfs
  ) > 0,
  "No transcription factors were prioritized."
)


write.csv(
  
  top_tfs,
  
  file.path(
    tables_dir,
    "NicheNet_top_TFs.csv"
  ),
  
  row.names =
    FALSE
  
)


cat(
  "\nTop NicheNet-associated transcription factors:\n"
)


print(
  top_tfs
)


# ==========================================================
# 34. Save reusable NicheNet analysis object
# ==========================================================

nichenet_output <- list(
  
  sender =
    sender_celltype,
  
  receiver =
    receiver_celltype,
  
  expression_pct =
    expression_pct,
  
  selected_modules =
    selected_modules,
  
  geneset_oi =
    geneset_oi,
  
  background_expressed_genes =
    background_expressed_genes,
  
  expressed_genes_sender =
    expressed_genes_sender,
  
  expressed_genes_receiver =
    expressed_genes_receiver,
  
  expressed_receptors =
    expressed_receptors,
  
  potential_ligands =
    potential_ligands,
  
  ligand_activity =
    ligand_activity,
  
  top_ligands =
    best_upstream_ligands,
  
  active_ligand_target_links =
    active_ligand_target_links_df,
  
  ligand_receptor_links =
    ligand_receptor_links_df,
  
  signaling_network =
    active_signaling_network,
  
  top_tfs =
    top_tfs
  
)


nichenet_output_file <- file.path(
  
  objects_dir,
  
  "NicheNet_TNK_EMT_GSE183904.rds"
  
)


saveRDS(
  
  nichenet_output,
  
  nichenet_output_file
  
)


# ==========================================================
# 35. Save session information
# ==========================================================

capture.output(
  
  sessionInfo(),
  
  file = file.path(
    results_dir,
    "sessionInfo_06_NicheNet.txt"
  )
  
)


# ==========================================================
# 36. Completion summary
# ==========================================================

cat(
  
  "\n",
  "=====================================================\n",
  " NICHENET ANALYSIS COMPLETED\n",
  "=====================================================\n",
  "\n",
  "Sender:\n",
  sender_celltype,
  "\n\n",
  "Receiver:\n",
  receiver_celltype,
  "\n\n",
  "Expression threshold:\n",
  expression_pct,
  "\n\n",
  "Potential sender-focused ligands:\n",
  length(
    potential_ligands
  ),
  "\n\n",
  "Selected trajectory modules:\n",
  paste(
    selected_modules,
    collapse = ", "
  ),
  "\n\n",
  "Gene set of interest:\n",
  length(
    geneset_oi
  ),
  " genes\n\n",
  "Prioritized ligands:\n",
  length(
    best_upstream_ligands
  ),
  "\n\n",
  "Prioritized TFs:\n",
  nrow(
    top_tfs
  ),
  "\n\n",
  "NicheNet output object:\n",
  nichenet_output_file,
  "\n",
  "=====================================================\n"
  
)


log_message(
  "06_nichenet_ligand_analysis.R completed."
)