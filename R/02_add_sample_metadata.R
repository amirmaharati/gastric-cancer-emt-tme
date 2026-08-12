############################################################
# Add Sample and Clinical Metadata to GSE183904
#
# Project:
#   Gastric Cancer EMT–TME Analysis
#
# Input:
#   objects/seurat_integrated_GSE183904.rds
#   data/metadata/Supp.csv
#
# Output:
#   objects/seurat_with_metadata_GSE183904.rds
#
#   results/tables/GSE183904_sample_metadata.csv
#   results/tables/GSE183904_metadata_matching_summary.csv
#   results/figures/UMAP_sample_status_GSE183904.pdf
#   results/figures/UMAP_stage_GSE183904.pdf
#   results/figures/UMAP_lauren_subtype_GSE183904.pdf
#
# Purpose:
#   Add sample-level clinical annotations to each cell
#   using GEO GSM sample identifiers.
############################################################


# ==========================================================
# 0. Required packages
# ==========================================================

required_packages <- c(
  "Seurat",
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
      "Missing required packages:\n",
      paste(
        missing_packages,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
  
}


suppressPackageStartupMessages({
  
  library(Seurat)
  library(ggplot2)
  
})


# ==========================================================
# 1. Reproducibility settings
# ==========================================================

set.seed(123)


# ==========================================================
# 2. Define project directories
# ==========================================================

base_dir <- "."


metadata_dir <- file.path(
  base_dir,
  "data",
  "metadata"
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
  metadata_dir,
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


# ==========================================================
# 4. Define input files
# ==========================================================

seurat_file <- file.path(
  objects_dir,
  "seurat_integrated_GSE183904.rds"
)


metadata_file <- file.path(
  metadata_dir,
  "Supp.csv"
)


checkpoint(
  file.exists(seurat_file),
  paste0(
    "Seurat object not found:\n",
    seurat_file,
    "\nRun 01_download_qc_integration.R first."
  )
)


checkpoint(
  file.exists(metadata_file),
  paste0(
    "Metadata file not found:\n",
    metadata_file
  )
)


# ==========================================================
# 5. Load integrated Seurat object
# ==========================================================

log_message(
  "Loading integrated GSE183904 Seurat object..."
)


seurat_obj <- readRDS(
  seurat_file
)


checkpoint(
  inherits(seurat_obj, "Seurat"),
  "Input object is not a valid Seurat object."
)


checkpoint(
  "sample" %in%
    colnames(
      seurat_obj@meta.data
    ),
  "Column 'sample' is missing from Seurat metadata."
)


log_message(
  "Cells loaded:",
  ncol(seurat_obj)
)


log_message(
  "Samples detected:",
  length(
    unique(
      seurat_obj$sample
    )
  )
)


# ==========================================================
# 6. Read supplementary sample metadata
# ==========================================================

log_message(
  "Reading Supp.csv..."
)


supp <- read.csv(
  metadata_file,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c(
    "",
    "NA",
    "N/A"
  )
)


checkpoint(
  "Sample ID" %in% colnames(supp),
  "'Sample ID' column not found in Supp.csv."
)


checkpoint(
  "Tissue type" %in% colnames(supp),
  "'Tissue type' column not found in Supp.csv."
)


# ==========================================================
# 7. Keep rows corresponding to actual scRNA-seq samples
# ==========================================================

sample_metadata <- supp[
  !is.na(supp[["Sample ID"]]) &
    trimws(
      supp[["Sample ID"]]
    ) != "",
]


rownames(sample_metadata) <- NULL


log_message(
  "Rows containing valid Sample IDs:",
  nrow(sample_metadata)
)


checkpoint(
  nrow(sample_metadata) > 0,
  "No valid Sample IDs found in Supp.csv."
)


# ==========================================================
# 8. Validate unique GSM identifiers
# ==========================================================

sample_metadata[["Sample ID"]] <- trimws(
  sample_metadata[["Sample ID"]]
)


duplicated_samples <- sample_metadata[
  duplicated(
    sample_metadata[["Sample ID"]]
  ),
  "Sample ID"
]


checkpoint(
  length(duplicated_samples) == 0,
  paste0(
    "Duplicated Sample IDs detected:\n",
    paste(
      unique(duplicated_samples),
      collapse = ", "
    )
  )
)


log_message(
  "Unique metadata samples:",
  length(
    unique(
      sample_metadata[["Sample ID"]]
    )
  )
)


# ==========================================================
# 9. Standardize metadata column names
# ==========================================================

rename_if_present <- function(
    data,
    old_name,
    new_name
) {
  
  if (old_name %in% colnames(data)) {
    
    colnames(data)[
      colnames(data) == old_name
    ] <- new_name
    
  }
  
  data
  
}


sample_metadata <- rename_if_present(
  sample_metadata,
  "Patient ID",
  "patient_id"
)


sample_metadata <- rename_if_present(
  sample_metadata,
  "Sample number",
  "sample_number"
)


sample_metadata <- rename_if_present(
  sample_metadata,
  "Sample ID",
  "sample"
)


sample_metadata <- rename_if_present(
  sample_metadata,
  "Tissue type",
  "tissue_type"
)


sample_metadata <- rename_if_present(
  sample_metadata,
  "Tumor/Normal",
  "tumor_normal"
)


sample_metadata <- rename_if_present(
  sample_metadata,
  "Age ",
  "age"
)


sample_metadata <- rename_if_present(
  sample_metadata,
  "Gender",
  "gender"
)


sample_metadata <- rename_if_present(
  sample_metadata,
  "Location",
  "location"
)


sample_metadata <- rename_if_present(
  sample_metadata,
  "Stage*",
  "stage"
)


sample_metadata <- rename_if_present(
  sample_metadata,
  "Laurens",
  "lauren_subtype"
)


sample_metadata <- rename_if_present(
  sample_metadata,
  "IM Status",
  "im_status"
)


sample_metadata <- rename_if_present(
  sample_metadata,
  "MMR",
  "mmr"
)


sample_metadata <- rename_if_present(
  sample_metadata,
  "MMR-Low/High",
  "mmr_status"
)


sample_metadata <- rename_if_present(
  sample_metadata,
  "EBV",
  "ebv"
)


sample_metadata <- rename_if_present(
  sample_metadata,
  "TCGA",
  "tcga_subtype"
)


sample_metadata <- rename_if_present(
  sample_metadata,
  "Primary Tumor",
  "primary_tumor"
)


sample_metadata <- rename_if_present(
  sample_metadata,
  "Primary Normal",
  "primary_normal"
)


sample_metadata <- rename_if_present(
  sample_metadata,
  "Peritoneal Tumor",
  "peritoneal_tumor"
)


sample_metadata <- rename_if_present(
  sample_metadata,
  "Peritoneal Normal",
  "peritoneal_normal"
)


sample_metadata <- rename_if_present(
  sample_metadata,
  "Organoid Tumor",
  "organoid_tumor"
)


sample_metadata <- rename_if_present(
  sample_metadata,
  "Organoid Normal",
  "organoid_normal"
)


# ==========================================================
# 10. Remove irrelevant empty import columns
# ==========================================================

unwanted_columns <- grep(
  "^Unnamed",
  colnames(sample_metadata),
  value = TRUE
)


if (length(unwanted_columns) > 0) {
  
  sample_metadata[
    unwanted_columns
  ] <- NULL
  
}


# ==========================================================
# 11. Create standardized biological sample status
# ==========================================================

sample_metadata$tissue_type <- trimws(
  sample_metadata$tissue_type
)


sample_metadata$sample_status <- NA_character_


sample_metadata$sample_status[
  tolower(
    sample_metadata$tissue_type
  ) == "normal"
] <- "Normal"


sample_metadata$sample_status[
  tolower(
    sample_metadata$tissue_type
  ) == "tumor"
] <- "Tumor"


sample_metadata$sample_status[
  tolower(
    sample_metadata$tissue_type
  ) == "metastasis"
] <- "Metastasis"


checkpoint(
  all(
    !is.na(
      sample_metadata$sample_status
    )
  ),
  paste0(
    "Unrecognized tissue types detected:\n",
    paste(
      unique(
        sample_metadata$tissue_type[
          is.na(
            sample_metadata$sample_status
          )
        ]
      ),
      collapse = ", "
    )
  )
)


sample_metadata$sample_status <- factor(
  sample_metadata$sample_status,
  levels = c(
    "Normal",
    "Tumor",
    "Metastasis"
  )
)


# ==========================================================
# 12. Metadata summary
# ==========================================================

cat(
  "\nSample status distribution:\n"
)


print(
  table(
    sample_metadata$sample_status
  )
)


cat(
  "\nStage distribution:\n"
)


print(
  table(
    sample_metadata$stage,
    useNA = "ifany"
  )
)


cat(
  "\nLauren subtype distribution:\n"
)


print(
  table(
    sample_metadata$lauren_subtype,
    useNA = "ifany"
  )
)


# ==========================================================
# 13. Compare samples between Seurat and Supp.csv
# ==========================================================

seurat_samples <- sort(
  unique(
    as.character(
      seurat_obj$sample
    )
  )
)


metadata_samples <- sort(
  unique(
    as.character(
      sample_metadata$sample
    )
  )
)


missing_in_metadata <- setdiff(
  seurat_samples,
  metadata_samples
)


missing_in_seurat <- setdiff(
  metadata_samples,
  seurat_samples
)


cat(
  "\nSamples in Seurat but absent from metadata:\n"
)


if (length(missing_in_metadata) == 0) {
  
  cat(
    "None\n"
  )
  
} else {
  
  print(
    missing_in_metadata
  )
  
}


cat(
  "\nSamples in metadata but absent from Seurat:\n"
)


if (length(missing_in_seurat) == 0) {
  
  cat(
    "None\n"
  )
  
} else {
  
  print(
    missing_in_seurat
  )
  
}


checkpoint(
  length(missing_in_metadata) == 0,
  paste0(
    "Some Seurat samples have no metadata:\n",
    paste(
      missing_in_metadata,
      collapse = ", "
    )
  )
)


# ==========================================================
# 14. Create metadata lookup table
# ==========================================================

rownames(sample_metadata) <- sample_metadata$sample


metadata_columns_to_add <- c(
  "patient_id",
  "sample_number",
  "tissue_type",
  "tumor_normal",
  "age",
  "gender",
  "location",
  "stage",
  "lauren_subtype",
  "im_status",
  "mmr",
  "mmr_status",
  "ebv",
  "tcga_subtype",
  "primary_tumor",
  "primary_normal",
  "peritoneal_tumor",
  "peritoneal_normal",
  "organoid_tumor",
  "organoid_normal",
  "sample_status"
)


metadata_columns_to_add <- intersect(
  metadata_columns_to_add,
  colnames(sample_metadata)
)


# ==========================================================
# 15. Match metadata to every single cell
# ==========================================================

log_message(
  "Matching sample metadata to individual cells..."
)


cell_sample_ids <- as.character(
  seurat_obj$sample
)


metadata_match <- match(
  cell_sample_ids,
  sample_metadata$sample
)


checkpoint(
  all(
    !is.na(
      metadata_match
    )
  ),
  "At least one cell could not be matched to sample metadata."
)


for (
  metadata_column
  in metadata_columns_to_add
) {
  
  seurat_obj@meta.data[
    [metadata_column]
  ] <- sample_metadata[
    metadata_match,
    metadata_column
  ]
  
}


# ==========================================================
# 16. Validate metadata addition
# ==========================================================

checkpoint(
  "sample_status" %in%
    colnames(
      seurat_obj@meta.data
    ),
  "sample_status was not added successfully."
)


checkpoint(
  all(
    !is.na(
      seurat_obj$sample_status
    )
  ),
  "Some cells have missing sample_status."
)


cat(
  "\nCell counts by sample status:\n"
)


print(
  table(
    seurat_obj$sample_status
  )
)


# ==========================================================
# 17. Save clean sample-level metadata table
# ==========================================================

sample_metadata_output <- file.path(
  tables_dir,
  "GSE183904_sample_metadata.csv"
)


write.csv(
  sample_metadata,
  sample_metadata_output,
  row.names = FALSE
)


log_message(
  "Saved cleaned sample metadata:",
  sample_metadata_output
)


# ==========================================================
# 18. Save matching summary
# ==========================================================

matching_summary <- data.frame(
  
  metric = c(
    "Seurat samples",
    "Metadata samples",
    "Matched samples",
    "Missing metadata samples",
    "Metadata-only samples"
  ),
  
  value = c(
    length(seurat_samples),
    length(metadata_samples),
    length(
      intersect(
        seurat_samples,
        metadata_samples
      )
    ),
    length(missing_in_metadata),
    length(missing_in_seurat)
  )
  
)


matching_summary_file <- file.path(
  tables_dir,
  "GSE183904_metadata_matching_summary.csv"
)


write.csv(
  matching_summary,
  matching_summary_file,
  row.names = FALSE
)


# ==========================================================
# 19. UMAP — Normal / Tumor / Metastasis
# ==========================================================

p_status <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "sample_status",
  raster = FALSE
) +
  ggtitle(
    "GSE183904 — Sample Status"
  )


ggsave(
  filename = file.path(
    plots_dir,
    "UMAP_sample_status_GSE183904.pdf"
  ),
  plot = p_status,
  width = 10,
  height = 8
)


# ==========================================================
# 20. UMAP — Stage
# ==========================================================

if (
  "stage" %in%
  colnames(
    seurat_obj@meta.data
  )
) {
  
  p_stage <- DimPlot(
    seurat_obj,
    reduction = "umap",
    group.by = "stage",
    raster = FALSE
  ) +
    ggtitle(
      "GSE183904 — Clinical Stage"
    )
  
  
  ggsave(
    filename = file.path(
      plots_dir,
      "UMAP_stage_GSE183904.pdf"
    ),
    plot = p_stage,
    width = 10,
    height = 8
  )
  
}


# ==========================================================
# 21. UMAP — Lauren subtype
# ==========================================================

if (
  "lauren_subtype" %in%
  colnames(
    seurat_obj@meta.data
  )
) {
  
  p_lauren <- DimPlot(
    seurat_obj,
    reduction = "umap",
    group.by = "lauren_subtype",
    raster = FALSE
  ) +
    ggtitle(
      "GSE183904 — Lauren Histologic Subtype"
    )
  
  
  ggsave(
    filename = file.path(
      plots_dir,
      "UMAP_lauren_subtype_GSE183904.pdf"
    ),
    plot = p_lauren,
    width = 10,
    height = 8
  )
  
}


# ==========================================================
# 22. Save Seurat object with metadata
# ==========================================================

output_file <- file.path(
  objects_dir,
  "seurat_with_metadata_GSE183904.rds"
)


saveRDS(
  seurat_obj,
  output_file
)


log_message(
  "Saved Seurat object with clinical metadata:",
  output_file
)


# ==========================================================
# 23. Save session information
# ==========================================================

capture.output(
  sessionInfo(),
  file = file.path(
    results_dir,
    "sessionInfo_02_metadata.txt"
  )
)


# ==========================================================
# 24. Completion message
# ==========================================================

cat(
  "\n",
  "=====================================================\n",
  " GSE183904 METADATA INTEGRATION COMPLETED\n",
  "=====================================================\n",
  "\n",
  "Input object:\n",
  seurat_file,
  "\n\n",
  "Metadata source:\n",
  metadata_file,
  "\n\n",
  "Output object:\n",
  output_file,
  "\n",
  "=====================================================\n"
)