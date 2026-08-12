############################################################
# TCGA-STAD Validation of EMT-Associated lncRNAs
#
# Project:
#   Gastric Cancer EMT–TME Analysis
#
# External cohort:
#   TCGA-STAD
#
# Input from scRNA-seq pipeline:
#   results/tables/EMT_lncRNA_candidates_for_TCGA.csv
#
# TCGA inputs downloaded automatically:
#   STAR - Counts, Primary Tumor
#   Clinical metadata
#
# Primary endpoint:
#   Overall survival (OS)
#
# Workflow:
#   01. Read candidate EMT-associated lncRNAs
#   02. Download TCGA-STAD primary-tumor STAR counts
#   03. Obtain TCGA-STAD clinical metadata
#   04. Construct overall-survival endpoint
#   05. Resolve duplicate tumor samples per patient
#   06. Transform raw counts using DESeq2 VST
#   07. Match candidate lncRNAs by gene symbol
#   08. Standardize expression (z score)
#   09. Run univariate Cox models
#   10. BH-adjust survival p-values
#   11. Generate per-gene Kaplan-Meier curves
#   12. Optionally fit an exploratory multigene Cox model
#
# IMPORTANT:
#   The multigene risk model is trained and evaluated in the
#   same TCGA-STAD cohort. It is therefore exploratory and
#   should NOT be described as an externally validated
#   prognostic signature.
############################################################


# ==========================================================
# 0. Required packages
# ==========================================================

required_packages <- c(
  "TCGAbiolinks",
  "SummarizedExperiment",
  "DESeq2",
  "survival",
  "survminer",
  "dplyr",
  "tibble",
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
  
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(DESeq2)
  library(survival)
  library(survminer)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  
})


# ==========================================================
# 1. Reproducibility and analysis parameters
# ==========================================================

set.seed(123)


tcga_project <- "TCGA-STAD"


survival_endpoint <- "OS"


# Minimum patients required for an individual Cox model
minimum_survival_samples <- 30


# BH-FDR threshold
cox_fdr_cutoff <- 0.05


# Median expression is used only for Kaplan-Meier
# visualization.
km_split_method <- "median"


# Maximum number of significant genes used in the
# exploratory multigene Cox model.
maximum_multigene_genes <- 5


# ==========================================================
# 2. Project directories
# ==========================================================

base_dir <- "."


data_dir <- file.path(
  base_dir,
  "data"
)


tcga_dir <- file.path(
  data_dir,
  "TCGA_STAD"
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


tcga_plots_dir <- file.path(
  plots_dir,
  "TCGA_STAD"
)


tcga_tables_dir <- file.path(
  tables_dir,
  "TCGA_STAD"
)


dir.create(
  tcga_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  objects_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  tcga_tables_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  tcga_plots_dir,
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


numify <- function(x) {
  
  suppressWarnings(
    as.numeric(
      as.character(
        x
      )
    )
  )
  
}


# ----------------------------------------------------------
# Safe filename conversion
# ----------------------------------------------------------

safe_filename <- function(x) {
  
  gsub(
    "[^A-Za-z0-9_.-]",
    "_",
    x
  )
  
}


# ==========================================================
# 4. Read lncRNA candidates from script 07
# ==========================================================

candidate_file <- file.path(
  tables_dir,
  "EMT_lncRNA_candidates_for_TCGA.csv"
)


checkpoint(
  file.exists(candidate_file),
  paste0(
    "Candidate lncRNA file not found:\n",
    candidate_file,
    "\n\nRun 07_emt_lncrna_analysis.R first."
  )
)


candidate_table <- read.csv(
  candidate_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


checkpoint(
  "gene" %in%
    colnames(
      candidate_table
    ),
  "Candidate file must contain a column named 'gene'."
)


candidate_lncRNAs <- unique(
  trimws(
    as.character(
      candidate_table$gene
    )
  )
)


candidate_lncRNAs <- candidate_lncRNAs[
  !is.na(candidate_lncRNAs) &
    candidate_lncRNAs != ""
]


checkpoint(
  length(candidate_lncRNAs) > 0,
  "No candidate lncRNAs were found."
)


log_message(
  "Candidate EMT-associated lncRNAs:",
  length(candidate_lncRNAs)
)


# ==========================================================
# 5. TCGA expression object cache
# ==========================================================

expression_se_file <- file.path(
  tcga_dir,
  "TCGA_STAD_primary_tumor_STAR_counts.rds"
)


# ==========================================================
# 6. Download TCGA-STAD primary-tumor RNA-seq
# ==========================================================

if (
  !file.exists(
    expression_se_file
  )
) {
  
  log_message(
    "Querying TCGA-STAD primary-tumor STAR counts..."
  )
  
  
  query_exp <- GDCquery(
    
    project =
      tcga_project,
    
    data.category =
      "Transcriptome Profiling",
    
    data.type =
      "Gene Expression Quantification",
    
    workflow.type =
      "STAR - Counts",
    
    sample.type =
      "Primary Tumor"
    
  )
  
  
  log_message(
    "Downloading TCGA-STAD RNA-seq data..."
  )
  
  
  GDCdownload(
    query_exp
  )
  
  
  log_message(
    "Preparing TCGA-STAD expression object..."
  )
  
  
  exp_se <- GDCprepare(
    query_exp
  )
  
  
  saveRDS(
    exp_se,
    expression_se_file
  )
  
  
} else {
  
  log_message(
    "Using cached TCGA-STAD expression object."
  )
  
  
  exp_se <- readRDS(
    expression_se_file
  )
  
}


checkpoint(
  inherits(
    exp_se,
    "SummarizedExperiment"
  ) ||
    inherits(
      exp_se,
      "RangedSummarizedExperiment"
    ),
  "TCGA expression object is not a SummarizedExperiment."
)


# ==========================================================
# 7. Inspect available expression assays
# ==========================================================

available_assays <- SummarizedExperiment::assayNames(
  exp_se
)


cat(
  "\nAvailable TCGA expression assays:\n"
)


print(
  available_assays
)


checkpoint(
  length(
    available_assays
  ) > 0,
  "No expression assays were found."
)


# ----------------------------------------------------------
# STAR-count GDC objects can contain several assays.
# Prefer the raw unstranded count assay.
# ----------------------------------------------------------

preferred_count_assays <- c(
  "unstranded",
  "counts"
)


count_assay <- intersect(
  preferred_count_assays,
  available_assays
)[1]


if (
  is.na(
    count_assay
  )
) {
  
  warning(
    paste0(
      "Could not find an assay explicitly named ",
      "'unstranded' or 'counts'. ",
      "Using the first available assay: ",
      available_assays[1]
    ),
    call. = FALSE
  )
  
  
  count_assay <- available_assays[1]
  
}


log_message(
  "Raw-count assay selected:",
  count_assay
)


raw_counts <- SummarizedExperiment::assay(
  exp_se,
  count_assay
)


checkpoint(
  nrow(raw_counts) > 10000 &&
    ncol(raw_counts) > 50,
  "Unexpected TCGA-STAD expression dimensions."
)


log_message(
  "TCGA raw matrix:",
  nrow(raw_counts),
  "genes x",
  ncol(raw_counts),
  "samples"
)


# ==========================================================
# 8. Obtain gene annotations
# ==========================================================

gene_info <- as.data.frame(
  SummarizedExperiment::rowData(
    exp_se
  )
)


checkpoint(
  nrow(gene_info) ==
    nrow(raw_counts),
  "rowData and expression matrix have different row counts."
)


# ----------------------------------------------------------
# Identify gene-symbol column
# ----------------------------------------------------------

symbol_candidates <- c(
  "gene_name",
  "external_gene_name",
  "gene_symbol",
  "symbol"
)


symbol_column <- intersect(
  symbol_candidates,
  colnames(
    gene_info
  )
)[1]


checkpoint(
  !is.na(
    symbol_column
  ),
  paste0(
    "No gene-symbol column was found in TCGA rowData.\n",
    "Available columns:\n",
    paste(
      colnames(gene_info),
      collapse = ", "
    )
  )
)


gene_symbols <- as.character(
  gene_info[
    [symbol_column]
  ]
)


# ==========================================================
# 9. Create full TCGA sample annotation
# ==========================================================

tcga_barcodes <- colnames(
  raw_counts
)


checkpoint(
  all(
    nchar(
      tcga_barcodes
    ) >= 12
  ),
  "Some TCGA barcodes are shorter than 12 characters."
)


patient_ids <- substr(
  tcga_barcodes,
  1,
  12
)


sample_annotation <- data.frame(
  
  full_barcode =
    tcga_barcodes,
  
  patient_id =
    patient_ids,
  
  library_size =
    colSums(
      raw_counts
    ),
  
  stringsAsFactors =
    FALSE
  
)


# ==========================================================
# 10. Resolve multiple primary-tumor samples per patient
#
# If a patient has >1 primary-tumor RNA-seq sample,
# retain the sample with the largest raw library size.
# ==========================================================

sample_annotation <- sample_annotation %>%
  
  dplyr::arrange(
    patient_id,
    dplyr::desc(
      library_size
    )
  )


selected_samples <- sample_annotation %>%
  
  dplyr::group_by(
    patient_id
  ) %>%
  
  dplyr::slice_head(
    n = 1
  ) %>%
  
  dplyr::ungroup()


duplicated_patient_count <-
  
  nrow(
    sample_annotation
  ) -
  nrow(
    selected_samples
  )


log_message(
  "Primary-tumor RNA-seq samples:",
  nrow(sample_annotation)
)


log_message(
  "Unique patients retained:",
  nrow(selected_samples)
)


log_message(
  "Duplicate patient samples removed:",
  duplicated_patient_count
)


write.csv(
  
  selected_samples,
  
  file.path(
    tcga_tables_dir,
    "TCGA_STAD_selected_primary_tumor_samples.csv"
  ),
  
  row.names =
    FALSE
  
)


# ==========================================================
# 11. Restrict count matrix to one sample per patient
# ==========================================================

raw_counts <- raw_counts[
  ,
  selected_samples$full_barcode,
  drop = FALSE
]


colnames(raw_counts) <- selected_samples$patient_id


checkpoint(
  !any(
    duplicated(
      colnames(
        raw_counts
      )
    )
  ),
  "Duplicated patient IDs remain after sample selection."
)


# ==========================================================
# 12. Download TCGA-STAD clinical metadata
# ==========================================================

clinical_file <- file.path(
  tcga_dir,
  "TCGA_STAD_clinical.rds"
)


if (
  !file.exists(
    clinical_file
  )
) {
  
  log_message(
    "Downloading TCGA-STAD clinical data..."
  )
  
  
  stad_clin <- GDCquery_clinic(
    
    project =
      tcga_project,
    
    type =
      "clinical"
    
  )
  
  
  saveRDS(
    stad_clin,
    clinical_file
  )
  
  
} else {
  
  log_message(
    "Using cached TCGA-STAD clinical data."
  )
  
  
  stad_clin <- readRDS(
    clinical_file
  )
  
}


checkpoint(
  nrow(
    stad_clin
  ) > 50,
  "Unexpectedly small TCGA-STAD clinical cohort."
)


# ==========================================================
# 13. Identify patient barcode column
# ==========================================================

barcode_candidates <- c(
  "bcr_patient_barcode",
  "submitter_id",
  "case_submitter_id"
)


clinical_barcode_column <- intersect(
  barcode_candidates,
  colnames(
    stad_clin
  )
)[1]


checkpoint(
  !is.na(
    clinical_barcode_column
  ),
  paste0(
    "Could not identify the TCGA patient barcode column.\n",
    "Available clinical columns include:\n",
    paste(
      head(
        colnames(
          stad_clin
        ),
        50
      ),
      collapse = ", "
    )
  )
)


stad_clin$patient_id <- substr(
  
  as.character(
    stad_clin[
      [clinical_barcode_column]
    ]
  ),
  
  1,
  12
  
)


# ==========================================================
# 14. Construct Overall Survival endpoint
# ==========================================================

required_survival_columns <- c(
  "days_to_death",
  "days_to_last_follow_up"
)


missing_survival_columns <- setdiff(
  required_survival_columns,
  colnames(
    stad_clin
  )
)


checkpoint(
  length(
    missing_survival_columns
  ) == 0,
  paste0(
    "Clinical survival columns missing: ",
    paste(
      missing_survival_columns,
      collapse = ", "
    )
  )
)


stad_clin$days_to_death <- numify(
  stad_clin$days_to_death
)


stad_clin$days_to_last_follow_up <- numify(
  stad_clin$days_to_last_follow_up
)


# ----------------------------------------------------------
# Death event
#
# Prefer vital_status when available.
# Otherwise fall back to presence of days_to_death.
# ----------------------------------------------------------

if (
  "vital_status" %in%
  colnames(
    stad_clin
  )
) {
  
  vital_status_clean <- tolower(
    trimws(
      as.character(
        stad_clin$vital_status
      )
    )
  )
  
  
  stad_clin$OS_event <- ifelse(
    
    vital_status_clean %in%
      c(
        "dead",
        "deceased"
      ),
    
    1,
    
    ifelse(
      
      vital_status_clean %in%
        c(
          "alive",
          "living"
        ),
      
      0,
      
      NA_real_
      
    )
    
  )
  
  
} else {
  
  stad_clin$OS_event <- ifelse(
    
    !is.na(
      stad_clin$days_to_death
    ) &
      stad_clin$days_to_death > 0,
    
    1,
    
    0
    
  )
  
}


# ----------------------------------------------------------
# Overall-survival time in days
# ----------------------------------------------------------

stad_clin$OS_time <- ifelse(
  
  stad_clin$OS_event == 1 &
    !is.na(
      stad_clin$days_to_death
    ),
  
  stad_clin$days_to_death,
  
  stad_clin$days_to_last_follow_up
  
)


# If vital_status was unavailable/ambiguous but death time is
# recorded, classify as event.
stad_clin$OS_event[
  is.na(
    stad_clin$OS_event
  ) &
    !is.na(
      stad_clin$days_to_death
    ) &
    stad_clin$days_to_death > 0
] <- 1


stad_survival <- stad_clin %>%
  
  dplyr::filter(
    
    !is.na(
      patient_id
    ),
    
    !is.na(
      OS_time
    ),
    
    OS_time > 0,
    
    !is.na(
      OS_event
    )
    
  )


# Remove duplicated clinical patient records if present
stad_survival <- stad_survival[
  !duplicated(
    stad_survival$patient_id
  ),
  ,
  drop = FALSE
]


rownames(
  stad_survival
) <- stad_survival$patient_id


log_message(
  "Patients with usable OS data:",
  nrow(
    stad_survival
  )
)


# ==========================================================
# 15. Match expression and survival data
# ==========================================================

common_patients <- intersect(
  
  colnames(
    raw_counts
  ),
  
  rownames(
    stad_survival
  )
  
)


checkpoint(
  length(
    common_patients
  ) > 50,
  paste0(
    "Only ",
    length(
      common_patients
    ),
    " TCGA patients matched expression and OS data."
  )
)


raw_counts <- raw_counts[
  ,
  common_patients,
  drop = FALSE
]


survival_data <- stad_survival[
  common_patients,
  ,
  drop = FALSE
]


checkpoint(
  identical(
    colnames(
      raw_counts
    ),
    rownames(
      survival_data
    )
  ),
  "Expression and survival samples are not aligned."
)


log_message(
  "Expression-survival matched cohort:",
  length(
    common_patients
  )
)


log_message(
  "OS events:",
  sum(
    survival_data$OS_event == 1
  )
)


# ==========================================================
# 16. Save survival cohort
# ==========================================================

survival_output <- data.frame(
  
  patient_id =
    common_patients,
  
  OS_time_days =
    survival_data$OS_time,
  
  OS_event =
    survival_data$OS_event,
  
  stringsAsFactors =
    FALSE
  
)


write.csv(
  
  survival_output,
  
  file.path(
    tcga_tables_dir,
    "TCGA_STAD_OS_cohort.csv"
  ),
  
  row.names =
    FALSE
  
)


# ==========================================================
# 17. Prepare integer count matrix for DESeq2
# ==========================================================

log_message(
  "Preparing raw counts for DESeq2 VST..."
)


raw_counts <- round(
  raw_counts
)


storage.mode(
  raw_counts
) <- "integer"


checkpoint(
  all(
    raw_counts >= 0
  ),
  "Negative values detected in the selected count assay."
)


# ==========================================================
# 18. Remove genes with no expression
# ==========================================================

keep_genes <- rowSums(
  raw_counts
) > 0


raw_counts_filtered <- raw_counts[
  keep_genes,
  ,
  drop = FALSE
]


gene_symbols_filtered <- gene_symbols[
  keep_genes
]


log_message(
  "Genes retained for VST:",
  nrow(
    raw_counts_filtered
  )
)


# ==========================================================
# 19. DESeq2 variance-stabilizing transformation
# ==========================================================

log_message(
  "Running DESeq2 variance-stabilizing transformation..."
)


deseq_coldata <- data.frame(
  
  patient_id =
    colnames(
      raw_counts_filtered
    ),
  
  row.names =
    colnames(
      raw_counts_filtered
    )
  
)


dds <- DESeqDataSetFromMatrix(
  
  countData =
    raw_counts_filtered,
  
  colData =
    deseq_coldata,
  
  design =
    ~ 1
  
)


# Remove extremely rare genes before dispersion estimation
dds <- dds[
  rowSums(
    counts(
      dds
    )
  ) >= 10,
]


checkpoint(
  nrow(
    dds
  ) > 10000,
  "Too few genes remain after DESeq2 count filtering."
)


vsd <- varianceStabilizingTransformation(
  
  dds,
  
  blind =
    TRUE
  
)


vst_matrix <- assay(
  vsd
)


log_message(
  "VST matrix:",
  nrow(vst_matrix),
  "genes x",
  ncol(vst_matrix),
  "patients"
)


# ==========================================================
# 20. Align gene annotation after DESeq2 filtering
# ==========================================================

original_gene_ids <- rownames(
  raw_counts_filtered
)


gene_symbol_lookup <- setNames(
  
  gene_symbols_filtered,
  
  original_gene_ids
  
)


vst_gene_symbols <- gene_symbol_lookup[
  rownames(
    vst_matrix
  )
]


checkpoint(
  length(
    vst_gene_symbols
  ) ==
    nrow(
      vst_matrix
    ),
  "Gene-symbol annotation did not align with VST matrix."
)


# ==========================================================
# 21. Map scRNA-seq candidate lncRNAs to TCGA genes
# ==========================================================

candidate_mapping_list <- lapply(
  
  candidate_lncRNAs,
  
  function(
    candidate_gene
  ) {
    
    matching_rows <- which(
      vst_gene_symbols ==
        candidate_gene
    )
    
    
    if (
      length(
        matching_rows
      ) == 0
    ) {
      
      return(
        
        data.frame(
          
          gene =
            candidate_gene,
          
          ensembl_id =
            NA_character_,
          
          found_in_TCGA =
            FALSE,
          
          n_matching_rows =
            0,
          
          stringsAsFactors =
            FALSE
          
        )
        
      )
      
    }
    
    
    # ------------------------------------------------------
    # In the unlikely event that multiple Ensembl rows map
    # to the same symbol, retain the row with the highest
    # mean raw expression.
    # ------------------------------------------------------
    
    if (
      length(
        matching_rows
      ) > 1
    ) {
      
      row_ids <- rownames(
        vst_matrix
      )[
        matching_rows
      ]
      
      
      row_means <- rowMeans(
        
        raw_counts_filtered[
          row_ids,
          ,
          drop = FALSE
        ]
        
      )
      
      
      chosen_index <- matching_rows[
        which.max(
          row_means
        )
      ]
      
      
    } else {
      
      chosen_index <- matching_rows
      
    }
    
    
    data.frame(
      
      gene =
        candidate_gene,
      
      ensembl_id =
        rownames(
          vst_matrix
        )[
          chosen_index
        ],
      
      found_in_TCGA =
        TRUE,
      
      n_matching_rows =
        length(
          matching_rows
        ),
      
      stringsAsFactors =
        FALSE
      
    )
    
  }
  
)


candidate_mapping <- dplyr::bind_rows(
  candidate_mapping_list
)


write.csv(
  
  candidate_mapping,
  
  file.path(
    tcga_tables_dir,
    "TCGA_STAD_candidate_lncRNA_mapping.csv"
  ),
  
  row.names =
    FALSE
  
)


candidate_mapping_found <- candidate_mapping %>%
  
  dplyr::filter(
    found_in_TCGA
  )


checkpoint(
  nrow(
    candidate_mapping_found
  ) > 0,
  "None of the scRNA-seq lncRNA candidates were found in TCGA-STAD."
)


log_message(
  "Candidate lncRNAs requested:",
  length(
    candidate_lncRNAs
  )
)


log_message(
  "Candidate lncRNAs found in TCGA:",
  nrow(
    candidate_mapping_found
  )
)


# ==========================================================
# 22. Construct candidate lncRNA VST expression matrix
# ==========================================================

candidate_expression <- vst_matrix[
  candidate_mapping_found$ensembl_id,
  ,
  drop = FALSE
]


rownames(
  candidate_expression
) <- candidate_mapping_found$gene


checkpoint(
  !any(
    duplicated(
      rownames(
        candidate_expression
      )
    )
  ),
  "Duplicated candidate symbols remain after TCGA mapping."
)


# ==========================================================
# 23. Remove zero-variance candidate genes
# ==========================================================

candidate_sd <- apply(
  
  candidate_expression,
  
  1,
  
  sd,
  
  na.rm =
    TRUE
  
)


candidate_expression <- candidate_expression[
  is.finite(
    candidate_sd
  ) &
    candidate_sd > 0,
  ,
  drop = FALSE
]


checkpoint(
  nrow(
    candidate_expression
  ) > 0,
  "All candidate lncRNAs have zero variance in TCGA-STAD."
)


# ==========================================================
# 24. Z-standardize expression
#
# Cox hazard ratios will therefore represent the effect of
# a one-standard-deviation increase in VST expression.
# ==========================================================

candidate_expression_z <- t(
  
  scale(
    
    t(
      candidate_expression
    )
    
  )
  
)


checkpoint(
  all(
    is.finite(
      candidate_expression_z
    )
  ),
  "Non-finite standardized expression values detected."
)


saveRDS(
  
  candidate_expression,
  
  file.path(
    objects_dir,
    "TCGA_STAD_candidate_lncRNA_VST_expression.rds"
  )
  
)


saveRDS(
  
  candidate_expression_z,
  
  file.path(
    objects_dir,
    "TCGA_STAD_candidate_lncRNA_z_expression.rds"
  )
  
)


############################################################
#
#                UNIVARIATE COX ANALYSIS
#
############################################################


# ==========================================================
# 25. Function for individual Cox models
# ==========================================================

run_univariate_cox <- function(
    gene,
    expression_matrix,
    clinical_data
) {
  
  expression_vector <- as.numeric(
    
    expression_matrix[
      gene,
      ,
      drop = TRUE
    ]
    
  )
  
  
  analysis_df <- data.frame(
    
    OS_time =
      clinical_data$OS_time,
    
    OS_event =
      clinical_data$OS_event,
    
    expression_z =
      expression_vector
    
  )
  
  
  analysis_df <- analysis_df[
    
    complete.cases(
      analysis_df
    ),
    
    ,
    
    drop = FALSE
    
  ]
  
  
  if (
    nrow(
      analysis_df
    ) <
    minimum_survival_samples
  ) {
    
    return(
      
      data.frame(
        
        gene =
          gene,
        
        n =
          nrow(
            analysis_df
          ),
        
        events =
          sum(
            analysis_df$OS_event == 1
          ),
        
        beta =
          NA_real_,
        
        HR =
          NA_real_,
        
        CI_low =
          NA_real_,
        
        CI_high =
          NA_real_,
        
        p_value =
          NA_real_,
        
        stringsAsFactors =
          FALSE
        
      )
      
    )
    
  }
  
  
  if (
    length(
      unique(
        analysis_df$OS_event
      )
    ) < 2
  ) {
    
    return(
      
      data.frame(
        
        gene =
          gene,
        
        n =
          nrow(
            analysis_df
          ),
        
        events =
          sum(
            analysis_df$OS_event == 1
          ),
        
        beta =
          NA_real_,
        
        HR =
          NA_real_,
        
        CI_low =
          NA_real_,
        
        CI_high =
          NA_real_,
        
        p_value =
          NA_real_,
        
        stringsAsFactors =
          FALSE
        
      )
      
    )
    
  }
  
  
  fit <- try(
    
    survival::coxph(
      
      survival::Surv(
        OS_time,
        OS_event
      ) ~ expression_z,
      
      data =
        analysis_df
      
    ),
    
    silent =
      TRUE
    
  )
  
  
  if (
    inherits(
      fit,
      "try-error"
    )
  ) {
    
    return(
      
      data.frame(
        
        gene =
          gene,
        
        n =
          nrow(
            analysis_df
          ),
        
        events =
          sum(
            analysis_df$OS_event == 1
          ),
        
        beta =
          NA_real_,
        
        HR =
          NA_real_,
        
        CI_low =
          NA_real_,
        
        CI_high =
          NA_real_,
        
        p_value =
          NA_real_,
        
        stringsAsFactors =
          FALSE
        
      )
      
    )
    
  }
  
  
  fit_summary <- summary(
    fit
  )
  
  
  coefficient <- fit_summary$coefficients[
    1,
    "coef"
  ]
  
  
  hazard_ratio <- fit_summary$coefficients[
    1,
    "exp(coef)"
  ]
  
  
  p_value <- fit_summary$coefficients[
    1,
    "Pr(>|z|)"
  ]
  
  
  confidence_interval <- fit_summary$conf.int[
    1,
    c(
      "lower .95",
      "upper .95"
    )
  ]
  
  
  data.frame(
    
    gene =
      gene,
    
    n =
      fit$n,
    
    events =
      fit$nevent,
    
    beta =
      coefficient,
    
    HR =
      hazard_ratio,
    
    CI_low =
      confidence_interval[1],
    
    CI_high =
      confidence_interval[2],
    
    p_value =
      p_value,
    
    stringsAsFactors =
      FALSE
    
  )
  
}


# ==========================================================
# 26. Run univariate Cox models
# ==========================================================

log_message(
  "Running TCGA-STAD univariate OS Cox models..."
)


cox_results <- lapply(
  
  rownames(
    candidate_expression_z
  ),
  
  run_univariate_cox,
  
  expression_matrix =
    candidate_expression_z,
  
  clinical_data =
    survival_data
  
)


cox_results <- dplyr::bind_rows(
  cox_results
)


# ==========================================================
# 27. BH multiple-testing correction
# ==========================================================

cox_results$p_adj_BH <- p.adjust(
  
  cox_results$p_value,
  
  method =
    "BH"
  
)


cox_results$endpoint <- "Overall survival"


cox_results$expression_scale <-
  
  "HR per 1 SD increase in VST expression"


cox_results <- cox_results %>%
  
  dplyr::arrange(
    p_adj_BH,
    p_value
  )


write.csv(
  
  cox_results,
  
  file.path(
    tcga_tables_dir,
    "TCGA_STAD_univariate_Cox_OS_all_candidates.csv"
  ),
  
  row.names =
    FALSE
  
)


cat(
  "\nUnivariate TCGA-STAD Cox results:\n"
)


print(
  cox_results
)


# ==========================================================
# 28. Significant/prognostic candidate lncRNAs
# ==========================================================

prognostic_lncRNAs <- cox_results %>%
  
  dplyr::filter(
    
    !is.na(
      p_adj_BH
    ),
    
    p_adj_BH <
      cox_fdr_cutoff
    
  )


write.csv(
  
  prognostic_lncRNAs,
  
  file.path(
    tcga_tables_dir,
    "TCGA_STAD_OS_lncRNAs_BH_FDR05.csv"
  ),
  
  row.names =
    FALSE
  
)


log_message(
  "lncRNAs associated with OS at BH-FDR < 0.05:",
  nrow(
    prognostic_lncRNAs
  )
)


############################################################
#
#               KAPLAN-MEIER ANALYSIS
#
############################################################


# ==========================================================
# 29. Function for individual KM plots
#
# Median dichotomization is used for visualization only.
# Continuous-expression Cox analysis above remains the
# primary survival model.
# ==========================================================

generate_km_plot <- function(
    gene,
    expression_matrix,
    clinical_data
) {
  
  expression_vector <- as.numeric(
    
    expression_matrix[
      gene,
      ,
      drop = TRUE
    ]
    
  )
  
  
  median_expression <- median(
    expression_vector,
    na.rm = TRUE
  )
  
  
  expression_group <- ifelse(
    
    expression_vector >
      median_expression,
    
    "High",
    
    "Low"
    
  )
  
  
  km_df <- data.frame(
    
    OS_time =
      clinical_data$OS_time,
    
    OS_event =
      clinical_data$OS_event,
    
    expression =
      expression_vector,
    
    expression_group =
      factor(
        expression_group,
        levels = c(
          "Low",
          "High"
        )
      )
    
  )
  
  
  km_df <- km_df[
    complete.cases(
      km_df
    ),
    ,
    drop = FALSE
  ]
  
  
  fit <- survival::survfit(
    
    survival::Surv(
      OS_time,
      OS_event
    ) ~ expression_group,
    
    data =
      km_df
    
  )
  
  
  plot_object <- survminer::ggsurvplot(
    
    fit,
    
    data =
      km_df,
    
    risk.table =
      TRUE,
    
    pval =
      TRUE,
    
    conf.int =
      FALSE,
    
    xlab =
      "Time (days)",
    
    ylab =
      "Overall survival probability",
    
    title =
      paste0(
        gene,
        " — TCGA-STAD"
      ),
    
    legend.title =
      gene,
    
    legend.labs =
      c(
        "Low",
        "High"
      )
    
  )
  
  
  list(
    
    fit =
      fit,
    
    data =
      km_df,
    
    plot =
      plot_object
    
  )
  
}


# ==========================================================
# 30. Generate KM figures for FDR-significant lncRNAs
# ==========================================================

if (
  nrow(
    prognostic_lncRNAs
  ) > 0
) {
  
  log_message(
    "Generating KM plots for FDR-significant lncRNAs..."
  )
  
  
  for (
    gene
    in prognostic_lncRNAs$gene
  ) {
    
    km_result <- generate_km_plot(
      
      gene =
        gene,
      
      expression_matrix =
        candidate_expression_z,
      
      clinical_data =
        survival_data
      
    )
    
    
    file_stem <- paste0(
      
      "TCGA_STAD_KM_OS_",
      
      safe_filename(
        gene
      )
      
    )
    
    
    pdf(
      
      file.path(
        tcga_plots_dir,
        paste0(
          file_stem,
          ".pdf"
        )
      ),
      
      width =
        7,
      
      height =
        7
      
    )
    
    
    print(
      km_result$plot
    )
    
    
    dev.off()
    
  }
  
  
} else {
  
  log_message(
    "No FDR-significant lncRNAs; individual KM plots skipped."
  )
  
}


# ==========================================================
# 31. Exploratory KM plots for top univariate candidates
#
# These are descriptive and are NOT significance-selected
# validation results.
# ==========================================================

top_exploratory_genes <- cox_results %>%
  
  dplyr::filter(
    !is.na(
      p_value
    )
  ) %>%
  
  dplyr::slice_head(
    n =
      min(
        5,
        sum(
          !is.na(
            cox_results$p_value
          )
        )
      )
  ) %>%
  
  dplyr::pull(
    gene
  )


if (
  length(
    top_exploratory_genes
  ) > 0
) {
  
  for (
    gene
    in top_exploratory_genes
  ) {
    
    km_result <- generate_km_plot(
      
      gene =
        gene,
      
      expression_matrix =
        candidate_expression_z,
      
      clinical_data =
        survival_data
      
    )
    
    
    pdf(
      
      file.path(
        tcga_plots_dir,
        paste0(
          "TCGA_STAD_exploratory_KM_",
          safe_filename(
            gene
          ),
          ".pdf"
        )
      ),
      
      width =
        7,
      
      height =
        7
      
    )
    
    
    print(
      km_result$plot
    )
    
    
    dev.off()
    
  }
  
}


############################################################
#
#       EXPLORATORY MULTIGENE RISK MODEL
#
############################################################


# ==========================================================
# 32. Select OS genes for exploratory multigene model
#
# Only BH-FDR significant candidates are eligible.
# ==========================================================

multigene_candidates <- prognostic_lncRNAs %>%
  
  dplyr::arrange(
    p_adj_BH
  ) %>%
  
  dplyr::pull(
    gene
  ) %>%
  
  unique()


multigene_candidates <- head(
  
  multigene_candidates,
  
  maximum_multigene_genes
  
)


# ==========================================================
# 33. Fit multivariable Cox model when >=2 genes available
# ==========================================================

multigene_model <- NULL

multigene_results <- NULL


if (
  length(
    multigene_candidates
  ) >= 2
) {
  
  log_message(
    "Fitting exploratory multigene Cox model with:",
    paste(
      multigene_candidates,
      collapse = ", "
    )
  )
  
  
  multigene_expression <- t(
    
    candidate_expression_z[
      multigene_candidates,
      ,
      drop = FALSE
    ]
    
  )
  
  
  multigene_expression <- as.data.frame(
    multigene_expression,
    check.names = FALSE
  )
  
  
  multigene_df <- data.frame(
    
    OS_time =
      survival_data$OS_time,
    
    OS_event =
      survival_data$OS_event,
    
    multigene_expression,
    
    check.names =
      FALSE
    
  )
  
  
  multigene_df <- multigene_df[
    complete.cases(
      multigene_df
    ),
    ,
    drop = FALSE
  ]
  
  
  # --------------------------------------------------------
  # Build formula safely even for non-syntactic gene names
  # --------------------------------------------------------
  
  quoted_genes <- paste0(
    "`",
    multigene_candidates,
    "`"
  )
  
  
  multigene_formula <- as.formula(
    
    paste0(
      
      "survival::Surv(OS_time, OS_event) ~ ",
      
      paste(
        quoted_genes,
        collapse = " + "
      )
      
    )
    
  )
  
  
  multigene_model <- survival::coxph(
    
    multigene_formula,
    
    data =
      multigene_df
    
  )
  
  
  multigene_summary <- summary(
    multigene_model
  )
  
  
  capture.output(
    
    multigene_summary,
    
    file = file.path(
      tcga_tables_dir,
      "TCGA_STAD_exploratory_multigene_Cox_summary.txt"
    )
    
  )
  
  
  multigene_results <- data.frame(
    
    gene =
      rownames(
        multigene_summary$coefficients
      ),
    
    beta =
      multigene_summary$coefficients[
        ,
        "coef"
      ],
    
    HR =
      multigene_summary$coefficients[
        ,
        "exp(coef)"
      ],
    
    CI_low =
      multigene_summary$conf.int[
        ,
        "lower .95"
      ],
    
    CI_high =
      multigene_summary$conf.int[
        ,
        "upper .95"
      ],
    
    p_value =
      multigene_summary$coefficients[
        ,
        "Pr(>|z|)"
      ],
    
    stringsAsFactors =
      FALSE
    
  )
  
  
  multigene_results$gene <- gsub(
    
    "^`|`$",
    
    "",
    
    multigene_results$gene
    
  )
  
  
  write.csv(
    
    multigene_results,
    
    file.path(
      tcga_tables_dir,
      "TCGA_STAD_exploratory_multigene_Cox_coefficients.csv"
    ),
    
    row.names =
      FALSE
    
  )
  
  
  # ========================================================
  # 34. Apparent risk score
  # ========================================================
  
  risk_score <- predict(
    
    multigene_model,
    
    type =
      "lp"
    
  )
  
  
  median_risk <- median(
    risk_score,
    na.rm = TRUE
  )
  
  
  risk_group <- ifelse(
    
    risk_score >=
      median_risk,
    
    "High risk",
    
    "Low risk"
    
  )
  
  
  risk_df <- data.frame(
    
    patient_id =
      rownames(
        multigene_df
      ),
    
    OS_time =
      multigene_df$OS_time,
    
    OS_event =
      multigene_df$OS_event,
    
    risk_score =
      risk_score,
    
    risk_group =
      factor(
        risk_group,
        levels = c(
          "Low risk",
          "High risk"
        )
      ),
    
    stringsAsFactors =
      FALSE
    
  )
  
  
  write.csv(
    
    risk_df,
    
    file.path(
      tcga_tables_dir,
      "TCGA_STAD_exploratory_multigene_risk_scores.csv"
    ),
    
    row.names =
      FALSE
    
  )
  
  
  # ========================================================
  # 35. Exploratory multigene KM curve
  # ========================================================
  
  risk_survfit <- survival::survfit(
    
    survival::Surv(
      OS_time,
      OS_event
    ) ~ risk_group,
    
    data =
      risk_df
    
  )
  
  
  risk_km <- survminer::ggsurvplot(
    
    risk_survfit,
    
    data =
      risk_df,
    
    risk.table =
      TRUE,
    
    pval =
      TRUE,
    
    xlab =
      "Time (days)",
    
    ylab =
      "Overall survival probability",
    
    legend.title =
      "Exploratory risk group",
    
    legend.labs =
      c(
        "Low risk",
        "High risk"
      ),
    
    title =
      "TCGA-STAD Exploratory Multi-lncRNA Risk Model"
    
  )
  
  
  pdf(
    
    file.path(
      tcga_plots_dir,
      "TCGA_STAD_exploratory_multigene_KM.pdf"
    ),
    
    width =
      7,
    
    height =
      7
    
  )
  
  
  print(
    risk_km
  )
  
  
  dev.off()
  
  
} else {
  
  log_message(
    paste0(
      "Fewer than two BH-FDR significant lncRNAs. ",
      "Exploratory multigene model skipped."
    )
  )
  
}


############################################################
#
#                    ANALYSIS SUMMARY
#
############################################################


# ==========================================================
# 36. Save analysis summary
# ==========================================================

analysis_summary <- data.frame(
  
  metric = c(
    
    "Candidate lncRNAs from scRNA-seq",
    
    "Candidates found in TCGA",
    
    "TCGA patients with expression and OS",
    
    "OS events",
    
    "Candidate lncRNAs tested by Cox",
    
    "OS-associated lncRNAs BH-FDR < 0.05",
    
    "Genes in exploratory multigene model"
    
  ),
  
  value = c(
    
    length(
      candidate_lncRNAs
    ),
    
    nrow(
      candidate_mapping_found
    ),
    
    length(
      common_patients
    ),
    
    sum(
      survival_data$OS_event == 1
    ),
    
    nrow(
      cox_results
    ),
    
    nrow(
      prognostic_lncRNAs
    ),
    
    length(
      multigene_candidates
    )
    
  ),
  
  stringsAsFactors =
    FALSE
  
)


write.csv(
  
  analysis_summary,
  
  file.path(
    tcga_tables_dir,
    "TCGA_STAD_survival_analysis_summary.csv"
  ),
  
  row.names =
    FALSE
  
)


# ==========================================================
# 37. Save reusable TCGA survival object
# ==========================================================

tcga_survival_object <- list(
  
  project =
    tcga_project,
  
  endpoint =
    survival_endpoint,
  
  candidate_lncRNAs =
    candidate_lncRNAs,
  
  candidate_mapping =
    candidate_mapping,
  
  candidate_expression_VST =
    candidate_expression,
  
  candidate_expression_z =
    candidate_expression_z,
  
  survival_data =
    survival_output,
  
  univariate_cox =
    cox_results,
  
  prognostic_lncRNAs =
    prognostic_lncRNAs,
  
  exploratory_multigene_genes =
    multigene_candidates,
  
  exploratory_multigene_results =
    multigene_results
  
)


tcga_survival_object_file <- file.path(
  
  objects_dir,
  
  "TCGA_STAD_lncRNA_survival_validation.rds"
  
)


saveRDS(
  
  tcga_survival_object,
  
  tcga_survival_object_file
  
)


# ==========================================================
# 38. Save session information
# ==========================================================

capture.output(
  
  sessionInfo(),
  
  file = file.path(
    results_dir,
    "sessionInfo_08_TCGA_STAD_survival.txt"
  )
  
)


# ==========================================================
# 39. Completion message
# ==========================================================

cat(
  
  "\n",
  "=====================================================\n",
  " TCGA-STAD lncRNA SURVIVAL ANALYSIS COMPLETED\n",
  "=====================================================\n",
  "\n",
  "Endpoint:\n",
  "Overall survival (OS)",
  "\n\n",
  "scRNA-seq candidate lncRNAs:\n",
  length(
    candidate_lncRNAs
  ),
  "\n\n",
  "Candidates detected in TCGA-STAD:\n",
  nrow(
    candidate_mapping_found
  ),
  "\n\n",
  "Matched TCGA-STAD patients:\n",
  length(
    common_patients
  ),
  "\n\n",
  "OS events:\n",
  sum(
    survival_data$OS_event == 1
  ),
  "\n\n",
  "BH-FDR significant OS-associated lncRNAs:\n",
  nrow(
    prognostic_lncRNAs
  ),
  "\n\n",
  "Results directory:\n",
  tcga_tables_dir,
  "\n\n",
  "Main survival object:\n",
  tcga_survival_object_file,
  "\n",
  "=====================================================\n"
  
)


log_message(
  "08_tcga_stad_survival.R completed."
)