############################################################
# Cell-Type Annotation
# Gastric Cancer scRNA-seq — GSE183904
#
# Input:
#   objects/seurat_integrated_GSE183904.rds
#
# Output:
#   objects/seurat_annotated_GSE183904.rds
#
# Purpose:
#   Annotate major gastric cancer and tumor-microenvironment
#   cell populations using canonical marker genes and
#   manually curated Seurat cluster assignments.
############################################################


# ==========================================================
# 1. Load packages
# ==========================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})


# ==========================================================
# 2. Define project directories
# ==========================================================

base_dir    <- "."
objects_dir <- file.path(base_dir, "objects")
plots_dir   <- file.path(base_dir, "results", "figures")

dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)


# ==========================================================
# 3. Load integrated Seurat object
# ==========================================================

input_file <- file.path(
  objects_dir,
  "seurat_with_metadata_GSE183904.rds"
)

if (!file.exists(input_file)) {
  stop(
    "Integrated Seurat object not found: ",
    input_file,
    "\nRun 01_download_qc_integration.R and 02_add_sample_metadata.R first."
  )
}

seurat_obj <- readRDS(input_file)

if (!inherits(seurat_obj, "Seurat")) {
  stop("Input object is not a valid Seurat object.")
}

if (!"seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
  stop("seurat_clusters not found in metadata.")
}


# ==========================================================
# 4. Define canonical marker sets
# ==========================================================

marker_sets <- list(
  
  Epithelial = c(
    "CDH1",
    "MUC5AC",
    "TFF1",
    "LIPF",
    "PGA3"
  ),
  
  TNK = c(
    "CD2",
    "CD3D",
    "CD3E",
    "NCAM1",
    "NKG7",
    "CD4",
    "CD8A",
    "CD8B",
    "PRF1",
    "IFNG",
    "GZMB"
  ),
  
  B_cell = c(
    "CD79A",
    "MS4A1",
    "CD79B",
    "IGHM",
    "CD37"
  ),
  
  Plasma_cell = c(
    "JCHAIN",
    "MZB1",
    "XBP1",
    "IGHG1",
    "IGHA1"
  ),
  
  Macrophages = c(
    "CD68",
    "CD163",
    "C1QA",
    "C1QB",
    "APOE",
    "MRC1"
  ),
  
  DC = c(
    "FCER1A",
    "CD74",
    "HLA-DRA",
    "HLA-DPA1",
    "CLEC10A"
  ),
  
  Mast_cell = c(
    "TPSB2",
    "TPSAB1",
    "CPA3",
    "KIT",
    "MS4A2"
  ),
  
  CAF = c(
    "FAP",
    "COL1A2",
    "THY1",
    "ACTA2",
    "TAGLN",
    "BGN"
  ),
  
  Endothelial = c(
    "VWF",
    "PECAM1",
    "KDR",
    "PLVAP",
    "FLT1"
  ),
  
  Pericytes = c(
    "RGS5",
    "PDGFRB",
    "MCAM",
    "CSPG4",
    "ACTA2"
  )
)


# ==========================================================
# 5. Keep only markers present in the dataset
# ==========================================================

marker_sets <- lapply(
  marker_sets,
  function(x) intersect(x, rownames(seurat_obj))
)

marker_sets <- marker_sets[
  lengths(marker_sets) > 0
]


# ==========================================================
# 6. Plot Seurat clusters before annotation
# ==========================================================

p_clusters <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  ggtitle("Seurat clusters")

ggsave(
  filename = file.path(
    plots_dir,
    "01_umap_seurat_clusters.pdf"
  ),
  plot = p_clusters,
  width = 10,
  height = 8
)


# ==========================================================
# 7. Generate marker DotPlots
# ==========================================================

pdf(
  file.path(
    plots_dir,
    "02_celltype_marker_dotplots.pdf"
  ),
  width = 12,
  height = 8
)

for (cell_type_name in names(marker_sets)) {
  
  print(
    DotPlot(
      seurat_obj,
      features = marker_sets[[cell_type_name]]
    ) +
      RotatedAxis() +
      ggtitle(cell_type_name)
  )
}

dev.off()


# ==========================================================
# 8. Generate marker FeaturePlots
# ==========================================================

pdf(
  file.path(
    plots_dir,
    "03_celltype_marker_featureplots.pdf"
  ),
  width = 12,
  height = 10
)

for (cell_type_name in names(marker_sets)) {
  
  genes_to_plot <- head(
    marker_sets[[cell_type_name]],
    9
  )
  
  print(
    FeaturePlot(
      seurat_obj,
      features = genes_to_plot,
      reduction = "umap",
      raster = TRUE
    )
  )
}

dev.off()


# ==========================================================
# 9. Manual cluster-to-cell-type annotation
# ==========================================================

cluster_map <- list(
  
  CAF = c(
    8, 11, 22
  ),
  
  Pericytes = c(
    13, 24
  ),
  
  Endothelial = c(
    7, 32
  ),
  
  B_cell = c(
    10
  ),
  
  Plasma_cell = c(
    4, 6, 16, 17, 23
  ),
  
  `T/NK` = c(
    0, 2, 3, 20, 26, 28, 31
  ),
  
  Epithelial = c(
    1, 9, 15, 19, 25, 27, 29, 33, 35
  ),
  
  Macrophages = c(
    5, 12
  ),
  
  DC = c(
    21
  ),
  
  Mast_cell = c(
    34
  )
)


# ==========================================================
# 10. Apply annotations
# ==========================================================

seurat_obj$cell_type <- "Not Known"

for (cell_type_name in names(cluster_map)) {
  
  clusters <- as.character(
    cluster_map[[cell_type_name]]
  )
  
  cells_to_annotate <-
    as.character(seurat_obj$seurat_clusters) %in% clusters
  
  seurat_obj$cell_type[cells_to_annotate] <-
    cell_type_name
}


# ==========================================================
# 11. Annotation summary
# ==========================================================

cat("\nCell-type counts:\n")

print(
  sort(
    table(seurat_obj$cell_type),
    decreasing = TRUE
  )
)


cluster_annotation_table <- data.frame(
  cluster = levels(seurat_obj$seurat_clusters)
)

cluster_annotation_table$cell_type <- sapply(
  cluster_annotation_table$cluster,
  function(cluster_id) {
    
    assigned_types <- unique(
      seurat_obj$cell_type[
        as.character(seurat_obj$seurat_clusters) ==
          cluster_id
      ]
    )
    
    paste(
      assigned_types,
      collapse = "; "
    )
  }
)

write.csv(
  cluster_annotation_table,
  file.path(
    objects_dir,
    "cluster_celltype_annotations.csv"
  ),
  row.names = FALSE
)


# ==========================================================
# 12. Plot final annotated UMAP
# ==========================================================

p_celltypes <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  ggtitle("GSE183904 — Major Cell Types")

ggsave(
  filename = file.path(
    plots_dir,
    "04_umap_cell_types.pdf"
  ),
  plot = p_celltypes,
  width = 10,
  height = 8
)


# ==========================================================
# 13. Save annotated object
# ==========================================================

output_file <- file.path(
  objects_dir,
  "seurat_annotated_GSE183904.rds"
)

saveRDS(
  seurat_obj,
  output_file
)

cat(
  "\nAnnotation completed successfully.\n",
  "Saved object:\n",
  output_file,
  "\n"
)