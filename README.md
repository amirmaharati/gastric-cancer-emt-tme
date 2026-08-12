# Gastric Cancer EMT–TME

Reproducible single-cell RNA-seq analysis of **epithelial–mesenchymal transition (EMT)** and **tumor-microenvironment (TME) signaling** in gastric cancer using **GSE183904**, followed by lncRNA prioritization and **TCGA-STAD** survival validation.

## Overview

This repository implements a focused reanalysis of GSE183904 to investigate:

- major gastric-cancer and tumor-microenvironment cell populations;
- epithelial EMT states;
- EMT-associated cellular trajectories;
- T/NK-to-epithelial ligand signaling with NicheNet;
- transcription-factor–lncRNA associations;
- EMT-associated lncRNA candidates;
- TCGA-STAD overall-survival associations.

This repository is a **focused reanalysis**, not a complete reproduction of every analysis in the original GSE183904 publication.

## Analysis workflow

```text
GSE183904 scRNA-seq
        ↓
Raw count matrices
        ↓
Quality control
        ↓
Sample-wise DoubletFinder
        ↓
Normalization + PCA
        ↓
Harmony batch correction
        ↓
UMAP + clustering
        ↓
Sample metadata integration
        ↓
Cell-type annotation
        ↓
Epithelial-cell isolation
        ↓
EMT scoring
        ↓
EMT-low / intermediate / high
        ↓
Monocle3 trajectory
        ↓
Trajectory-associated gene modules
        ↓
T/NK → EMT-low NicheNet analysis
        ↓
TF–lncRNA analysis
        ↓
EMT-associated lncRNA prioritization
        ↓
TCGA-STAD overall-survival validation
```

## Repository structure

```text
gastric-cancer-emt-tme/
│
├── README.md
├── LICENSE
├── CITATION.cff
├── .gitignore
│
├── R/
│   ├── 01_download_qc_integration.R
│   ├── 02_add_sample_metadata.R
│   ├── 03_cell_type_annotation.R
│   ├── 04_emt_scoring.R
│   ├── 05_emt_trajectory_monocle3.R
│   ├── 06_nichenet_ligand_analysis.R
│   ├── 07_emt_lncrna_analysis.R
│   ├── 08_tcga_stad_survival.R
│   └── 09_generate_results_summary.R
│
├── data/
│   ├── README.md
│   └── metadata/
│       └── Supp.csv
│
├── results/
│   ├── figures/
│   ├── tables/
│   └── summary/
│
└── objects/                     # generated locally; gitignored
```

Large raw datasets, downloaded reference resources, TCGA files, and generated `.rds` objects are intentionally excluded from version control.

## Data sources

### GSE183904

The primary single-cell dataset is **GEO accession GSE183904**.

The pipeline downloads the raw count archive automatically in:

```text
R/01_download_qc_integration.R
```

Sample-level annotations are read from:

```text
data/metadata/Supp.csv
```

See [`data/README.md`](data/README.md) for full data provenance and local storage details.

### NicheNet

Human NicheNet prior models are downloaded automatically by:

```text
R/06_nichenet_ligand_analysis.R
```

The analysis is implemented as a sender-focused:

```text
T/NK cells → EMT-low epithelial cells
```

signaling analysis.

### GENCODE

Long non-coding RNAs are annotated using **GENCODE human release 49**.

The annotation is obtained automatically by:

```text
R/07_emt_lncrna_analysis.R
```

### TCGA-STAD

Candidate lncRNAs are evaluated in **TCGA Stomach Adenocarcinoma (TCGA-STAD)** using primary-tumor RNA-seq and clinical overall-survival data.

TCGA acquisition and analysis are implemented in:

```text
R/08_tcga_stad_survival.R
```

## Pipeline

### 1. Preprocessing and integration

```text
R/01_download_qc_integration.R
```

Main steps:

- GSE183904 download and extraction;
- Seurat object construction;
- mitochondrial, ribosomal, and hemoglobin QC metrics;
- cell-level filtering;
- sample-wise DoubletFinder;
- normalization;
- 2,000 variable features;
- PCA;
- Harmony integration;
- UMAP;
- graph-based clustering.

Current QC thresholds:

```text
nFeature_RNA > 300
nCount_RNA   > 1,000
nCount_RNA   < 50,000
percent.mt   < 30
percent.ribo < 45
percent.hb   < 1
```

The reconstructed DoubletFinder section currently contains an explicit expected-doublet-rate parameter. This should be updated if the original experiment-specific loading rate is recovered.

### 2. Sample metadata

```text
R/02_add_sample_metadata.R
```

Adds sample-level metadata, including available:

- patient ID;
- Normal / Tumor / Metastasis status;
- stage;
- Lauren subtype;
- MMR information;
- EBV status;
- TCGA subtype.

### 3. Cell-type annotation

```text
R/03_cell_type_annotation.R
```

Major annotated populations include:

- Epithelial;
- T/NK;
- B cells;
- Plasma cells;
- Macrophages;
- Dendritic cells;
- Mast cells;
- Cancer-associated fibroblasts;
- Endothelial cells;
- Pericytes.

### 4. EMT scoring

```text
R/04_emt_scoring.R
```

Epithelial cells are reprocessed independently and clustered at:

```text
resolution = 2.0
```

An EMT composite score is calculated from mesenchymal and epithelial gene programs.

A five-component Gaussian mixture model defines:

```text
EMT_rank1
EMT_rank2
EMT_rank3
EMT_rank4
EMT_rank5
```

which are collapsed into:

```text
EMT_rank1 + EMT_rank2 → EMT_low
EMT_rank3             → EMT_intermediate
EMT_rank4 + EMT_rank5 → EMT_high
```

### 5. EMT trajectory

```text
R/05_emt_trajectory_monocle3.R
```

Monocle3 is used to infer an epithelial trajectory rooted programmatically in **EMT-low cells**.

Trajectory-associated genes are selected using:

```text
q-value < 0.05
Moran's I > 0.1
```

Trajectory-dependent gene modules are identified at:

```text
resolution = 0.001
```

### 6. T/NK-to-epithelial signaling

```text
R/06_nichenet_ligand_analysis.R
```

Sender-focused NicheNet analysis requires:

1. ligand expression in T/NK cells;
2. cognate receptor expression in EMT-low epithelial cells;
3. compatibility with the NicheNet prior model.

Ligands are ranked by corrected AUPR against EMT-associated trajectory genes.

### 7. EMT-associated lncRNAs

```text
R/07_emt_lncrna_analysis.R
```

Candidate lncRNAs are prioritized by integrating:

- GENCODE v49 lncRNA annotation;
- NicheNet-associated transcription factors;
- TF–lncRNA Spearman correlations;
- EMT-high versus EMT-low differential expression.

Historical TF–lncRNA criterion:

```text
|rho| >= 0.5
raw p < 0.05
```

A BH-adjusted sensitivity analysis is also generated.

EMT differential-expression criterion:

```text
adjusted p < 0.05
|log2FC| > 0.25
```

### 8. TCGA-STAD validation

```text
R/08_tcga_stad_survival.R
```

The survival workflow uses:

```text
Primary Tumor STAR counts
        ↓
DESeq2 variance-stabilizing transformation
        ↓
candidate lncRNA expression
        ↓
z-standardization
        ↓
univariate Cox proportional-hazards regression
        ↓
Benjamini–Hochberg correction
```

The primary endpoint is **overall survival (OS)**.

Hazard ratios are interpreted per **1 standard-deviation increase in transformed expression**.

Any multigene risk model produced by the script is exploratory because model fitting and evaluation are performed in the same TCGA-STAD cohort.

### 9. Results summary

```text
R/09_generate_results_summary.R
```

This script reads the standardized outputs from the preceding analyses and generates repository-ready summary tables and a Markdown analysis summary under:

```text
results/summary/
```

## Running the analysis

Run all scripts from the **repository root**.

```r
source("R/01_download_qc_integration.R")
source("R/02_add_sample_metadata.R")
source("R/03_cell_type_annotation.R")
source("R/04_emt_scoring.R")
source("R/05_emt_trajectory_monocle3.R")
source("R/06_nichenet_ligand_analysis.R")
source("R/07_emt_lncrna_analysis.R")
source("R/08_tcga_stad_survival.R")
source("R/09_generate_results_summary.R")
```

Because the complete pipeline downloads and processes large single-cell and TCGA datasets, substantial memory, storage, and runtime may be required.

## Reproducibility

The analysis scripts use fixed random seeds where applicable.

Each major stage also writes `sessionInfo()` output to `results/`.

A dependency lockfile (`renv.lock`) can be added after the pipeline is successfully tested in a clean R environment.

## Results

Selected lightweight figures and tables can be committed under:

```text
results/figures/
results/tables/
results/summary/
```

Large intermediate objects should remain local.

After a complete reproducibility run, the main computed overview will be available at:

```text
results/summary/analysis_summary.md
```

## Important interpretation notes

- EMT categories are computationally inferred states, not experimentally proven cell-state transitions.
- Monocle pseudotime represents an inferred transcriptional trajectory and should not be interpreted as direct temporal lineage tracing.
- NicheNet identifies ligand–target relationships supported by expression and prior knowledge; it does not establish causal signaling.
- TF–lncRNA correlations are associative.
- Candidate lncRNAs should be described as **EMT-associated** rather than as proven EMT-inducing regulators unless experimentally validated.
- The TCGA-STAD multigene model, if generated, is exploratory unless validated in an independent cohort.

## Data and version-control policy

See:

```text
data/README.md
.gitignore
```

for details.

The repository intentionally excludes:

```text
objects/
data/GSE183904_RAW/
data/TCGA_STAD/
data/nichenet/
data/gencode/*.gtf*
*.rds
```

## Citation

If you use this repository, please cite the original GSE183904 study and the relevant software/data resources.

A repository-level `CITATION.cff` will be included to provide a standardized citation for this analysis code.

## License

A repository license will be provided in `LICENSE`.

The license for this repository's analysis code does not override the licenses or usage terms of the original datasets, external databases, supplementary metadata, or third-party software.
