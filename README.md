# TCGA-BRCA Multi-Omics Factor Analysis

[![Live App](https://img.shields.io/badge/Live%20App-Shiny-1976D2?logo=rstudio)](https://arcturex.shinyapps.io/brca-navigator/)

**Discovering multi-omics latent factors linked to breast cancer survival** — integrating RNA-seq, DNA methylation, CNV, and RPPA proteomics from TCGA-BRCA (n=636) using MOFA2.

---

## Key Findings

- **15 latent factors** discovered; Factors 1 + 2 explain 47% of coordinated variation across omics
- **Factor 5** (anti-tumor immune response — PLA2G2D, TIGIT, ICOS, ZAP70): HR=0.63 [0.46–0.87], p=0.004 after age/stage adjustment
- **Factor 7**: +95.6 days restricted mean survival time (RMST) in high vs low group — a "delayed protection" signal invisible to Cox models
- **Factor 14**: strongest RSF variable importance (3× next best), top risk factor
- **Time-dependent AUC**: 0.705 at 3 years; cross-method consistency (Cox, RSF, RMST, KM, TimeROC)

[See full results →](https://arcturex.shinyapps.io/brca-navigator/)

---

## Pipeline

| Step | Script | Description |
|------|--------|-------------|
| **1. Download** | `step1_download_brca.R` | TCGA-BRCA data from UCSC Xena (RNA-seq, 450K methylation, Gistic2 CNV, RPPA, clinical) |
| **2. Preprocess** | `step2_build_mae.R` | KNN imputation, top-feature selection (RNA: 8k, Methyl: 8k, CNV: 3k), patient alignment, MultiAssayExperiment |
| **3. MOFA** | `step3_train_mofa.R` | MOFA2 with 15 factors, seed=42, maxiter=1000, fast convergence |
| **4. Exploration** | `step4_exploration.R` | Variance explained, top weights, PAM50 subtype distribution, GO/KEGG enrichment |
| **5. Survival** | `step5_survival.R` | Univariate Cox OS/RFS, Kaplan-Meier curves |
| **6. Advanced** | `step6.*.R` | Multivariable Cox (LASSO + age/stage), Random Survival Forest (1000 trees), TimeROC (1/3/5yr), RMST (5yr), GSEA, external GEO validation |
| **Summary** | `generate_all_outputs.R` + `summary_figures.R` | CSV results, analysis summary, publication-ready figures |

> **Note:** CNV was dropped during QC (0 features survived filtering). External validation (GSE20685, Affymetrix microarray) could not map probes to TCGA RNA-seq gene symbols.

---

## Interactive Dashboard

**[brca-navigator](https://arcturex.shinyapps.io/brca-navigator/)** — a Shiny app built with bs4Dash + plotly + visNetwork, deployed on shinyapps.io. Seven tabs:

| Tab | Purpose |
|-----|---------|
| The Story | Scrollytelling narrative of the analysis |
| Factor Explorer | Browse all 15 factors, variance, top weights, survival curves |
| Patient Map | UMAP embedding colored by factors, survival, PAM50 subtype |
| Gene Explorer | Search genes, view weights across factors, expression distributions, survival impact |
| Survival Simulator | Adjust factor scores interactively and see predicted survival curves |
| Analysis Flow | Visual pipeline diagram |
| About | Methods, reproducibility, references |

---

## Reproducibility

### Requirements

- **R 4.2+** with: MultiAssayExperiment, MOFA2, ggplot2, ComplexHeatmap, survival, survminer, glmnet, ranger, timeROC, survRM2, clusterProfiler, UCSCXenaTools, impute, data.table, patchwork, GEOquery
- **Python 3.12+** with: mofapy2 (via pip, used by MOFA2 via reticulate)
- **Bioconductor** packages installed via `BiocManager::install()`

### Run the Pipeline

```r
setwd("path/to/PO3")
# Set Python path for reticulate:
# Sys.setenv(RETICULATE_PYTHON = "path/to/python.exe")

source("scripts/step1_download_brca.R")   # ~10 min
source("scripts/step2_build_mae.R")        # ~5 min
source("scripts/step3_train_mofa.R")       # ~30–60 min
source("scripts/step4_exploration.R")
source("scripts/step5_survival.R")
source("scripts/step6.1_continuous_cox.R")
source("scripts/step6.2_rsf.R")
source("scripts/step6.3_timeroc.R")
source("scripts/step6.4_rmst.R")
source("scripts/step6.5_gsea.R")
source("scripts/step6.6_validation.R")     # may fail (see note)
source("scripts/generate_all_outputs.R")
source("scripts/summary_figures.R")
```

### Dashboard

```r
source("dashboard_app/data_prep.R")  # once
# Open dashboard_app/app.R in RStudio → Run App
# Or deploy:
# rsconnect::deployApp(appDir = "dashboard_app", appName = "brca-navigator")
```

---

## Directory Structure

```
.
├── Data/                  # Raw data (UCSC Xena download), MAE object
├── Data/GEO/              # GSE20685 for external validation
├── models/                # MOFA trained model (.hdf5 + .rds)
├── scripts/               # All pipeline R scripts
├── Results/               # CSVs: variance, survival, RSF, TimeROC, RMST, enrichment
├── logs/                  # Preprocessing summaries, Cox results
├── figures/               # All publication-ready figures
├── dashboard_app/         # Shiny dashboard (app.R, data_prep.R, www/custom.css)
│   └── data/              # Dashboard-ready datasets
└── docs/                  # Pipeline plan, analysis report, deployment notes
```

---

## Key Parameters

| Parameter | Value |
|-----------|-------|
| MOFA latent factors | 15 |
| Random seeds | 42 (MOFA, RSF, control selection) |
| RNA top features | 8,000 |
| Methylation top features | 8,000 |
| Methylation imputation | KNN k=10, drop >20% missing |
| RSF trees | 1,000 |
| LASSO CV folds | 10 |
| TimeROC timepoints | 365 / 1095 / 1825 days |
| RMST tau | 1,825 days (5 years) |
| GEO validation | GSE20685 (Affymetrix GPL570) |

---

## Results Summary

| Metric | Value |
|--------|-------|
| Patients | 636 (BRCA) |
| Omics integrated | 3 (RNA + Methylation + RPPA) |
| Factors discovered | 15 |
| Best protective factor (adjusted) | Factor 5, HR=0.63 [0.46–0.87], p=0.004 |
| Best risk factor (adjusted) | Factor 14, HR=1.44 [1.07–1.95], p=0.017 |
| RSF C-index | 0.594 |
| TimeROC AUC at 3 years | 0.705 |
| RMST best | Factor 7: +95.6 days (p=0.049) |

---

## License

This project is for educational and research purposes. TCGA data is available under [TCGA publication guidelines](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga/using-tcga/citing-tcga).
