## TCGA-BRCA Multi-Omics MOFA2 Project — Complete Plan

**Goal:** Integrate 4 omics layers (RNA-seq, DNA Methylation, CNV, RPPA Proteomics) from TCGA-BRCA using MOFA2 to discover multi-omics factors linked to survival.

**Why breast cancer instead of lung?** BRCA has ~1100 patients (vs 370 in LUSC) = more statistical power. Breast cancer has well-defined molecular subtypes (LumA, LumB, HER2-E, Basal) that are clinically actionable. BRCA also has proteomics (RPPA) data not available for LUSC.

**Working Directory:** `C:/Users/hp/Dev/playground/PO3` — all code assumes this is your working directory

---

## DAY 1 — Data Acquisition

---

## Step 0: Install Required R Packages

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "MultiAssayExperiment", "MOFA2", "ggplot2", "ComplexHeatmap",
    "dplyr", "tidyr", "tibble", "readr", "stringr",
    "clusterProfiler", "org.Hs.eg.db", "enrichplot",
    "survival", "survminer", "impute", "data.table",
    "corrplot", "RColorBrewer", "UCSCXenaTools", "patchwork",
    "glmnet", "ranger", "rms", "GEOquery"
))
```

**🧠 What you're doing:** Installing all the R packages needed for the entire analysis — from data download to survival analysis.

**🎯 Why:** These packages handle every step: `UCSCXenaTools` downloads data, `MOFA2` does the multi-omics integration, `survival/survminer` runs the survival models, and `glmnet/ranger/rms` are for advanced survival methods in Part 2.

---

## Step 1: Download TCGA-BRCA Data

```
Sys.setenv(RETICULATE_PYTHON = "C:/Users/hp/AppData/Local/Programs/Python/Python312/python.exe")
setwd("C:/Users/hp/Dev/playground/PO3")
library(data.table)
library(UCSCXenaTools)
library(stringr)
library(MultiAssayExperiment)

data_dir <- "Data/TCGA.BRCA.sampleMap"
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== Download each dataset ===\n")

# RNA-seq
downloadTCGA("BRCA", "Gene Expression RNASeq", "IlluminaHiSeq RNASeqV2 pancan normalized")

# DNA Methylation  
downloadTCGA("BRCA", "DNA Methylation", "Methylation450K")

# Copy Number
downloadTCGA("BRCA", "Gene Level Copy Number", "Gistic2")

# Clinical
downloadTCGA("BRCA", "Phenotype", "Clinical Information")

# RPPA — try different names since "Protein Expression" failed
cat("\nTrying RPPA download with alternative names...\n")
tryCatch({
    downloadTCGA("BRCA", "Protein Expression", "RPPA")
}, error = function(e) cat("RPPA not found via downloadTCGA:", e$message, "\n"))

# === Copy files from temp to project directory ===
cat("\n=== Finding downloaded files ===\n")
tmp_dirs <- list.dirs(tempdir(), recursive = TRUE)
samplemap_dirs <- tmp_dirs[grepl("TCGA.BRCA.sampleMap$", tmp_dirs)]
if (length(samplemap_dirs) > 0) {
    src <- samplemap_dirs[1]
    files <- list.files(src, full.names = TRUE)
    file.copy(files, data_dir, overwrite = TRUE)
    cat("Copied", length(files), "files from temp:\n", paste(basename(files), collapse = "\n "), "\n")
} else {
    cat("No temp directory found — files may already be in Data/TCGA.BRCA.sampleMap/\n")
}

# === List what we have ===
cat("\n=== Files in project directory ===\n")
local_files <- list.files(data_dir, full.names = TRUE)
if (length(local_files) == 0) {
    stop("No files found. Try downloading manually or check internet connection.")
}
print(basename(local_files))

# === Load by pattern ===
cat("\n=== Loading data ===\n")
find_file <- function(pattern) {
    f <- list.files(data_dir, pattern = pattern, full.names = TRUE)
    if (length(f) == 0) stop("No file matching: ", pattern)
    cat("  Found:", basename(f[1]), "\n")
    f[1]
}

rna_mat <- NULL; meth_mat <- NULL; cnv_mat <- NULL; rppa_mat <- NULL

if (length(list.files(data_dir, pattern = "HiSeq")) > 0) {
    d <- fread(find_file("HiSeq"))
    rna_mat <- as.matrix(d[, -1]); rownames(rna_mat) <- d$sample
}
if (length(list.files(data_dir, pattern = "Methyl")) > 0) {
    d <- fread(find_file("Methyl"))
    meth_mat <- as.matrix(d[, -1]); rownames(meth_mat) <- d$sample
}
if (length(list.files(data_dir, pattern = "Gistic2")) > 0) {
    d <- fread(find_file("Gistic2"))
    cnv_mat <- as.matrix(d[, -1]); rownames(cnv_mat) <- d$sample
}
rppa_files <- list.files(data_dir, pattern = "RPPA|protein", ignore.case = TRUE)
if (length(rppa_files) > 0) {
    d <- fread(file.path(data_dir, rppa_files[1]))
    rppa_mat <- as.matrix(d[, -1]); rownames(rppa_mat) <- d$sample
    cat("  RPPA loaded:", dim(rppa_mat), "\n")
} else {
    cat("  RPPA: not available — continuing with 3 omics\n")
}

clinical <- fread(find_file("clinicalMatrix"))

cat("\nLoaded:\n")
if (!is.null(rna_mat)) cat("  RNA-seq:", nrow(rna_mat), "x", ncol(rna_mat), "\n")
if (!is.null(meth_mat)) cat("  Methyl:", nrow(meth_mat), "x", ncol(meth_mat), "\n")
if (!is.null(cnv_mat)) cat("  CNV:", nrow(cnv_mat), "x", ncol(cnv_mat), "\n")
if (!is.null(rppa_mat)) cat("  RPPA:", nrow(rppa_mat), "x", ncol(rppa_mat), "\n")
cat("  Clinical:", nrow(clinical), "x", ncol(clinical), "\n")

cat("\n=== Finding overlapping patients ===\n")
extract_patient_id <- function(barcodes) str_sub(barcodes, 1, 12)

available <- list()
if (!is.null(rna_mat)) available$RNA <- rna_mat
if (!is.null(meth_mat)) available$Methyl <- meth_mat
if (!is.null(cnv_mat)) available$CNV <- cnv_mat
if (!is.null(rppa_mat)) available$RPPA <- rppa_mat

common <- Reduce(intersect, lapply(available, function(m) unique(extract_patient_id(colnames(m)))))
cat("Patients with all available omics:", length(common), "\n")

saveRDS(list(rna = rna_mat, meth = meth_mat, cnv = cnv_mat, rppa = rppa_mat,
             clinical = clinical, common = common), "Data/brca_raw.rds")
cat("Saved raw data to Data/brca_raw.rds\n")
```

**🧠 What you're doing:** Downloading 5 types of data for ~1100 breast cancer patients from UCSC Xena — gene expression (RNA), DNA methylation, copy number alterations (CNV), protein expression (RPPA), and clinical information (survival, stage, etc.).

**🎯 Why:** MOFA2 needs multiple data layers from the SAME patients to find coordinated patterns. Each omics layer captures a different biological angle: RNA = what genes are active, Methylation = which genes are silenced, CNV = which genes are gained/lost, RPPA = which proteins are actually present.

---

## DAY 2 — Data Preprocessing

---

## Step 2.1: Handle Duplicate Samples

```
extract_patient_id <- function(barcodes) str_sub(barcodes, 1, 12)
deduplicate <- function(mat) {
    pids <- extract_patient_id(colnames(mat))
    dup <- duplicated(pids)
    if (sum(dup) > 0) cat("Removing", sum(dup), "duplicate samples\n")
    mat[, !dup, drop = FALSE]
}

top_n <- function(mat, n) {
    vars <- apply(mat, 1, var, na.rm = TRUE)
    n <- min(n, length(vars))
    mat[head(order(vars, decreasing = TRUE), n), , drop = FALSE]
}
rna_top <- top_n(rna_mat, 8000)
meth_top <- top_n(meth_imputed, 8000)
cnv_top <- top_n(cnv_mat, 3000)

# Deduplicate + subset
rna_f <- deduplicate(rna_top)
meth_f <- deduplicate(meth_top)
cnv_f <- deduplicate(cnv_top)
rppa_f <- deduplicate(rppa_mat)
```

**🧠 What you're doing:** Some patients have multiple tumor samples (e.g., two biopsies from different sites). We keep one sample per patient.

**🎯 Why:** Duplicate samples would bias the analysis — MOFA2 assumes one column = one patient. Multiple samples from the same patient would look like independent patients and create false structure.

---

## Step 2.2: Impute Missing Methylation Values and Select Top Variable Features

```
missing_pct <- rowSums(is.na(meth_mat)) / ncol(meth_mat)
meth_imputed <- impute.knn(as.matrix(meth_mat[missing_pct < 0.2, ]), k=10)$data

# Top variable features (robust — uses integer indexing, avoids name mismatch)
top_n <- function(mat, n) {
    vars <- apply(mat, 1, var, na.rm = TRUE)
    n <- min(n, length(vars))
    mat[head(order(vars, decreasing = TRUE), n), , drop = FALSE]
}
rna_top <- top_n(rna_mat, 8000)
meth_top <- top_n(meth_imputed, 8000)
cnv_top <- top_n(cnv_mat, 3000)
```

**🧠 What you're doing:** Picking the most informative features from each data layer. Genes/probes that barely change across patients carry little information. RPPA has only ~200 proteins, so we keep them all. RNA and methylation are reduced to 8000, CNV to 3000. Methylation is also imputed to fill in missing values using k-NN.

**🎯 Why:** MOFA2 training time scales with the number of features. 8000 genes + 8000 probes + 3000 CNV + 200 RPPA = ~19,200 features, which is manageable. Keeping everything (e.g., 20,000 genes) would be slow and noisy. Methylation imputation is needed because MOFA2 cannot handle NAs well across the whole matrix.

---

## Step 2.3: Align All Patients and Build MAE

```
subset_and_align <- function(mat, ref) {
    pids <- extract_patient_id(colnames(mat))
    mat <- mat[, pids %in% ref, drop=FALSE]
    pids <- extract_patient_id(colnames(mat))
    mat[, order(match(pids, ref)), drop=FALSE]
}
common <- sort(common)
rna_f <- subset_and_align(rna_f, common)
meth_f <- subset_and_align(meth_f, common)
cnv_f <- subset_and_align(cnv_f, common)
rppa_f <- subset_and_align(rppa_f, common)

# Remove duplicated feature names (drop=FALSE prevents vector collapse)
rna_f <- rna_f[!duplicated(rownames(rna_f)), , drop = FALSE]
meth_f <- meth_f[!duplicated(rownames(meth_f)), , drop = FALSE]
cnv_f <- cnv_f[!duplicated(rownames(cnv_f)), , drop = FALSE]
rppa_f <- rppa_f[!duplicated(rownames(rppa_f)), , drop = FALSE]

colnames(rna_f) <- colnames(meth_f) <- colnames(cnv_f) <- colnames(rppa_f) <- common

stopifnot(all(colnames(rna_f) == common))
stopifnot(all(colnames(meth_f) == common))
stopifnot(all(colnames(cnv_f) == common))
stopifnot(all(colnames(rppa_f) == common))

# Build MAE
clinical_mat <- as.data.frame(clinical)
cid <- extract_patient_id(clinical_mat[[1]])
dup <- duplicated(cid)
clinical_mat <- clinical_mat[!dup, , drop=FALSE]
rownames(clinical_mat) <- cid[!dup]
clinical_mat <- clinical_mat[common, , drop=FALSE]

mae <- MultiAssayExperiment(
    experiments = ExperimentList(list(
        RNA = as.matrix(rna_f), Methyl = as.matrix(meth_f),
        CNV = as.matrix(cnv_f), RPPA = as.matrix(rppa_f))),
    colData = DataFrame(clinical_mat))

saveRDS(mae, "Data/brca_mae.rds")
summary(mae)
cat("Saved MAE with", ncol(mae), "patients\n")
```

**🧠 What you're doing:** Building a single organized R object (MultiAssayExperiment) that contains all 4 omics layers + clinical data, with every column matching the same patient. It's like a binder with 4 tabs, each tab sorted the same way.

**🎯 Why:** The MAE is the standard format MOFA2 expects. By building it once and saving to `Data/brca_mae.rds`, you can reload it in any later session without repeating the alignment.

---

## Step 2.4: Preprocessing Summary

```
cat("=== MAE Summary ===\n")
print(mae)
cat("\n=== Patients per Omics ===\n")
cat("RNA + Methyl + CNV + RPPA:", ncol(mae), "\n")
cat("Total features:", sum(sapply(experiments(mae), nrow)), "\n")
cat("Clinical columns:", ncol(colData(mae)), "\n")
```

**🧠 What you're doing:** Checking the final preprocessed data — how many patients, how many features per omics, what clinical data is available.

**🎯 Why:** Before training MOFA, you need to confirm the numbers look right. Expect ~400-600 patients with all 4 omics. If < 300, consider dropping RPPA.

---

## DAY 3 — MOFA Training

---

## Step 3.1: Install MOFA2 Python Backend (Run Once)

```
reticulate::py_install("mofapy2", pip = TRUE)
```

**🧠 What you're doing:** Installing the Python package `mofapy2` that MOFA2's R package uses as its computational engine.

**🎯 Why:** The R package MOFA2 is a wrapper that calls a Python backend (`mofapy2`) for the actual model fitting. This step only needs to run once.

---

## Step 3.2: Train MOFA Model

```
Sys.setenv(RETICULATE_PYTHON = "C:/Users/hp/AppData/Local/Programs/Python/Python312/python.exe")
setwd("C:/Users/hp/Dev/playground/PO3")
library(MOFA2)
library(MultiAssayExperiment)

mae <- readRDS("Data/brca_mae.rds")

# Build MOFA data list
mofa_data <- list()
for (nm in names(experiments(mae))) {
    mat <- as.matrix(assay(experiments(mae)[[nm]]))
    dup <- duplicated(rownames(mat))
    if (any(dup)) mat <- mat[!dup, ]
    mofa_data[[nm]] <- mat
}
cat("Views:", names(mofa_data), "\nFeatures:", sapply(mofa_data, nrow), "\nSamples:", unique(sapply(mofa_data, ncol)), "\n")

# Drop empty views (e.g. CNV with 0 features)
mofa_data <- mofa_data[sapply(mofa_data, nrow) > 0]
cat("Keeping views:", names(mofa_data), "\n")

mofa_obj <- create_mofa(mofa_data)

model_opts <- get_default_model_options(mofa_obj)
model_opts$num_factors <- 15
model_opts$spikeslab_factors <- FALSE
model_opts$ard_factors <- TRUE

train_opts <- get_default_training_options(mofa_obj)
train_opts$convergence_mode <- "fast"
train_opts$maxiter <- 1000
train_opts$seed <- 42

mofa_obj <- prepare_mofa(mofa_obj,
    data_options = get_default_data_options(mofa_obj),
    model_options = model_opts,
    training_options = train_opts)

dir.create("models", showWarnings=FALSE)
mofa_trained <- run_mofa(mofa_obj, outfile="models/brca_mofa2.hdf5")
saveRDS(mofa_trained, "models/brca_mofa_trained.rds")
cat("MOFA training complete.\n")
```

**🧠 What you're doing:** Training the MOFA model — this is the core computation. MOFA learns 15 latent factors that explain variation across ALL 4 omics simultaneously. Each factor captures a pattern of coordinated molecular changes (e.g., "Factor 1 might capture Luminal vs Basal differences").

**🎯 Why:** This is the central analysis. The factors become your new features — each patient gets 15 scores (one per factor) that summarize their multi-omics profile. These scores are used in all downstream analyses.

**⏱ Timing:** ~30-60 minutes depending on CPU. BRCA has more patients and more omics than LUSC.

---

## Step 3.3: Quick Quality Check

```
mofa_object <- readRDS("models/brca_mofa_trained.rds")

# Variance explained per factor
var_plot <- plot_variance_explained(mofa_object)
print(var_plot)

# Total variance explained per view
tot_variance <- plot_variance_explained(mofa_object, plot_total = TRUE)
print(tot_variance)

ve <- get_variance_explained(mofa_object)
cat("Factors with >2% variance (total across views):\n",
    sum(rowSums(ve$r2_per_factor$group1) > 0.02), "/ 15\n")
```

**🧠 What you're doing:** Checking if the factors actually capture meaningful variation. Each factor should explain at least a few percent of variance. Factors near 0% are noise and can be dropped.

**🎯 Why:** Not all 15 factors will be useful. Typically 5-10 factors capture most of the signal. We use this to decide which factors to test in survival analysis.

---

## DAY 4 — Survival Analysis

---

## Step 4.1: Prepare Clinical Data with PAM50 Subtypes (Exploration)

```
mofa_object <- readRDS("models/brca_mofa_trained.rds")
clinical <- readRDS("Data/brca_mae.rds")

# Extract factor scores
factor_scores <- get_factors(mofa_object, as.data.frame = TRUE)
head(factor_scores)

# Prepare clinical data
clinical_data <- colData(clinical) %>% as.data.frame()
colnames(clinical_data)

# Check for PAM50 subtypes
grep("PAM50|Subtype|BRCA_Subtype", colnames(clinical_data), value = TRUE, ignore.case = TRUE)

# Check for survival columns
grep("days_to_death|days_to_last_followup|vital_status|OS|recurrence|progression",
    colnames(clinical_data), value = TRUE, ignore.case = TRUE)
```

**🧠 What you're doing:** Peeking at the clinical data to find PAM50 subtype and survival columns needed for downstream analysis.

**🎯 Why:** BRCA has well-defined subtypes (LumA, LumB, HER2-E, Basal, Normal-like) that correlate strongly with survival. These are important covariates.

---

## Step 4.2: PAM50 Subtype Distribution (Exploration)

```
# Identify PAM50 column (adjust name if needed)
pam50_col <- grep("PAM50|Subtype|BRCA_Subtype", colnames(clinical_data), value = TRUE, ignore.case = TRUE)[1]
cat("PAM50 column:", pam50_col, "\n")

pam50 <- clinical_data[[pam50_col]]
print(table(pam50, useNA = "ifany"))

# Barplot
barplot(table(pam50), main = "PAM50 Subtype Distribution", las = 2, col = "steelblue")
```

**🧠 What you're doing:** Visualizing the distribution of breast cancer subtypes. LumA is typically most common (~40-50%), followed by LumB and Basal.

**🎯 Why:** Different subtypes have very different survival outcomes. Subtype is a critical covariate in Cox models to isolate factor effects from known clinical effects.

---

## Step 5.0: Load MOFA and Create Survival Dataset

```
Sys.setenv(RETICULATE_PYTHON = "C:/Users/hp/AppData/Local/Programs/Python/Python312/python.exe")
setwd("C:/Users/hp/Dev/playground/PO3")
library(stringr)
library(MultiAssayExperiment)
library(MOFA2)
library(survival)
library(survminer)

mae <- readRDS("Data/brca_mae.rds")
mofa_trained <- readRDS("models/brca_mofa_trained.rds")
clinical_mat <- as.data.frame(colData(mae))

factor_scores <- get_factors(mofa_trained)[[1]]
colnames(factor_scores) <- paste0("Factor", 1:ncol(factor_scores))

# === OS time/event (used by Steps 5.1-5.3) ===
os_time <- as.numeric(ifelse(clinical_mat$Vital_Status_nature2012=="DECEASED", clinical_mat$Days_to_date_of_Death_nature2012, clinical_mat$Days_to_Date_of_Last_Contact_nature2012))
os_event <- ifelse(clinical_mat$Vital_Status_nature2012=="DECEASED", 1, 0)
cat("OS:", sum(os_event), "deaths,", sum(!os_event), "censored\n")
```

**🧠 What you're doing:** Loading the trained MOFA model and creating a clean dataset with factor scores + clinical variables ready for survival analysis.

**🎯 Why:** The factor scores are the key MOFA output — each patient has 15 numbers representing their multi-omics profile. These become predictors in survival models. Run once before Steps 5.1-5.3.

---

## Step 5.1: Overall Survival — Univariate Cox Regression

```
# Run after Step 5.0
cox_os <- data.frame()
for (f in colnames(factor_scores)) {
    cox <- coxph(as.formula(paste0("Surv(os_time, os_event) ~ `", f, "`")),
                 data=data.frame(os_time, os_event, factor_scores))
    s <- summary(cox)
    cox_os <- rbind(cox_os, data.frame(Factor=f, Endpoint="OS",
        HR=exp(coef(cox)), p_value=coef(s)[,"Pr(>|z|)"],
        CI_lower=exp(confint(cox)[1]), CI_upper=exp(confint(cox)[2])))
}
cox_os$p_adj <- p.adjust(cox_os$p_value, method="BH")
write.csv(cox_os, "logs/cox_os.csv", row.names=FALSE)
print(cox_os)
```

**🧠 What you're doing:** Testing each MOFA factor one-at-a-time for association with overall survival. HR > 1 = higher factor score = worse survival. HR < 1 = higher score = better survival.

**🎯 Why:** This tells us which multi-omics patterns are clinically relevant. BH correction adjusts for testing 15 factors.

---

## Step 5.2: Overall Survival — Significant Factor KM Curves

```
# Run after Step 5.0 + Step 5.1
best <- cox_os$Factor[which.min(cox_os$p_value)]
grp <- ifelse(factor_scores[,best] > median(factor_scores[,best], na.rm=TRUE), "High", "Low")
fit <- survfit(Surv(os_time, os_event) ~ grp)
p <- ggsurvplot(fit, pval=TRUE, risk.table=TRUE, title=paste(best, "(OS)"))
dir.create("figures", showWarnings=FALSE)
png(paste0("figures/km_", best, "_OS.png"), width=8, height=6, units="in", res=300)
print(p); dev.off()
```

**🧠 What you're doing:** Creating Kaplan-Meier survival curves for each significant factor, splitting patients into High vs Low groups by the median factor score.

**🎯 Why:** KM curves visualize the survival difference. The p-value from the log-rank test confirms the Cox result.

---

## Step 5.3: Recurrence-Free Survival (RFS)

```
# Run after Step 5.0
rfs_time <- os_time; rfs_event <- os_event
rfs_cols <- grep("recurrence|progression|recur|progres", names(clinical_mat), value=TRUE, ignore.case=TRUE)
for (col in rfs_cols) {
    v <- as.numeric(clinical_mat[[col]])
    if (sum(!is.na(v)) > 100 && median(v, na.rm=TRUE) < 5000) { rfs_time <- v; break }
}
for (col in rfs_cols) {
    v <- suppressWarnings(as.numeric(clinical_mat[[col]]))
    if (sum(!is.na(v)) > 50 && all(v %in% c(0,1,NA))) { rfs_event <- v; break }
    if (is.character(clinical_mat[[col]])) {
        v <- as.numeric(grepl("recur|progress|dead|deceased", clinical_mat[[col]], ignore.case=TRUE))
        rfs_event <- v; break
    }
}
na_rfs <- which(is.na(rfs_time)); rfs_time[na_rfs] <- os_time[na_rfs]; rfs_event[na_rfs] <- os_event[na_rfs]

cox_rfs <- data.frame()
for (f in colnames(factor_scores)) {
    cox <- coxph(as.formula(paste0("Surv(rfs_time, rfs_event) ~ `", f, "`")),
                 data=data.frame(rfs_time, rfs_event, factor_scores))
    s <- summary(cox)
    cox_rfs <- rbind(cox_rfs, data.frame(Factor=f, Endpoint="RFS",
        HR=exp(coef(cox)), p_value=coef(s)[,"Pr(>|z|)"],
        CI_lower=exp(confint(cox)[1]), CI_upper=exp(confint(cox)[2])))
}
cox_rfs$p_adj <- p.adjust(cox_rfs$p_value, method="BH")
write.csv(cox_rfs, "logs/cox_rfs.csv", row.names=FALSE)
print(cox_rfs)
```

**🧠 What you're doing:** Same analysis but for recurrence-free survival — testing which factors predict tumor recurrence after treatment.

**🎯 Why:** RFS is clinically important — it tells us which molecular patterns drive recurrence, not just overall death. Less events = lower power, so interpret with caution.

---

## DAY 6 — Multivariable Cox and Advanced Methods

---

## Step 6.1: Multivariable Cox with PAM50 + LASSO Selection

```
Sys.setenv(RETICULATE_PYTHON = "C:/Users/hp/AppData/Local/Programs/Python/Python312/python.exe")
setwd("C:/Users/hp/Dev/playground/PO3")
library(stringr)
library(MultiAssayExperiment)
library(MOFA2)
library(survival)
library(glmnet)
library(ggplot2)
library(survminer)

mae <- readRDS("Data/brca_mae.rds")
mofa_trained <- readRDS("models/brca_mofa_trained.rds")
clinical_mat <- as.data.frame(colData(mae))

factor_scores <- get_factors(mofa_trained)[[1]]
colnames(factor_scores) <- paste0("Factor", 1:ncol(factor_scores))

os_time <- as.numeric(ifelse(clinical_mat$Vital_Status_nature2012=="DECEASED", clinical_mat$Days_to_date_of_Death_nature2012, clinical_mat$Days_to_Date_of_Last_Contact_nature2012))
os_event <- ifelse(clinical_mat$Vital_Status_nature2012=="DECEASED", 1, 0)
age <- as.numeric(clinical_mat$Age_at_Initial_Pathologic_Diagnosis_nature2012)
stage <- ifelse(grepl("III|IV", clinical_mat$AJCC_Stage_nature2012), "Late", "Early")
cat("OS:", sum(os_event), "deaths\n")

# A: LASSO
set.seed(42)
cv <- cv.glmnet(scale(factor_scores), Surv(os_time, os_event), family="cox", alpha=1, nfolds=10)
sel <- names(which(as.matrix(coef(cv, s="lambda.1se")) != 0))
if (length(sel)==0) sel <- names(which(as.matrix(coef(cv, s="lambda.min")) != 0))
cat("LASSO:", if(length(sel)>0) paste(sel, collapse=", ") else "none\n")

# B: Continuous univariate + C: Adjusted
results <- data.frame()
results_adj <- data.frame()
clin <- data.frame(os_time, os_event, age=scale(age), stage=factor(stage), scale(factor_scores))
clin <- clin[complete.cases(clin), ]

for (f in paste0("Factor", 1:15)) {
    cox <- coxph(as.formula(paste0("Surv(os_time, os_event) ~ `", f, "`")),
                 data=data.frame(os_time, os_event, scale(factor_scores)))
    s <- summary(cox)
    results <- rbind(results, data.frame(Factor=f, HR=exp(coef(cox)), p_value=coef(s)[,"Pr(>|z|)"],
        CI_lower=exp(confint(cox)[1]), CI_upper=exp(confint(cox)[2])))
    cox2 <- tryCatch(coxph(as.formula(paste0("Surv(os_time, os_event) ~ `", f, "` + age + stage")), data=clin), error=function(e) NULL)
    if (!is.null(cox2)) {
        s2 <- summary(cox2)
        results_adj <- rbind(results_adj, data.frame(Factor=f, HR=exp(coef(cox2))[1], p_value=coef(s2)[1,"Pr(>|z|)"],
            CI_lower=exp(confint(cox2))[f,1], CI_upper=exp(confint(cox2))[f,2]))
    }
}
results$p_adj <- p.adjust(results$p_value, method="BH")
results_adj$p_adj <- p.adjust(results_adj$p_value, method="BH")

cat("\nUnadjusted:\n"); for (i in 1:nrow(results)) cat(sprintf("%-10s HR=%6.3f [%5.3f-%5.3f] p=%6.4f adj=%6.4f\n", results$Factor[i], results$HR[i], results$CI_lower[i], results$CI_upper[i], results$p_value[i], results$p_adj[i]))
cat("\nAdjusted for age+stage:\n"); for (i in 1:nrow(results_adj)) cat(sprintf("%-10s HR=%6.3f [%5.3f-%5.3f] p=%6.4f adj=%6.4f\n", results_adj$Factor[i], results_adj$HR[i], results_adj$CI_lower[i], results_adj$CI_upper[i], results_adj$p_value[i], results_adj$p_adj[i]))

write.csv(results, "logs/cox_continuous_os.csv", row.names=FALSE)
write.csv(results_adj, "logs/cox_adjusted_clinical_os.csv", row.names=FALSE)

# D: KM curves (median + tertile for top 3)
top3 <- results$Factor[order(results$p_value)][1:3]
for (f in top3) {
    med <- median(factor_scores[,f], na.rm=TRUE)
    grp <- ifelse(factor_scores[,f] > med, "High", "Low")
    fit <- survfit(Surv(os_time, os_event) ~ grp)
    p <- ggsurvplot(fit, pval=TRUE, risk.table=TRUE, title=paste(f, "- median (p=", round(results$p_value[results$Factor==f],4), ")"))
    png(paste0("figures/km_median_", f, ".png"), width=8, height=6, units="in", res=300); print(p); dev.off()

    q <- quantile(factor_scores[,f], c(1/3, 2/3), na.rm=TRUE)
    grp <- factor(ifelse(factor_scores[,f] > q[2], "High", ifelse(factor_scores[,f] <= q[1], "Low", "Mid")), levels=c("Low","Mid","High"))
    idx <- grp != "Mid"
    fit <- survfit(Surv(os_time[idx], os_event[idx]) ~ droplevels(grp[idx]))
    p <- ggsurvplot(fit, pval=TRUE, risk.table=TRUE, title=paste(f, "- extreme tertiles"))
    png(paste0("figures/km_tertile_", f, ".png"), width=8, height=6, units="in", res=300); print(p); dev.off()
}
cat("Saved figures.\n")
```

**🧠 What you're doing:** Testing MOFA factors while controlling for age and stage. LASSO selects which factors are most predictive.

**🎯 Why:** Univariate Cox can be misleading if a factor is just correlated with age/stage. This tells us if factors add independent prognostic value beyond standard clinical variables. Full standalone script at `scripts/step6.1_continuous_cox.R`.

---

## Step 6.2: Random Survival Forests

```
Sys.setenv(RETICULATE_PYTHON = "C:/Users/hp/AppData/Local/Programs/Python/Python312/python.exe")
setwd("C:/Users/hp/Dev/playground/PO3")
library(MultiAssayExperiment)
library(MOFA2)
library(survival)
library(ranger)

mae <- readRDS("Data/brca_mae.rds")
mofa_trained <- readRDS("models/brca_mofa_trained.rds")

factor_scores <- get_factors(mofa_trained)[[1]]
colnames(factor_scores) <- paste0("Factor", 1:ncol(factor_scores))

clinical_mat <- as.data.frame(colData(mae))
os_time <- as.numeric(ifelse(clinical_mat$Vital_Status_nature2012=="DECEASED",
    clinical_mat$Days_to_date_of_Death_nature2012,
    clinical_mat$Days_to_Date_of_Last_Contact_nature2012))
os_event <- ifelse(clinical_mat$Vital_Status_nature2012=="DECEASED", 1, 0)

surv_df <- data.frame(time=os_time, status=os_event, factor_scores)
surv_df <- surv_df[complete.cases(surv_df), ]

rsf_model <- ranger(
    Surv(time, status) ~ .,
    data = surv_df,
    importance = "permutation",
    num.trees = 1000,
    seed = 42
)

vimp <- data.frame(
    Factor = names(rsf_model$variable.importance),
    Importance = rsf_model$variable.importance
)
vimp <- vimp[order(vimp$Importance, decreasing = TRUE), ]
print(vimp)

cat("RSF C-index:", round(rsf_model$prediction.error, 3), "\n")

dir.create("Results", showWarnings=FALSE)
write.csv(vimp, "Results/brca_rsf_importance.csv", row.names = FALSE)

dir.create("figures", showWarnings=FALSE)
png("figures/brca_rsf_importance.png", width = 8, height = 5, units = "in", res = 150)
par(mar = c(4, 8, 3, 2))
barplot(vimp$Importance, names.arg = vimp$Factor, horiz = TRUE, las = 1,
    col = "steelblue", main = "Random Survival Forest — Variable Importance")
dev.off()
cat("Done.\n")
```

**🧠 What you're doing:** Running a Random Survival Forest — a machine learning method that doesn't assume linear relationships between factors and survival.

**🎯 Why:** Cox regression assumes the log-hazard changes linearly with each factor. RSF can capture non-linear effects and interactions, giving a more complete picture of which factors matter. Full standalone script at `scripts/step6.2_rsf.R`.

---

## Step 6.3: Time-Dependent ROC Curves

```
Sys.setenv(RETICULATE_PYTHON = "C:/Users/hp/AppData/Local/Programs/Python/Python312/python.exe")
setwd("C:/Users/hp/Dev/playground/PO3")
library(MultiAssayExperiment)
library(MOFA2)
library(survival)
library(timeROC)

mae <- readRDS("Data/brca_mae.rds")
mofa_trained <- readRDS("models/brca_mofa_trained.rds")

factor_scores <- get_factors(mofa_trained)[[1]]
colnames(factor_scores) <- paste0("Factor", 1:ncol(factor_scores))

clinical_mat <- as.data.frame(colData(mae))
os_time <- as.numeric(ifelse(clinical_mat$Vital_Status_nature2012=="DECEASED",
    clinical_mat$Days_to_date_of_Death_nature2012,
    clinical_mat$Days_to_Date_of_Last_Contact_nature2012))
os_event <- ifelse(clinical_mat$Vital_Status_nature2012=="DECEASED", 1, 0)

surv_df <- data.frame(time=os_time, status=os_event, factor_scores)
surv_df <- surv_df[complete.cases(surv_df), ]

times <- c(365, 1095, 1825)

roc_result <- timeROC(
    T = surv_df$time,
    delta = surv_df$status,
    marker = as.matrix(surv_df[, paste0("Factor", 1:15)]),
    cause = 1,
    times = times,
    iid = TRUE
)

cat("Time-dependent AUC:\n")
for (i in seq_along(times)) {
    cat("Year", c(1, 3, 5)[i], ":", round(roc_result$AUC[i], 3), "\n")
}

dir.create("figures", showWarnings=FALSE)
png("figures/brca_timeroc_auc.png", width = 8, height = 6, units = "in", res = 150)
plot(roc_result, time = times[1], col = "red", title = "")
plot(roc_result, time = times[2], col = "blue", add = TRUE)
plot(roc_result, time = times[3], col = "green", add = TRUE)
legend("bottomright", legend = c("1 Year", "3 Years", "5 Years"),
    col = c("red", "blue", "green"), lty = 1, lwd = 2)
title(main = "Time-Dependent ROC — MOFA Factors")
dev.off()
cat("Done.\n")
```

**🧠 What you're doing:** Calculating how well MOFA factors predict survival at specific time points (1, 3, 5 years).

**🎯 Why:** AUC tells us the discriminative power at clinically relevant timepoints. 5-year AUC > 0.7 would be clinically useful. Full standalone script at `scripts/step6.3_timeroc.R`.

---

## Step 6.4: Restricted Mean Survival Time (RMST)

```
Sys.setenv(RETICULATE_PYTHON = "C:/Users/hp/AppData/Local/Programs/Python/Python312/python.exe")
setwd("C:/Users/hp/Dev/playground/PO3")
library(MultiAssayExperiment)
library(MOFA2)
library(survival)
library(survRM2)

mae <- readRDS("Data/brca_mae.rds")
mofa_trained <- readRDS("models/brca_mofa_trained.rds")

factor_scores <- get_factors(mofa_trained)[[1]]
colnames(factor_scores) <- paste0("Factor", 1:ncol(factor_scores))

clinical_mat <- as.data.frame(colData(mae))
os_time <- as.numeric(ifelse(clinical_mat$Vital_Status_nature2012=="DECEASED",
    clinical_mat$Days_to_date_of_Death_nature2012,
    clinical_mat$Days_to_Date_of_Last_Contact_nature2012))
os_event <- ifelse(clinical_mat$Vital_Status_nature2012=="DECEASED", 1, 0)

surv_df <- data.frame(time=os_time, status=os_event, factor_scores)
surv_df <- surv_df[complete.cases(surv_df), ]

rmst_result <- rmst2(
    time = surv_df$time,
    status = surv_df$status,
    arm = rowMeans(surv_df[, paste0("Factor", 1:15)]),
    tau = 1825
)
print(rmst_result)

cat("\nPer-factor RMST:\n")
for (f in paste0("Factor", 1:15)) {
    arm_high <- ifelse(surv_df[[f]] > median(surv_df[[f]]), 1, 0)
    rmst_f <- rmst2(surv_df$time, surv_df$status, arm = arm_high, tau = 1825)
    cat(f, "- RMST difference:", round(rmst_f$unadjusted.result[1, 1], 1),
        "days, p:", format(rmst_f$unadjusted.result[1, 2], digits = 3), "\n")
}
cat("Done.\n")
```

**🧠 What you're doing:** RMST estimates the mean survival time within a restricted window (5 years), comparing high vs low factor groups.

**🎯 Why:** RMST is more interpretable than HR — it tells "how many more days of life" the low group has. Much easier to explain to clinicians. Full standalone script at `scripts/step6.4_rmst.R`.

---

## Step 6.5: Gene Set Enrichment of Top Factor Weights

```
Sys.setenv(RETICULATE_PYTHON = "C:/Users/hp/AppData/Local/Programs/Python/Python312/python.exe")
setwd("C:/Users/hp/Dev/playground/PO3")
library(MultiAssayExperiment)
library(MOFA2)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)

mofa_trained <- readRDS("models/brca_mofa_trained.rds")
mae_obj <- readRDS("Data/brca_mae.rds")

fs <- get_factors(mofa_trained, as.data.frame = TRUE)
fw <- pivot_wider(fs, names_from = "factor", values_from = "value")
fm <- as.data.frame(fw)
rownames(fm) <- fm$sample
fm$sample <- NULL

rna_data <- as.matrix(assay(mae_obj, "RNA"))

cat("Correlating Factor 1 with all genes...\n")
correlations <- apply(rna_data, 1, function(g) cor(g, fm$Factor1, use = "complete.obs"))
top_genes <- names(sort(abs(correlations), decreasing = TRUE)[1:500])

gene_df <- bitr(top_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

go_result <- enrichGO(
    gene = gene_df$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
)

dir.create("Results", showWarnings=FALSE)
dir.create("figures", showWarnings=FALSE)

if (!is.null(go_result) && nrow(go_result) > 0) {
    png("figures/brca_factor1_go_enrichment.png", width = 10, height = 8, units = "in", res = 150)
    print(dotplot(go_result, showCategory = 15))
    dev.off()
    write.csv(as.data.frame(go_result), "Results/brca_factor1_go_enrichment.csv", row.names = FALSE)
    cat("GSEA complete —", nrow(go_result), "enriched terms found.\n")
} else {
    cat("No enriched terms found at p<0.05.\n")
}
```

**🧠 What you're doing:** Finding which biological pathways (GO terms) are associated with Factor 1 by correlating factor scores with gene expression and testing enriched pathways.

**🎯 Why:** This gives biological meaning to statistical factors. Factor 1 might represent "immune response" or "cell cycle" — understanding this helps interpret what MOFA discovered. Full standalone script at `scripts/step6.5_gsea.R`.

---

## Step 6.6: External Validation in Independent GEO Dataset

```
Sys.setenv(RETICULATE_PYTHON = "C:/Users/hp/AppData/Local/Programs/Python/Python312/python.exe")
setwd("C:/Users/hp/Dev/playground/PO3")
library(MultiAssayExperiment)
library(MOFA2)
library(tidyr)
library(GEOquery)
library(survival)

mofa_trained <- readRDS("models/brca_mofa_trained.rds")
mae_obj <- readRDS("Data/brca_mae.rds")

fs <- get_factors(mofa_trained, as.data.frame = TRUE)
fw <- pivot_wider(fs, names_from = "factor", values_from = "value")
fm <- as.data.frame(fw)
rownames(fm) <- fm$sample
fm$sample <- NULL

rna_data <- as.matrix(assay(mae_obj, "RNA"))

correlations <- apply(rna_data, 1, function(g) cor(g, fm$Factor1, use = "complete.obs"))
top500 <- names(sort(abs(correlations), decreasing = TRUE)[1:500])

cat("Downloading GEO dataset GSE20685...\n")
gse <- getGEO("GSE20685", GSEMatrix = TRUE)[[1]]
geo_expr <- exprs(gse)
geo_clin <- pData(gse)

cat("GEO clinical columns:\n")
print(grep("survival|death|follow|os|days", colnames(geo_clin), value = TRUE, ignore.case = TRUE))

geo_surv <- data.frame(
    time = as.numeric(geo_clin$`survival time (months)`),
    status = as.numeric(geo_clin$`overall survival event (1=dead, 0=alive)`)
)
geo_surv <- geo_surv[complete.cases(geo_surv), ]
cat("GEO patients with survival data:", nrow(geo_surv), "\n")

overlap <- intersect(top500, rownames(geo_expr))
cat("Overlapping genes with TCGA top 500:", length(overlap), "\n")

geo_scores <- colMeans(geo_expr[overlap, , drop=FALSE], na.rm = TRUE)
geo_surv$score <- geo_scores[rownames(geo_surv)]
geo_surv <- geo_surv[complete.cases(geo_surv$score), ]

cox_geo <- tryCatch(
    coxph(Surv(time, status) ~ score, data = geo_surv),
    error = function(e) NULL
)

if (!is.null(cox_geo)) {
    s <- summary(cox_geo)
    cat("\nGEO validation:\n")
    cat("HR:", round(s$conf.int[1, 1], 3),
        "(", round(s$conf.int[1, 3], 3), "-", round(s$conf.int[1, 4], 3), ")",
        "p:", format(s$coefficients[1, 5], digits = 3), "\n")
} else {
    cat("GEO validation Cox model failed.\n")
}
cat("Done.\n")
```

**🧠 What you're doing:** Testing if the Factor 1 signature (top correlated genes) predicts survival in an independent breast cancer cohort (GSE20685).

**🎯 Why:** External validation is the gold standard — it proves your MOFA factors aren't just fitting noise in the TCGA data. A significant result in GEO means your findings are generalizable. Full standalone script at `scripts/step6.6_validation.R`.

---

## Troubleshooting Guide (BRCA-Specific)

| Problem | Solution |
| --- | --- |
| RPPA file not found | Check `Data/TCGA.BRCA.sampleMap/` — may need manual download. File should match `*RPPA*` |
| PAM50 subtype missing | Try `clinical_mat$PAM50Call_RNAseq` or `clinical_mat$PAM50_mRNA_nature2012`. BRCA clinicalMatrix uses `_nature2012` suffixed names |
| OS column name mismatch | Scripts updated: use `Vital_Status_nature2012`, `Days_to_date_of_Death_nature2012`, `Days_to_Date_of_Last_Contact_nature2012` |
| Age/stage column name mismatch | Scripts updated: use `Age_at_Initial_Pathologic_Diagnosis_nature2012` and `AJCC_Stage_nature2012` |
| Breast cancer is hormone-driven. Should I add ER/PR/HER2 status to Cox? | Yes — add `ER_Status_nature2012`, `PR_Status_nature2012`, `HER2_Final_Status_nature2012` as covariates if available |
| Recurrence data not found | BRCA in Xena may use `days_to_recurrence` or `days_to_tumor_progression`. Check column names with `grep("recur"` |
| Too many patients with no RPPA data | Filter to only patients with all 4 omics (Step 1.2). If < 400, consider dropping RPPA for more power |
| MOFA training slower than LUSC | 4 views + more patients = 2-3x slower. Allow 30-60 min. Use `convergence_mode = "fast"` |
| GEO validation fails | Try GSE20685 (largest BRCA with survival) manually first. Download from GEO website if `getGEO()` fails |
| Step 2.2 error: subscript out of bounds | Use `head(order(...), min(n, length(vars)))` instead of `names(sort(...))` |
| Step 2.3 error: colnames on object with <2 dims | `drop = FALSE` was missing in rowname deduplication lines. **FIXED** — all 4 now use `drop = FALSE` |
| create\_mofa error with empty views | **FIXED** — Step 3.2 now drops views with 0 features (e.g. CNV) before calling `create_mofa` |
| Step 3.3/4.1 file not found | **FIXED** — all references updated from `Data/brca_mofa.rds` to `models/brca_mofa_trained.rds` |