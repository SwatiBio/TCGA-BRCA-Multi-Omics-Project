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
