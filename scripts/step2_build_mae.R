setwd("C:/Users/hp/Dev/playground/PO3")
library(data.table)
library(stringr)
library(MultiAssayExperiment)
library(impute)

extract_patient_id <- function(barcodes) str_sub(barcodes, 1, 12)
deduplicate <- function(mat) {
    pids <- extract_patient_id(colnames(mat))
    dup <- duplicated(pids)
    if (sum(dup) > 0) cat("Removing", sum(dup), "duplicate samples\n")
    mat[, !dup, drop = FALSE]
}

# Load raw
rna_mat <- as.matrix(fread("Data/TCGA.BRCA.sampleMap/HiSeqV2_PANCAN.gz")[, -1])
rownames(rna_mat) <- fread("Data/TCGA.BRCA.sampleMap/HiSeqV2_PANCAN.gz")$sample
meth_mat <- as.matrix(fread("Data/TCGA.BRCA.sampleMap/HumanMethylation450.gz")[, -1])
rownames(meth_mat) <- fread("Data/TCGA.BRCA.sampleMap/HumanMethylation450.gz")$sample
cnv_mat <- as.matrix(fread("Data/TCGA.BRCA.sampleMap/Gistic2_CopyNumber_Gistic2_all_data_by_genes.gz")[, -1])
rownames(cnv_mat) <- fread("Data/TCGA.BRCA.sampleMap/Gistic2_CopyNumber_Gistic2_all_data_by_genes.gz")$sample
rppa_file <- list.files("Data/TCGA.BRCA.sampleMap/", pattern="RPPA", full.names=TRUE)
tmp <- fread(rppa_file[1])
rppa_mat <- as.matrix(tmp[, -1]); rownames(rppa_mat) <- tmp$sample
clinical <- fread("Data/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix")

# Common patients
common <- Reduce(intersect, lapply(list(rna_mat, meth_mat, cnv_mat, rppa_mat),
    function(m) unique(extract_patient_id(colnames(m)))))

# Filter & impute methylation
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

# Deduplicate + subset
rna_f <- deduplicate(rna_top)
meth_f <- deduplicate(meth_top)
cnv_f <- deduplicate(cnv_top)
rppa_f <- deduplicate(rppa_mat)

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
