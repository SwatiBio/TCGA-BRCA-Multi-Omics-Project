Sys.setenv(RETICULATE_PYTHON = "C:/Users/hp/AppData/Local/Programs/Python/Python312/python.exe")
setwd("C:/Users/hp/Dev/playground/PO3")
library(MultiAssayExperiment)
library(MOFA2)
library(tidyr)
if (!requireNamespace("GEOquery", quietly=TRUE)) {
    BiocManager::install("GEOquery")
}
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
# Align RNA data columns to factor score rows
common_patients <- intersect(colnames(rna_data), rownames(fm))
rna_data <- rna_data[, common_patients, drop=FALSE]
fm <- fm[common_patients, , drop=FALSE]

correlations <- apply(rna_data, 1, function(g) cor(g, fm$Factor1, use = "complete.obs"))
top500 <- names(sort(abs(correlations), decreasing = TRUE)[1:500])

cat("Downloading GEO dataset GSE20685...\n")
dir.create("Data/GEO", showWarnings=FALSE, recursive=TRUE)
gzfile <- "Data/GEO/GSE20685_series_matrix.txt.gz"

if (!file.exists(gzfile) || file.info(gzfile)$size < 100e6) {
    url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE20nnn/GSE20685/matrix/GSE20685_series_matrix.txt.gz"
    options(timeout=600)
    
    dl_ok <- FALSE
    for (m in c("libcurl", "auto", "wget", "curl")) {
        tryCatch({
            download.file(url, gzfile, mode="wb", method=m)
            if (file.exists(gzfile) && file.info(gzfile)$size > 100e6) { dl_ok <- TRUE; break }
        }, error=function(e) cat("Method", m, "failed, trying next...\n"))
    }
    
    if (!dl_ok) {
        msg <- paste0(
            "Could not download GSE20685 automatically.\n",
            "Please manually download:\n",
            "  ", url, "\n",
            "  -> ", normalizePath("Data/GEO"), "/GSE20685_series_matrix.txt.gz\n",
            "Then re-run this script.\n"
        )
        stop(msg)
    }
} else {
    cat("Using cached file:", gzfile, "\n")
}
gse <- getGEO(filename=gzfile, AnnotGPL=TRUE)
geo_expr <- exprs(gse)
geo_clin <- pData(gse)

cat("GEO clinical columns:\n")
print(colnames(geo_clin)[grep("event|death|survival|follow|time", colnames(geo_clin), ignore.case=TRUE)])
cat("Sample clinical row:\n")
print(sapply(geo_clin[1, grep("event|death|survival|follow|time", colnames(geo_clin), ignore.case=TRUE)], as.character))

# Map survival columns for GSE20685 (breast cancer)
geo_surv <- data.frame(
    time = as.numeric(geo_clin$`follow_up_duration (years):ch1`) * 365.25,  # years -> days
    status = as.numeric(ifelse(geo_clin$`event_death:ch1` == "1", 1, 0))
)
geo_surv <- geo_surv[complete.cases(geo_surv) & geo_surv$time > 0, ]
cat("GEO patients with survival data:", nrow(geo_surv), "\n")

# Map probe IDs to gene symbols using platform annotation
cat("Mapping probes to genes...\n")
gpl <- getGEO(annotation(gse), destdir="Data/GEO", AnnotGPL=TRUE)
gpl_tab <- Table(gpl)
cat("GPL columns:", paste(colnames(gpl_tab), collapse=", "), "\n")
# Find the gene symbol column (varies by platform)
symbol_col <- grep("symbol|gene|gene_symbol", colnames(gpl_tab), value=TRUE, ignore.case=TRUE)[1]
if (is.na(symbol_col)) symbol_col <- colnames(gpl_tab)[grep("Symbol|GENE", colnames(gpl_tab))[1]]
cat("Using column:", symbol_col, "\n")
# Clean gene symbols: split on ///, trim whitespace, take first symbol
clean_symbols <- function(x) {
    x <- trimws(gsub("///.*", "", x))
    x[x %in% c("", "---", "NA")] <- NA
    x
}
probe2gene <- data.frame(
    probe = gpl_tab$ID,
    symbol = clean_symbols(gpl_tab[[symbol_col]]),
    stringsAsFactors=FALSE
)
probe2gene <- probe2gene[!is.na(probe2gene$symbol), ]

# Aggregate expression from probes to genes (take max probe per gene)
common_probes <- intersect(probe2gene$probe, rownames(geo_expr))
cat("Probes matching annotation:", length(common_probes), "\n")
geo_expr_genes <- t(sapply(split(common_probes, probe2gene$symbol[match(common_probes, probe2gene$probe)]),
    function(probes) apply(geo_expr[probes, , drop=FALSE], 2, max, na.rm=TRUE)))

overlap <- intersect(top500, rownames(geo_expr_genes))
cat("Overlapping genes with TCGA top 500:", length(overlap), "\n")

cox_geo <- NULL
if (length(overlap) >= 10) {
    geo_scores <- colMeans(geo_expr_genes[overlap, , drop=FALSE], na.rm=TRUE)
    geo_surv$score <- geo_scores[match(rownames(geo_surv), names(geo_scores))]
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
} else {
    cat("Too few overlapping genes (<10) for validation.\n")
}

# Save results
dir.create("Results", showWarnings=FALSE)
sink("Results/geo_validation.txt")
cat("GEO Validation Results\n")
cat("======================\n\n")
cat("Dataset: GSE20685 (breast cancer)\n")
cat("Platform:", annotation(gse), "\n")
cat("Probes in annotation:", nrow(probe2gene), "\n")
cat("Unique genes mapped:", length(unique(probe2gene$symbol)), "\n")
cat("Top 500 TCGA genes used\n")
cat("Overlap with GEO:", length(overlap), "\n")
if (!is.null(cox_geo) && length(overlap) >= 10) {
    s <- summary(cox_geo)
    cat("\nCox Regression:\n")
    cat("HR:", round(s$conf.int[1, 1], 3),
        "(", round(s$conf.int[1, 3], 3), "-", round(s$conf.int[1, 4], 3), ")",
        "\np:", format(s$coefficients[1, 5], digits=3), "\n")
    cat("N patients:", nrow(geo_surv), "\n")
}
sink()
cat("\nResults saved to Results/geo_validation.txt\n")
cat("Done.\n")
