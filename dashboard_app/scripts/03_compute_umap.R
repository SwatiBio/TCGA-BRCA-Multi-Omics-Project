# 03_compute_umap.R
# Computes UMAP coordinates from the 15 MOFA factor scores
# and adds them to patient_data.csv for the 2D patient map.
# Run this AFTER data_prep.R (needs data/patient_data.csv)
#
# Usage: source("dashboard_app/scripts/03_compute_umap.R")

library(umap)

data_dir <- "C:/Users/hp/Dev/playground/PO3/dashboard_app/data"
results_dir <- "C:/Users/hp/Dev/playground/PO3/Results"

# Load patient factor scores
patient_data <- read.csv(file.path(data_dir, "patient_data.csv"), row.names = 1)
factor_scores <- read.csv(file.path(results_dir, "patient_factor_scores.csv"), row.names = 1)

# Find columns matching Factor pattern
factor_cols <- grep("^Factor[0-9]", names(factor_scores), value = TRUE)
cat("Using", length(factor_cols), "factors for UMAP.\n")

# Run UMAP (15D -> 2D)
umap_result <- umap(as.matrix(factor_scores[, factor_cols]))

# Add coordinates to patient_data (subset to match patient_data rows)
layout <- umap_result$layout
patient_data$umap_x <- layout[rownames(patient_data), 1]
patient_data$umap_y <- layout[rownames(patient_data), 2]

# Save
write.csv(patient_data, file.path(data_dir, "patient_data.csv"))

cat("UMAP coordinates computed and saved.\n")
cat("  Patients in patient_data:", nrow(patient_data), "\n")
cat("  UMAP1 range:", round(range(patient_data$umap_x), 3), "\n")
cat("  UMAP2 range:", round(range(patient_data$umap_y), 3), "\n")
