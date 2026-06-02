# ============================================================
# data_prep.R — Prepare all data for the interactive dashboard
# Run this ONCE in RStudio before launching the app
# ============================================================

results_dir <- "C:/Users/hp/Dev/playground/PO3/Results"
setwd("C:/Users/hp/Dev/playground/PO3")

library(dplyr)
library(survival)
library(MultiAssayExperiment)

cat("========== PREPARING DASHBOARD DATA ==========\n\n")

# 1. Load factor scores ---------------------------------
cat("1. Loading factor scores...\n")
factor_scores <- read.csv("Results/patient_factor_scores.csv", row.names = 1)
cat("   Loaded", nrow(factor_scores), "patients,", ncol(factor_scores), "factors\n")

# 2. Load clinical data from MAE ------------------------
cat("2. Loading clinical data from MAE...\n")
mae <- readRDS("Data/brca_mae.rds")
clinical <- as.data.frame(colData(mae))
cat("   Loaded", nrow(clinical), "clinical records\n")

# 3. Extract survival data ------------------------------
cat("3. Extracting survival data...\n")
os_time <- as.numeric(ifelse(
  clinical$Vital_Status_nature2012 == "DECEASED",
  clinical$Days_to_date_of_Death_nature2012,
  clinical$Days_to_Date_of_Last_Contact_nature2012
))
os_event <- ifelse(clinical$Vital_Status_nature2012 == "DECEASED", 1, 0)
cat("   Deaths:", sum(os_event), "Censored:", sum(!os_event), "\n")

# 4. Merge scores + survival ----------------------------
cat("4. Merging factor scores with survival data...\n")
common_ids <- intersect(rownames(factor_scores), rownames(clinical))
cat("   Common patients:", length(common_ids), "\n")

patient_data <- factor_scores[common_ids, , drop = FALSE]
patient_data$os_time <- os_time[match(common_ids, rownames(clinical))]
patient_data$os_event <- os_event[match(common_ids, rownames(clinical))]
patient_data <- na.omit(patient_data)
cat("   After removing NAs:", nrow(patient_data), "patients\n")

# Also get PAM50 subtypes for coloring
pam50 <- clinical$PAM50Call_RNAseq[match(rownames(patient_data), rownames(clinical))]
patient_data$subtype <- as.character(pam50)
patient_data$subtype[is.na(patient_data$subtype) | patient_data$subtype == ""] <- "Unknown"

# Get age and stage for risk score model
patient_data$age <- as.numeric(clinical$Age_at_Initial_Pathologic_Diagnosis_nature2012)[match(rownames(patient_data), rownames(clinical))]
patient_data$stage <- clinical$AJCC_Stage_nature2012[match(rownames(patient_data), rownames(clinical))]

# 5. Fit multivariate Cox model for simulator -----------
cat("5. Fitting multivariate Cox model...\n")
factor_cols <- paste0("`", colnames(factor_scores), "`", collapse = " + ")
cox_formula <- as.formula(paste("Surv(os_time, os_event) ~", factor_cols))
cox_multi <- coxph(cox_formula, data = patient_data)
cat("   Model converged:", !any(is.na(coef(cox_multi))), "\n")

# 6. Baseline survival ----------------------------------
cat("6. Computing baseline survival curve...\n")
baseline_surv <- survfit(cox_multi)
baseline_df <- data.frame(
  time = baseline_surv$time,
  surv = baseline_surv$surv
)
cat("   Baseline curve has", nrow(baseline_df), "time points\n")

# 7. Create subtype-specific KM data for display --------
cat("7. Computing overall KM curve for reference...\n")
ref_km <- survfit(Surv(os_time, os_event) ~ 1, data = patient_data)
ref_km_df <- data.frame(
  time = ref_km$time,
  surv = ref_km$surv,
  upper = ref_km$upper,
  lower = ref_km$lower
)

# 8. Factor descriptions for the mind map ---------------
cat("8. Creating factor descriptions...\n")
factor_desc <- data.frame(
  factor = paste0("Factor", 1:15),
  top_genes = c(
    "GABRP, C4orf7, STAC2", "PTPRT, TUSC5, ATP1A2",
    "C8orf42, C3orf30, HTR2C", "TMEM179, ZNF692, BEX4",
    "PLA2G2D, SIRPG, TIGIT, ICOS, ZAP70", "IL6ST, GAS2L3, DDI2",
    "PPAPDC1A, COL10A1, EPYC, FN1", "CFB, ADRA2C, TTC22",
    "TFAP2B, CPNE7, GGT6", "C19orf20, COL7A1, WNT9A",
    "TBC1D1, EIF2AK2, MGA", "RPL13AP3, PPIAL4C, ZNF205",
    "TEX28, OR2T10, SNORA42", "SPNS1, DMXL2, NACA2",
    "RAB43, MTRF1L, C11orf10"
  ),
  category = c(
    rep("Luminal Biology", 2), rep("Methylation", 2),
    "Immune", rep("Proliferation", 2),
    rep("Risk Factors", 3), rep("Mixed", 3),
    rep("Undefined", 2)
  )
)

# 9. Save everything ------------------------------------
cat("9. Saving to dashboard_app/data/...\n")
dir.create("dashboard_app/data", recursive = TRUE, showWarnings = FALSE)

write.csv(patient_data, "dashboard_app/data/patient_data.csv", row.names = TRUE)
saveRDS(cox_multi, "dashboard_app/data/cox_multi.rds")
write.csv(baseline_df, "dashboard_app/data/baseline_survival.csv", row.names = FALSE)
write.csv(ref_km_df, "dashboard_app/data/reference_km.csv", row.names = FALSE)
write.csv(factor_desc, "dashboard_app/data/factor_descriptions.csv", row.names = FALSE)

# Copy result CSVs (skip any that don't exist)
result_files <- c(
  "factor_variance_summary.csv", "factor_variance_explained.csv",
  "rmst_per_factor.csv", "timeroc_auc.csv",
  "brca_rsf_importance.csv", "rsf_cindex.csv",
  "top_features_per_factor.csv", "factor_clinical_corrs.csv",
  "total_variance_explained.csv"
)
for (f in result_files) {
  src <- file.path("Results", f)
  if (file.exists(src)) {
    file.copy(src, "dashboard_app/data/", overwrite = TRUE)
    cat("   Copied:", f, "\n")
  } else {
    cat("   Skipped (not found):", f, "\n")
  }
}

# 10. Summary -------------------------------------------
cat("\n========== DONE ==========\n")
cat("Files saved to: dashboard_app/data/\n")
cat("  - patient_data.csv (", nrow(patient_data), "patients × ", ncol(patient_data), "columns)\n")
cat("  - cox_multi.rds (multivariate Cox model)\n")
cat("  - baseline_survival.csv (baseline survival curve)\n")
cat("  - reference_km.csv (overall KM curve)\n")
cat("  - factor_descriptions.csv (descriptions for mind map)\n")
cat("  + 9 result CSVs from analysis\n")
cat("\nNext step: Open app.R and click 'Run App' in RStudio!\n")
