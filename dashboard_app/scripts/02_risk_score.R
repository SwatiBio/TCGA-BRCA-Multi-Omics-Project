# 02_risk_score.R
# Builds the combined Cox risk score model, computes percentile rankings,
# and 3-year survival probabilities for every patient.
# Run this AFTER data_prep.R (needs data/patient_data.csv)
#
# Usage: source("dashboard_app/scripts/02_risk_score.R")

library(survival)
library(dplyr)

data_dir <- "C:/Users/hp/Dev/playground/PO3/dashboard_app/data"

patient_data <- read.csv(file.path(data_dir, "patient_data.csv"), row.names = 1)

# ---- Build combined Cox model using top factors + clinical variables ----
cox_combined <- coxph(
  Surv(os_time, os_event) ~ Factor5 + Factor14 + Factor8 + age + stage,
  data = patient_data
)

# ---- Compute risk score (linear predictor) ----
patient_data$risk_score <- predict(cox_combined, type = "lp")

# ---- Compute percentile ranking ----
patient_data$risk_percentile <- percent_rank(patient_data$risk_score) * 100

# ---- Compute 3-year survival probability ----
sf <- survfit(cox_combined, newdata = patient_data)
sf_sum <- summary(sf, times = 1095)
patient_data$surv_3yr <- if (is.matrix(sf_sum$surv)) sf_sum$surv[1, ] else sf_sum$surv

# ---- Save ----
write.csv(patient_data, file.path(data_dir, "patient_data.csv"))
saveRDS(cox_combined, file.path(data_dir, "cox_combined.rds"))

cat("Risk scores computed and saved.\n")
cat("  Patients:", nrow(patient_data), "\n")
cat("  Risk score range:", round(range(patient_data$risk_score, na.rm = TRUE), 3), "\n")
cat("  3-year survival range:", round(range(patient_data$surv_3yr, na.rm = TRUE), 3), "\n")
