Sys.setenv(RETICULATE_PYTHON = "C:/Users/hp/AppData/Local/Programs/Python/Python312/python.exe")
setwd("C:/Users/hp/Dev/playground/PO3")
library(MultiAssayExperiment)
library(MOFA2)
library(survival)
if (!requireNamespace("survRM2", quietly=TRUE)) install.packages("survRM2")
library(survRM2)

mae <- readRDS("Data/brca_mae.rds")
mofa_trained <- readRDS("models/brca_mofa_trained.rds")

factor_scores <- get_factors(mofa_trained)[[1]]
colnames(factor_scores) <- paste0("Factor", 1:ncol(factor_scores))

# Align clinical with factor_scores patients
common_ids <- intersect(rownames(colData(mae)), rownames(factor_scores))
clinical_mat <- as.data.frame(colData(mae))[common_ids, , drop=FALSE]
factor_scores <- factor_scores[common_ids, , drop=FALSE]

os_time <- as.numeric(ifelse(clinical_mat$Vital_Status_nature2012=="DECEASED",
    clinical_mat$Days_to_date_of_Death_nature2012,
    clinical_mat$Days_to_Date_of_Last_Contact_nature2012))
os_event <- ifelse(clinical_mat$Vital_Status_nature2012=="DECEASED", 1, 0)

valid <- !is.na(os_time) & os_time > 0
surv_df <- na.omit(data.frame(time=os_time[valid], status=os_event[valid], factor_scores[valid, ]))
cat("Patients for RMST:", nrow(surv_df), "\n")

cat("Per-factor RMST (5-year restricted mean survival time):\n")
for (f in paste0("Factor", 1:15)) {
    arm_high <- ifelse(surv_df[[f]] > median(surv_df[[f]]), 1, 0)
    rmst_f <- rmst2(surv_df$time, surv_df$status, arm = arm_high, tau = 1825)
    res <- rmst_f$unadjusted.result[1, ]
    # Columns: Estimate, lower .95, upper .95, p
    cat(sprintf("%-10s diff=%6.1f days [%6.1f-%6.1f] p=%6.4f\n",
        f, res[1], res[2], res[3], res[4]))
}
cat("Done.\n")
