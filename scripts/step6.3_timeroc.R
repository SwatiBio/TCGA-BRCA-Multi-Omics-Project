Sys.setenv(RETICULATE_PYTHON = "C:/Users/hp/AppData/Local/Programs/Python/Python312/python.exe")
setwd("C:/Users/hp/Dev/playground/PO3")
library(MultiAssayExperiment)
library(MOFA2)
library(survival)
if (!requireNamespace("timeROC", quietly=TRUE)) install.packages("timeROC")
library(timeROC)

mae <- readRDS("Data/brca_mae.rds")
mofa_trained <- readRDS("models/brca_mofa_trained.rds")

factor_scores <- get_factors(mofa_trained)[[1]]
colnames(factor_scores) <- paste0("Factor", 1:ncol(factor_scores))

# Subset clinical to match factor_scores patients
clinical_mat <- as.data.frame(colData(mae))
common_ids <- intersect(rownames(clinical_mat), rownames(factor_scores))
clinical_mat <- clinical_mat[common_ids, , drop=FALSE]
factor_scores <- factor_scores[common_ids, , drop=FALSE]

os_time <- as.numeric(ifelse(clinical_mat$Vital_Status_nature2012=="DECEASED",
    clinical_mat$Days_to_date_of_Death_nature2012,
    clinical_mat$Days_to_Date_of_Last_Contact_nature2012))
os_event <- ifelse(clinical_mat$Vital_Status_nature2012=="DECEASED", 1, 0)

valid <- !is.na(os_time) & os_time > 0
surv_df <- na.omit(data.frame(time=os_time[valid], status=os_event[valid], factor_scores[valid, ]))
cat("Patients for ROC:", nrow(surv_df), "\n")

# Use Cox linear predictor as the marker (combines all factors)
cox_all <- coxph(Surv(time, status) ~ ., data=surv_df)
lp <- stats::predict(cox_all)

times <- c(365, 1095, 1825)

roc_result <- timeROC(
    T = surv_df$time,
    delta = surv_df$status,
    marker = lp,
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
cat("Done. AUC values printed above; figure saved.\n")
