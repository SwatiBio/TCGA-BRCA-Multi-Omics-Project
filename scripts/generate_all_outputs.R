Sys.setenv(RETICULATE_PYTHON = "C:/Users/hp/AppData/Local/Programs/Python/Python312/python.exe")
setwd("C:/Users/hp/Dev/playground/PO3")
library(MultiAssayExperiment)
library(MOFA2)
library(survival)
library(ggplot2)

mae <- readRDS("Data/brca_mae.rds")
mofa_trained <- readRDS("models/brca_mofa_trained.rds")
clinical_mat <- as.data.frame(colData(mae))

factor_scores <- get_factors(mofa_trained)[[1]]
colnames(factor_scores) <- paste0("Factor", 1:ncol(factor_scores))

common_ids <- intersect(rownames(clinical_mat), rownames(factor_scores))
clinical_mat <- clinical_mat[common_ids, , drop=FALSE]
factor_scores <- factor_scores[common_ids, , drop=FALSE]

os_time <- as.numeric(ifelse(clinical_mat$Vital_Status_nature2012=="DECEASED",
    clinical_mat$Days_to_date_of_Death_nature2012,
    clinical_mat$Days_to_Date_of_Last_Contact_nature2012))
os_event <- ifelse(clinical_mat$Vital_Status_nature2012=="DECEASED", 1, 0)

dir.create("Results", showWarnings=FALSE)
dir.create("figures", showWarnings=FALSE)

# ============================================================
# 1. Variance Explained (Step 3.3)
# ============================================================
cat("=== Variance Explained ===\n")
ve <- get_variance_explained(mofa_trained)
ve_df <- as.data.frame(ve$r2_per_factor$group1)
ve_df$Factor <- paste0("Factor", 1:nrow(ve_df))
ve_df <- ve_df[, c("Factor", setdiff(names(ve_df), "Factor"))]
write.csv(ve_df, "Results/factor_variance_explained.csv", row.names=FALSE)

r2_total <- as.numeric(ve$r2_total$group1)
if (is.null(names(r2_total))) names(r2_total) <- colnames(ve_df)[-1]
ve_total <- data.frame(View = names(r2_total), Total_VE = round(r2_total, 4))
write.csv(ve_total, "Results/total_variance_explained.csv", row.names=FALSE)

png("figures/variance_explained.png", width=10, height=7, units="in", res=150)
print(plot_variance_explained(mofa_trained))
dev.off()

png("figures/variance_explained_total.png", width=8, height=5, units="in", res=150)
print(plot_variance_explained(mofa_trained, plot_total=TRUE))
dev.off()
cat("  Saved: figures/variance_explained.png, Results/factor_variance_explained.csv\n")

# ============================================================
# 2. Time-dependent ROC AUC (Step 6.3)
# ============================================================
cat("\n=== Time-dependent ROC AUC ===\n")
if (!requireNamespace("timeROC", quietly=TRUE)) install.packages("timeROC")
library(timeROC)

valid <- !is.na(os_time) & os_time > 0
surv_df <- na.omit(data.frame(time=os_time[valid], status=os_event[valid], factor_scores[valid, ]))
cat("  Patients for ROC:", nrow(surv_df), "\n")

cox_all <- coxph(Surv(time, status) ~ ., data=surv_df)
lp <- stats::predict(cox_all)
times <- c(365, 1095, 1825)

roc_result <- timeROC(
    T = surv_df$time, delta = surv_df$status,
    marker = lp, cause = 1, times = times, iid = TRUE
)

cat("  ROC SE structure:\n")
str(roc_result$se)
se_vals <- NULL
if (is.list(roc_result$se) && !is.null(roc_result$se$se)) {
    se_vals <- as.numeric(roc_result$se$se)
} else if (is.numeric(roc_result$se)) {
    se_vals <- as.numeric(roc_result$se)
} else if (is.list(roc_result$se) && length(roc_result$se) > 0) {
    se_vals <- as.numeric(roc_result$se[[1]])
}
if (length(se_vals) != length(times)) se_vals <- rep(NA_real_, length(times))
ci_lower <- roc_result$AUC - 1.96 * se_vals
ci_upper <- roc_result$AUC + 1.96 * se_vals
roc_df <- data.frame(
    Timepoint = c("1 Year", "3 Years", "5 Years"),
    Days = times,
    AUC = round(roc_result$AUC, 4),
    CI_lower = round(ci_lower, 4),
    CI_upper = round(ci_upper, 4)
)
write.csv(roc_df, "Results/timeroc_auc.csv", row.names=FALSE)
print(roc_df)
cat("  Saved: Results/timeroc_auc.csv\n")

# ============================================================
# 3. RMST per factor (Step 6.4)
# ============================================================
cat("\n=== RMST per Factor ===\n")
if (!requireNamespace("survRM2", quietly=TRUE)) install.packages("survRM2")
library(survRM2)

rmst_results <- data.frame()
for (f in paste0("Factor", 1:15)) {
    arm_high <- ifelse(surv_df[[f]] > median(surv_df[[f]]), 1, 0)
    rmst_f <- rmst2(surv_df$time, surv_df$status, arm=arm_high, tau=1825)
    res <- rmst_f$unadjusted.result[1, ]
    # Columns: Estimate, lower .95, upper .95, p
    rmst_results <- rbind(rmst_results, data.frame(
        Factor = f,
        RMST_Diff_Days = round(res[1], 1),
        CI_lower = round(res[2], 1),
        CI_upper = round(res[3], 1),
        P_value = round(res[4], 5)
    ))
    cat(sprintf("  %-10s diff=%6.1f days [%6.1f-%6.1f] p=%6.4f\n",
        f, res[1], res[2], res[3], res[4]))
}
write.csv(rmst_results, "Results/rmst_per_factor.csv", row.names=FALSE)
cat("  Saved: Results/rmst_per_factor.csv\n")

# ============================================================
# 4. RSF C-index (Step 6.2)
# ============================================================
cat("\n=== RSF C-index ===\n")
library(ranger)
surv_rf <- data.frame(time=os_time, status=os_event, factor_scores)
surv_rf <- surv_rf[complete.cases(surv_rf), ]

rsf_model <- ranger(
    Surv(time, status) ~ .,
    data = surv_rf,
    importance = "permutation",
    num.trees = 1000,
    seed = 42
)

c_index <- 1 - rsf_model$prediction.error
cat("  RSF C-index:", round(c_index, 4), "\n")
write.csv(data.frame(RSF_C_index = round(c_index, 4)), "Results/rsf_cindex.csv", row.names=FALSE)
cat("  Saved: Results/rsf_cindex.csv\n")

# ============================================================
# 5. Combined Summary Report
# ============================================================
cat("\n=== Generating Combined Summary Report ===\n")
sink("Results/analysis_summary.txt")
cat("TCGA-BRCA MOFA2 — Analysis Summary\n")
cat("====================================\n\n")
cat("Patients:", ncol(mae), "\n")
cat("Features:", sum(sapply(experiments(mae), nrow)), "\n")
cat("Omics:", paste(names(experiments(mae)), collapse=", "), "\n")
cat("MOFA Factors:", ncol(factor_scores), "\n\n")

cat("--- Variance Explained ---\n")
for (i in 1:nrow(ve_df)) {
    total <- sum(as.numeric(ve_df[i, -1]))
    cat(ve_df$Factor[i], ": ", paste0(round(as.numeric(ve_df[i, -1]) * 100, 1), "%", collapse=", "),
        " (total: ", round(total * 100, 1), "%)\n", sep="")
}

cat("\n--- Significant Factors (OS, p.adj<0.05) ---\n")
cox_res <- read.csv("logs/cox_os.csv")
sig <- cox_res[cox_res$p_adj < 0.05, ]
if (nrow(sig) > 0) {
    for (i in 1:nrow(sig)) {
        cat(sprintf("%s: HR=%.3f [%.3f-%.3f] p=%.4f p.adj=%.4f\n",
            sig$Factor[i], sig$HR[i], sig$CI_lower[i], sig$CI_upper[i],
            sig$p_value[i], sig$p_adj[i]))
    }
} else {
    cat("  None at p.adj<0.05\n")
}

cat("\n--- Time-dependent AUC (all factors combined) ---\n")
for (i in 1:nrow(roc_df)) {
    cat(sprintf("%s: AUC=%.4f [%.4f-%.4f]\n",
        roc_df$Timepoint[i], roc_df$AUC[i], roc_df$CI_lower[i], roc_df$CI_upper[i]))
}

cat("\n--- Top 3 RSF Variable Importance ---\n")
rsf_vimp <- read.csv("Results/brca_rsf_importance.csv")
for (i in 1:min(3, nrow(rsf_vimp))) {
    cat(rsf_vimp$Factor[i], ": ", round(rsf_vimp$Importance[i], 4), "\n", sep="")
}

cat("\n--- RMST (5-year, per factor) ---\n")
for (i in 1:nrow(rmst_results)) {
    if (rmst_results$P_value[i] < 0.05) {
        cat(sprintf("%s: diff=%.1f days [%.1f-%.1f] p=%.4f *\n",
            rmst_results$Factor[i], rmst_results$RMST_Diff_Days[i],
            rmst_results$CI_lower[i], rmst_results$CI_upper[i],
            rmst_results$P_value[i]))
    }
}
sink()
cat("  Saved: Results/analysis_summary.txt\n")
cat("\nAll outputs generated successfully.\n")
