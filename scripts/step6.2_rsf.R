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
cat("Done. Figures saved to figures/.\n")
