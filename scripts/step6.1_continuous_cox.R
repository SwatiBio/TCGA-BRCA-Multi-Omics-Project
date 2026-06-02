Sys.setenv(RETICULATE_PYTHON = "C:/Users/hp/AppData/Local/Programs/Python/Python312/python.exe")
setwd("C:/Users/hp/Dev/playground/PO3")
library(stringr)
library(MultiAssayExperiment)
library(MOFA2)
library(survival)
library(glmnet)
library(ggplot2)
library(survminer)

mae <- readRDS("Data/brca_mae.rds")
mofa_trained <- readRDS("models/brca_mofa_trained.rds")
clinical_mat <- as.data.frame(colData(mae))

factor_scores <- get_factors(mofa_trained)[[1]]
colnames(factor_scores) <- paste0("Factor", 1:ncol(factor_scores))

os_time <- as.numeric(ifelse(clinical_mat$Vital_Status_nature2012=="DECEASED", clinical_mat$Days_to_date_of_Death_nature2012, clinical_mat$Days_to_Date_of_Last_Contact_nature2012))
os_event <- ifelse(clinical_mat$Vital_Status_nature2012=="DECEASED", 1, 0)
age <- as.numeric(clinical_mat$Age_at_Initial_Pathologic_Diagnosis_nature2012)
stage <- ifelse(grepl("III|IV", clinical_mat$AJCC_Stage_nature2012), "Late", "Early")
cat("OS:", sum(os_event), "deaths,", sum(os_event==0), "censored\n")

# Remove NAs and non-positive times (glmnet requires positive times)
valid <- !is.na(os_time) & os_time > 0
fs <- factor_scores[valid, ]
ot <- os_time[valid]
oe <- os_event[valid]
cat("Valid patients:", sum(valid), "\n")

# A: LASSO
set.seed(42)
cv <- cv.glmnet(scale(fs), Surv(ot, oe), family="cox", alpha=1, nfolds=10, cox.ties="breslow")
sel <- names(which(as.matrix(coef(cv, s="lambda.1se")) != 0))
if (length(sel)==0) sel <- names(which(as.matrix(coef(cv, s="lambda.min")) != 0))
cat("LASSO:", if(length(sel)>0) paste(sel, collapse=", ") else "none\n")

# B: Continuous univariate + C: Adjusted
results <- data.frame()
results_adj <- data.frame()
clin <- data.frame(ot, oe, age=scale(age[valid]), stage=factor(stage[valid]), scale(fs))
clin <- clin[complete.cases(clin), ]

for (f in paste0("Factor", 1:15)) {
    # Unadjusted
    cox <- coxph(as.formula(paste0("Surv(ot, oe) ~ `", f, "`")),
                 data=data.frame(ot, oe, scale(fs)))
    s <- summary(cox)
    results <- rbind(results, data.frame(Factor=f, HR=exp(coef(cox)), p_value=coef(s)[,"Pr(>|z|)"],
        CI_lower=exp(confint(cox)[1]), CI_upper=exp(confint(cox)[2])))
    # Adjusted
    cox2 <- tryCatch(coxph(as.formula(paste0("Surv(ot, oe) ~ `", f, "` + age + stage")), data=clin), error=function(e) NULL)
    if (!is.null(cox2)) {
        s2 <- summary(cox2)
        results_adj <- rbind(results_adj, data.frame(Factor=f, HR=exp(coef(cox2))[1], p_value=coef(s2)[1,"Pr(>|z|)"],
            CI_lower=exp(confint(cox2))[f,1], CI_upper=exp(confint(cox2))[f,2]))
    }
}
results$p_adj <- p.adjust(results$p_value, method="BH")
results_adj$p_adj <- p.adjust(results_adj$p_value, method="BH")

cat("\nUnadjusted:\n"); for (i in 1:nrow(results)) cat(sprintf("%-10s HR=%6.3f [%5.3f-%5.3f] p=%6.4f adj=%6.4f\n", results$Factor[i], results$HR[i], results$CI_lower[i], results$CI_upper[i], results$p_value[i], results$p_adj[i]))
cat("\nAdjusted for age+stage:\n"); for (i in 1:nrow(results_adj)) cat(sprintf("%-10s HR=%6.3f [%5.3f-%5.3f] p=%6.4f adj=%6.4f\n", results_adj$Factor[i], results_adj$HR[i], results_adj$CI_lower[i], results_adj$CI_upper[i], results_adj$p_value[i], results_adj$p_adj[i]))

write.csv(results, "logs/cox_continuous_os.csv", row.names=FALSE)
write.csv(results_adj, "logs/cox_adjusted_clinical_os.csv", row.names=FALSE)

# D: KM curves (median + tertile for top 3)
top3 <- results$Factor[order(results$p_value)][1:3]
for (f in top3) {
    med <- median(fs[,f], na.rm=TRUE)
    grp <- ifelse(fs[,f] > med, "High", "Low")
    kmd <- data.frame(time=ot, event=oe, group=grp)
    fit <- survfit(Surv(time, event) ~ group, data=kmd)
    p <- ggsurvplot(fit, data=kmd, pval=TRUE, risk.table=TRUE, title=paste(f, "- median (p=", round(results$p_value[results$Factor==f],4), ")"))
    png(paste0("figures/km_median_", f, ".png"), width=8, height=6, units="in", res=300); print(p); dev.off()

    q <- quantile(fs[,f], c(1/3, 2/3), na.rm=TRUE)
    grp <- factor(ifelse(fs[,f] > q[2], "High", ifelse(fs[,f] <= q[1], "Low", "Mid")), levels=c("Low","Mid","High"))
    idx <- grp != "Mid"
    kmd <- data.frame(time=ot[idx], event=oe[idx], group=droplevels(grp[idx]))
    fit <- survfit(Surv(time, event) ~ group, data=kmd)
    p <- ggsurvplot(fit, data=kmd, pval=TRUE, risk.table=TRUE, title=paste(f, "- extreme tertiles"))
    png(paste0("figures/km_tertile_", f, ".png"), width=8, height=6, units="in", res=300); print(p); dev.off()
}
cat("Saved figures.\n")
