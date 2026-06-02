Sys.setenv(RETICULATE_PYTHON = "C:/Users/hp/AppData/Local/Programs/Python/Python312/python.exe")
setwd("C:/Users/hp/Dev/playground/PO3")
library(stringr)
library(MultiAssayExperiment)
library(MOFA2)
library(survival)
library(survminer)

mae <- readRDS("Data/brca_mae.rds")
mofa_trained <- readRDS("models/brca_mofa_trained.rds")
clinical_mat <- as.data.frame(colData(mae))

factor_scores <- get_factors(mofa_trained)[[1]]
colnames(factor_scores) <- paste0("Factor", 1:ncol(factor_scores))

# === OS Cox ===
os_time <- as.numeric(ifelse(clinical_mat$Vital_Status_nature2012=="DECEASED", clinical_mat$Days_to_date_of_Death_nature2012, clinical_mat$Days_to_Date_of_Last_Contact_nature2012))
os_event <- ifelse(clinical_mat$Vital_Status_nature2012=="DECEASED", 1, 0)
cat("OS:", sum(os_event), "deaths,", sum(!os_event), "censored\n")

cox_os <- data.frame()
for (f in colnames(factor_scores)) {
    cox <- coxph(as.formula(paste0("Surv(os_time, os_event) ~ `", f, "`")),
                 data=data.frame(os_time, os_event, factor_scores))
    s <- summary(cox)
    cox_os <- rbind(cox_os, data.frame(Factor=f, Endpoint="OS",
        HR=exp(coef(cox)), p_value=coef(s)[,"Pr(>|z|)"],
        CI_lower=exp(confint(cox)[1]), CI_upper=exp(confint(cox)[2])))
}
cox_os$p_adj <- p.adjust(cox_os$p_value, method="BH")
write.csv(cox_os, "logs/cox_os.csv", row.names=FALSE)
print(cox_os)

# === RFS Cox ===
rfs_cols <- grep("recurrence|progression|recur|progres", names(clinical_mat), value=TRUE, ignore.case=TRUE)
cat("RFS columns:", paste(rfs_cols, collapse=", "), "\n")

rfs_time <- os_time; rfs_event <- os_event  # fallback = OS
for (col in rfs_cols) {
    v <- as.numeric(clinical_mat[[col]])
    if (sum(!is.na(v)) > 100 && median(v, na.rm=TRUE) < 5000) {
        rfs_time <- v
        cat("Using", col, "for RFS time\n"); break
    }
}
for (col in rfs_cols) {
    v <- suppressWarnings(as.numeric(clinical_mat[[col]]))
    if (sum(!is.na(v)) > 50 && all(v %in% c(0,1,NA))) {
        rfs_event <- v
        cat("Using", col, "for RFS event\n"); break
    }
    if (is.character(clinical_mat[[col]])) {
        v <- as.numeric(grepl("recur|progress|dead|deceased", clinical_mat[[col]], ignore.case=TRUE))
        rfs_event <- v; cat("Using", col, "for RFS (text parsed)\n"); break
    }
}

# Use OS time for patients with missing RFS time
na_rfs <- which(is.na(rfs_time)); rfs_time[na_rfs] <- os_time[na_rfs]; rfs_event[na_rfs] <- os_event[na_rfs]
cat("RFS:", sum(rfs_event, na.rm=TRUE), "events,", sum(rfs_event==0, na.rm=TRUE), "censored\n")

cox_rfs <- data.frame()
for (f in colnames(factor_scores)) {
    cox <- coxph(as.formula(paste0("Surv(rfs_time, rfs_event) ~ `", f, "`")),
                 data=data.frame(rfs_time, rfs_event, factor_scores))
    s <- summary(cox)
    cox_rfs <- rbind(cox_rfs, data.frame(Factor=f, Endpoint="RFS",
        HR=exp(coef(cox)), p_value=coef(s)[,"Pr(>|z|)"],
        CI_lower=exp(confint(cox)[1]), CI_upper=exp(confint(cox)[2])))
}
cox_rfs$p_adj <- p.adjust(cox_rfs$p_value, method="BH")
write.csv(cox_rfs, "logs/cox_rfs.csv", row.names=FALSE)
print(cox_rfs)

# === KM curves ===
for (label in c("OS","RFS")) {
    res <- if (label=="OS") cox_os else cox_rfs
    best <- res$Factor[which.min(res$p_value)]
    time <- if (label=="OS") os_time else rfs_time
    event <- if (label=="OS") os_event else rfs_event

    grp <- ifelse(factor_scores[,best] > median(factor_scores[,best], na.rm=TRUE), "High", "Low")
    km_data <- data.frame(time=time, event=event, group=grp)
    fit <- survfit(Surv(time, event) ~ group, data=km_data)
    p <- ggsurvplot(fit, data=km_data, pval=TRUE, risk.table=TRUE, title=paste(best, "(", label, ")"))
    dir.create("figures", showWarnings=FALSE)
    png(paste0("figures/km_", best, "_", label, ".png"), width=8, height=6, units="in", res=300)
    print(p); dev.off()
}
cat("Done.\n")
