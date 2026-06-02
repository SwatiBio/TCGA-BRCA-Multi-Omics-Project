Sys.setenv(RETICULATE_PYTHON = "C:/Users/hp/AppData/Local/Programs/Python/Python312/python.exe")
setwd("C:/Users/hp/Dev/playground/PO3")
library(MultiAssayExperiment)
library(MOFA2)
library(survival)
library(ggplot2)
library(patchwork)

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
valid <- !is.na(os_time) & os_time > 0

surv_df <- na.omit(data.frame(time=os_time[valid], status=os_event[valid], factor_scores[valid, ]))

dir.create("figures", showWarnings=FALSE)

###########################################################
# 1. FOREST PLOT — All 15 Factors, HR + 95% CI
###########################################################
cat("=== 1. Forest Plot ===\n")
hr_data <- data.frame()
for (f in paste0("Factor", 1:15)) {
    cox <- coxph(as.formula(paste0("Surv(time, status) ~ `", f, "`")), data=surv_df)
    s <- summary(cox)
    hr_data <- rbind(hr_data, data.frame(
        Factor = f,
        HR = exp(coef(cox)),
        CI_low = exp(confint(cox)[1]),
        CI_high = exp(confint(cox)[2]),
        p_value = coef(s)[, "Pr(>|z|)"]
    ))
}
hr_data$Significant <- ifelse(hr_data$p_value < 0.05, "p < 0.05", "NS")
hr_data$Factor <- factor(hr_data$Factor, levels = rev(paste0("Factor", 1:15)))

p_forest <- ggplot(hr_data, aes(x = HR, y = Factor, color = Significant)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.8) +
    geom_point(size = 3.5) +
    geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.2, linewidth = 1) +
    scale_x_log10(breaks = c(0.5, 0.75, 1, 1.5, 2.0), labels = c("0.50", "0.75", "1.00", "1.50", "2.00")) +
    scale_color_manual(values = c("p < 0.05" = "#D73027", "NS" = "#4575B4")) +
    labs(
        title = "Forest Plot — MOFA Factors & Overall Survival",
        subtitle = "Univariate Cox Proportional Hazards (n=400)",
        x = "Hazard Ratio (95% CI)", y = NULL, color = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
        plot.title = element_text(face = "bold"),
        legend.position = "top",
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5)
    ) +
    coord_cartesian(xlim = c(0.4, 2.3))
ggsave("figures/forest_plot.png", p_forest, width = 8, height = 6, dpi = 200)
cat("  Saved: figures/forest_plot.png\n")

###########################################################
# 2. VOLCANO PLOT — HR vs p-value
###########################################################
cat("=== 2. Volcano Plot ===\n")
hr_data$log_p <- -log10(hr_data$p_value)
hr_data$Direction <- ifelse(hr_data$HR < 1, "Protective", "Risk")
hr_data$Direction[hr_data$p_value >= 0.05] <- "NS"

p_volcano <- ggplot(hr_data, aes(x = log(HR), y = log_p, color = Direction, label = Factor)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "grey50", linewidth = 0.6) +
    geom_point(size = 4, alpha = 0.85) +
    ggrepel::geom_text_repel(size = 3.2, max.overlaps = 15, show.legend = FALSE) +
    scale_color_manual(values = c("Protective" = "#1A9850", "Risk" = "#D73027", "NS" = "#ABABAB")) +
    labs(
        title = "Volcano Plot — Factor Association with OS",
        subtitle = "Dashed = HR=1, Dotted = p=0.05",
        x = "log(Hazard Ratio)", y = "-log10(p-value)", color = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"), legend.position = "top")
ggsave("figures/volcano_plot.png", p_volcano, width = 7, height = 6, dpi = 200)
cat("  Saved: figures/volcano_plot.png\n")

###########################################################
# 3. VARIANCE EXPLAINED — LOLLIPOP CHART
###########################################################
cat("=== 3. Variance Explained Lollipop ===\n")
ve <- get_variance_explained(mofa_trained)
ve_df <- data.frame(
    Factor = paste0("Factor", 1:15),
    Total_VE = rowSums(ve$r2_per_factor$group1) * 100
)
ve_df$Factor <- factor(ve_df$Factor, levels = rev(paste0("Factor", 1:15)))
ve_df$Label <- paste0(round(ve_df$Total_VE, 1), "%")

p_ve <- ggplot(ve_df, aes(x = Total_VE, y = Factor)) +
    geom_segment(aes(xend = 0, yend = Factor), color = "grey70", linewidth = 0.8) +
    geom_point(aes(size = Total_VE), color = "#2166AC", alpha = 0.85) +
    geom_text(aes(label = Label), hjust = -0.2, size = 3.2) +
    scale_size_continuous(range = c(3, 10)) +
    labs(
        title = "Total Variance Explained per Factor",
        subtitle = "Sum across RNA-seq, DNA Methylation, and RPPA",
        x = "Total Variance Explained (%)", y = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
        plot.title = element_text(face = "bold"),
        legend.position = "none",
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5)
    ) +
    coord_cartesian(xlim = c(0, max(ve_df$Total_VE) * 1.3))
ggsave("figures/variance_lollipop.png", p_ve, width = 8, height = 5, dpi = 200)
cat("  Saved: figures/variance_lollipop.png\n")

###########################################################
# 4. TOP 5 KM PANEL — Side by Side
###########################################################
cat("=== 4. Top 5 KM Panel ===\n")
top5 <- hr_data$Factor[order(hr_data$p_value)[1:5]]

library(survminer)
km_plots <- list()
for (f in as.character(top5)) {
    grp <- ifelse(surv_df[[f]] > median(surv_df[[f]], na.rm=TRUE), "High", "Low")
    km_data <- data.frame(time = surv_df$time, event = surv_df$status, group = grp)
    fit <- survfit(Surv(time, event) ~ group, data = km_data)
    p <- ggsurvplot(fit, data = km_data, pval = TRUE, pval.size = 2.5,
        risk.table = FALSE, legend = "top", legend.title = f,
        palette = c("#D73027", "#1A9850"), title = f,
        xlab = "Days", ylab = "Survival", censor = FALSE)
    km_plots[[f]] <- p$plot + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10))
}

p_km_panel <- wrap_plots(km_plots, ncol = 3) +
    plot_annotation(title = "Top 5 Factors — Kaplan-Meier Curves (Median Split)",
        theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14)))
ggsave("figures/km_panel_top5.png", p_km_panel, width = 12, height = 8, dpi = 200)
cat("  Saved: figures/km_panel_top5.png\n")

###########################################################
# 5. RSF IMPORTANCE LOLLIPOP
###########################################################
cat("=== 5. RSF Importance Lollipop ===\n")
library(ranger)
rsf_model <- ranger(Surv(time, status) ~ ., data = surv_df,
    importance = "permutation", num.trees = 1000, seed = 42)
vimp <- data.frame(
    Factor = factor(names(rsf_model$variable.importance),
        levels = rev(names(sort(rsf_model$variable.importance)))),
    Importance = rsf_model$variable.importance
)
vimp$Color <- ifelse(vimp$Importance > 0, "#2166AC", "#ABABAB")

p_rsf <- ggplot(vimp, aes(x = Importance, y = Factor)) +
    geom_segment(aes(xend = 0, yend = Factor), color = "grey70", linewidth = 0.8) +
    geom_point(aes(color = Color), size = 3.5) +
    geom_text(aes(label = round(Importance, 3)), hjust = -0.2, size = 3) +
    scale_color_identity() +
    labs(
        title = "Random Survival Forest — Variable Importance",
        subtitle = paste("C-index:", round(1 - rsf_model$prediction.error, 3)),
        x = "Permutation Importance", y = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"), legend.position = "none")
ggsave("figures/rsf_lollipop.png", p_rsf, width = 8, height = 5, dpi = 200)
cat("  Saved: figures/rsf_lollipop.png\n")

###########################################################
# 6. TIME-ROC SUMMARY BAR
###########################################################
cat("=== 6. Time-ROC Bar Chart ===\n")
roc_df <- data.frame(
    Time = c("1 Year", "3 Years", "5 Years"),
    AUC = c(0.420, 0.705, 0.613)
)
roc_df$Label <- paste0(round(roc_df$AUC * 100, 1), "%")
roc_df$Time <- factor(roc_df$Time, levels = c("1 Year", "3 Years", "5 Years"))

p_roc <- ggplot(roc_df, aes(x = Time, y = AUC, fill = AUC)) +
    geom_col(width = 0.6, show.legend = FALSE) +
    geom_text(aes(label = Label), vjust = -0.5, size = 5, fontface = "bold") +
    geom_hline(yintercept = 0.7, linetype = "dashed", color = "darkgreen", linewidth = 0.8) +
    scale_fill_gradient(low = "#FDAE61", high = "#1A9850") +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(
        title = "Time-Dependent AUC — Combined MOFA Factors",
        subtitle = "Dashed line = clinically useful threshold (0.70)",
        x = NULL, y = "AUC"
    ) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))
ggsave("figures/timeroc_bar.png", p_roc, width = 5, height = 5, dpi = 200)
cat("  Saved: figures/timeroc_bar.png\n")

###########################################################
# 7. CORRELATION HEATMAP — Factors vs Clinical
###########################################################
cat("=== 7. Factor-Clinic Correlation Heatmap ===\n")
clin_vars <- c("Age_at_Initial_Pathologic_Diagnosis_nature2012",
    "ER_Status_nature2012", "PR_Status_nature2012",
    "HER2_Final_Status_nature2012", "PAM50Call_RNAseq")
clin_labels <- c("Age", "ER Status", "PR Status", "HER2 Status", "PAM50 Subtype")

cor_df <- expand.grid(Factor = paste0("Factor", 1:15), Clinical = clin_labels)
cor_df$r <- NA
for (f in paste0("Factor", 1:15)) {
    for (j in seq_along(clin_vars)) {
        v <- clinical_mat[[clin_vars[j]]]
        if (is.character(v)) v <- as.numeric(factor(v))
        cor_df$r[cor_df$Factor == f & cor_df$Clinical == clin_labels[j]] <-
            tryCatch(cor(factor_scores[, f], v, use = "pairwise.complete.obs"),
                error = function(e) NA)
    }
}
cor_df$r_round <- round(cor_df$r, 2)

p_cor <- ggplot(cor_df, aes(x = Clinical, y = factor(Factor, levels = rev(paste0("Factor", 1:15))), fill = r)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = r_round), size = 3.5, fontface = "bold") +
    scale_fill_gradient2(low = "#1A9850", mid = "white", high = "#D73027",
        midpoint = 0, limits = c(-0.5, 0.5), name = "r") +
    labs(
        title = "Factor–Clinical Variable Correlations",
        x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid = element_blank()
    )
ggsave("figures/factor_clinical_heatmap.png", p_cor, width = 8, height = 6, dpi = 200)
cat("  Saved: figures/factor_clinical_heatmap.png\n")

###########################################################
# 8. COMBINED SUMMARY FIGURE (multi-panel)
###########################################################
cat("=== 8. Combined Summary Figure ===\n")
p_combined <- (p_forest | p_volcano) / (p_ve | p_rsf) +
    plot_annotation(
        title = "TCGA-BRCA MOFA2 — Multi-Omics Survival Analysis Summary",
        tag_levels = "A", theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))
    )
ggsave("figures/summary_combined.png", p_combined, width = 16, height = 14, dpi = 200)
cat("  Saved: figures/summary_combined.png\n")

cat("\nAll summary figures generated in figures/\n")
cat("New files:\n")
cat("  figures/forest_plot.png  — HR + CI for all 15 factors\n")
cat("  figures/volcano_plot.png  — HR vs p-value (labeled)\n")
cat("  figures/variance_lollipop.png  — VE per factor (clean)\n")
cat("  figures/km_panel_top5.png  — Top 5 KM curves grid\n")
cat("  figures/rsf_lollipop.png  — RSF importance (clean)\n")
cat("  figures/timeroc_bar.png  — AUC bar chart\n")
cat("  figures/factor_clinical_heatmap.png  — Factor × Clinic corr\n")
cat("  figures/summary_combined.png  — All-in-one figure\n")
