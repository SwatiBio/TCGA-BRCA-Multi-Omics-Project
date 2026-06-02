Sys.setenv(RETICULATE_PYTHON = "C:/Users/hp/AppData/Local/Programs/Python/Python312/python.exe")
setwd("C:/Users/hp/Dev/playground/PO3")
library(MOFA2)
library(MultiAssayExperiment)

mofa_object <- readRDS("models/brca_mofa_trained.rds")
clinical <- readRDS("Data/brca_mae.rds")

# Extract factor scores
factor_scores <- get_factors(mofa_object, as.data.frame = TRUE)
head(factor_scores)

# Prepare clinical data
clinical_data <- colData(clinical) %>% as.data.frame()
colnames(clinical_data)

# Check for PAM50 subtypes
grep("PAM50|Subtype|BRCA_Subtype", colnames(clinical_data), value = TRUE, ignore.case = TRUE)

# Check for survival columns
grep("days_to_death|days_to_last_followup|vital_status|OS|recurrence|progression",
    colnames(clinical_data), value = TRUE, ignore.case = TRUE)

# --- Step 4.2: PAM50 Subtype Distribution ---
pam50_col <- grep("PAM50|Subtype|BRCA_Subtype", colnames(clinical_data), value = TRUE, ignore.case = TRUE)[1]
cat("PAM50 column:", pam50_col, "\n")

pam50 <- clinical_data[[pam50_col]]
print(table(pam50, useNA = "ifany"))

dir.create("figures", showWarnings = FALSE)

png("figures/pam50_subtype_distribution.png", width = 8, height = 6, units = "in", res = 300)
barplot(table(pam50),
    main = "PAM50 Subtype Distribution (TCGA-BRCA)",
    xlab = "PAM50 Subtype",
    ylab = "Number of Patients",
    las = 2,
    col = "steelblue")
dev.off()
cat("Saved: figures/pam50_subtype_distribution.png\n")
