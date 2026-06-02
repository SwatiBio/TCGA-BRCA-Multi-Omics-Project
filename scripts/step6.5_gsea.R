Sys.setenv(RETICULATE_PYTHON = "C:/Users/hp/AppData/Local/Programs/Python/Python312/python.exe")
setwd("C:/Users/hp/Dev/playground/PO3")
library(MultiAssayExperiment)
library(MOFA2)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)

mofa_trained <- readRDS("models/brca_mofa_trained.rds")
mae_obj <- readRDS("Data/brca_mae.rds")

fs <- get_factors(mofa_trained, as.data.frame = TRUE)
fw <- pivot_wider(fs, names_from = "factor", values_from = "value")
fm <- as.data.frame(fw)
rownames(fm) <- fm$sample
fm$sample <- NULL

rna_data <- as.matrix(assay(mae_obj, "RNA"))

cat("Correlating Factor 1 with all genes...\n")
correlations <- apply(rna_data, 1, function(g) cor(g, fm$Factor1, use = "complete.obs"))
top_genes <- names(sort(abs(correlations), decreasing = TRUE)[1:500])

gene_df <- bitr(top_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

go_result <- enrichGO(
    gene = gene_df$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
)

dir.create("Results", showWarnings=FALSE)
dir.create("figures", showWarnings=FALSE)

if (!is.null(go_result) && nrow(go_result) > 0) {
    png("figures/brca_factor1_go_enrichment.png", width = 10, height = 8, units = "in", res = 150)
    print(dotplot(go_result, showCategory = 15))
    dev.off()
    write.csv(as.data.frame(go_result), "Results/brca_factor1_go_enrichment.csv", row.names = FALSE)
    cat("GSEA complete —", nrow(go_result), "enriched terms found.\n")
} else {
    cat("No enriched terms found at p<0.05.\n")
}
