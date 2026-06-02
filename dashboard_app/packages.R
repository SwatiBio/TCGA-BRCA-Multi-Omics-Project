pkgs <- c("bs4Dash", "shinyjs", "plotly", "DT", "survminer", "umap")
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org", dependencies = TRUE)
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
  BiocManager::install("MultiAssayExperiment", update = FALSE, ask = FALSE)
}
