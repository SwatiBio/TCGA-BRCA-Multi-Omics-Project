Sys.setenv(RETICULATE_PYTHON = "C:/Users/hp/AppData/Local/Programs/Python/Python312/python.exe")
setwd("C:/Users/hp/Dev/playground/PO3")
library(MOFA2)
library(MultiAssayExperiment)

mae <- readRDS("Data/brca_mae.rds")

# Build MOFA data list
mofa_data <- list()
for (nm in names(experiments(mae))) {
    mat <- as.matrix(assay(experiments(mae)[[nm]]))
    dup <- duplicated(rownames(mat))
    if (any(dup)) mat <- mat[!dup, ]
    mofa_data[[nm]] <- mat
}
cat("Views:", names(mofa_data), "\nFeatures:", sapply(mofa_data, nrow), "\nSamples:", unique(sapply(mofa_data, ncol)), "\n")

# Drop empty views (e.g. CNV with 0 features)
mofa_data <- mofa_data[sapply(mofa_data, nrow) > 0]
cat("Keeping views:", names(mofa_data), "\n")

mofa_obj <- create_mofa(mofa_data)

model_opts <- get_default_model_options(mofa_obj)
model_opts$num_factors <- 15
model_opts$spikeslab_factors <- FALSE
model_opts$ard_factors <- TRUE

train_opts <- get_default_training_options(mofa_obj)
train_opts$convergence_mode <- "fast"
train_opts$maxiter <- 1000
train_opts$seed <- 42

mofa_obj <- prepare_mofa(mofa_obj,
    data_options = get_default_data_options(mofa_obj),
    model_options = model_opts,
    training_options = train_opts)

dir.create("models", showWarnings=FALSE)
mofa_trained <- run_mofa(mofa_obj, outfile="models/brca_mofa2.hdf5")
saveRDS(mofa_trained, "models/brca_mofa_trained.rds")
cat("MOFA training complete.\n")
