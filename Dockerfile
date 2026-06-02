FROM rocker/shiny-verse:4.3.2

RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libfontconfig1-dev \
    libcairo2-dev \
    libxt-dev \
    && rm -rf /var/lib/apt/lists/*

COPY dashboard_app /srv/shiny-server

RUN Rscript -e "
options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/__linux__/jammy/latest'), install.packages.check.source = 'no');

pkgs <- c('bs4Dash', 'shinyjs', 'plotly', 'DT', 'survminer', 'umap');
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

if (!requireNamespace('BiocManager', quietly = TRUE)) {
  install.packages('BiocManager', repos = 'https://cloud.r-project.org')
}
if (!requireNamespace('MultiAssayExperiment', quietly = TRUE)) {
  BiocManager::install('MultiAssayExperiment', update = FALSE, ask = FALSE)
}
"

EXPOSE 7860

CMD ["R", "-e", "shiny::runApp('/srv/shiny-server', host='0.0.0.0', port=7860)"]
