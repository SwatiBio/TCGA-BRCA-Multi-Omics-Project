FROM rocker/shiny:4.3.2

RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libfontconfig1-dev \
    libcairo2-dev \
    libxt-dev \
    && rm -rf /var/lib/apt/lists/*

# Install dev httpuv for WebSocket reliability
RUN install2.r --error --skipinstalled remotes && installGithub.r rstudio/httpuv

# Install CRAN packages (no --skipinstalled: force-consistency with pre-installed pkgs)
RUN install2.r --error \
    shiny \
    bs4Dash \
    shinyjs \
    plotly \
    DT \
    survminer \
    umap \
    dplyr \
    survival \
    ggplot2

# Verify all packages load correctly
RUN Rscript -e "for(pkg in c('shiny','bs4Dash','shinyjs','plotly','DT','dplyr','survival','survminer','ggplot2','umap')) { cat('Loading',pkg,'...'); suppressPackageStartupMessages(library(pkg,character.only=TRUE)); cat('OK\n') }" && echo "All packages loaded successfully"

# Install Bioc packages
RUN Rscript -e "install.packages('BiocManager',repos='https://cloud.r-project.org');BiocManager::install('MultiAssayExperiment',update=FALSE,ask=FALSE)"

COPY dashboard_app /srv/shiny-server
COPY shiny-server.conf /etc/shiny-server/shiny-server.conf

RUN chown -R shiny:shiny /srv/shiny-server

# Inject debug logging into server function (minimal, reversible change in COPY only)
RUN Rscript -e "\
  src <- readLines('/srv/shiny-server/app.R'); \
  server_start <- grep(\"^server <- function\", src); \
  cat('=== Server function at line', server_start, '===\n'); \
  debug_lines <- c( \
    '  try(cat(format(Sys.time()), \"- SERVER STARTED\\n\", file = \"/srv/shiny-server/www/debug.txt\"))', \
    '  options(shiny.error = function() {', \
    '  try(cat(format(Sys.time()), \"- REACTIVE ERROR:\", geterrmessage(), \"\\n\",', \
    '        file = \"/srv/shiny-server/www/debug.txt\", append = TRUE))', \
    '  })' \
  ); \
  src <- append(src, debug_lines, after = server_start); \
  writeLines(src, '/srv/shiny-server/app.R'); \
  cat('=== DEBUG LINES INJECTED ===\n'); \
  cat('Resulting server line:', grep(\"^server <- function\", src, value=TRUE), '\n')"

RUN Rscript -e "\
  setwd('/srv/shiny-server'); \
  cat('CWD:', getwd(), '\n'); \
  cat('data dir exists:', dir.exists('data'), '\n'); \
  cat('data files:', paste(dir('data'), collapse=', '), '\n'); \
  options(warn=2); \
  tryCatch({ \
    source('app.R', local=TRUE); \
    cat('=== APP SOURCED SUCCESSFULLY ===\n'); \
    obj <- try(shinyApp(ui, server)); \
    cat('shinyApp class:', paste(class(obj), collapse=', '), '\n'); \
  }, error = function(e) { \
    cat('!!! SOURCE ERROR:', conditionMessage(e), '\n'); \
    q(status=1, save='no'); \
  })"

EXPOSE 7860

RUN Rscript -e "sink('/srv/shiny-server/www/build_info.txt'); cat('BUILD TIME:', date(), '\n'); cat('R version:', R.version.string, '\n'); cat('Packages:\n'); for(pkg in c('shiny','bs4Dash','shinyjs','plotly','DT','dplyr','survival','survminer','ggplot2','umap')) cat(' ', pkg, '-', as.character(packageVersion(pkg)), '\n'); sink()"

CMD ["/usr/bin/shiny-server"]
