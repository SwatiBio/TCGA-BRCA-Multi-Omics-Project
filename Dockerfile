FROM rocker/shiny:4.3.2

RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libfontconfig1-dev \
    libcairo2-dev \
    libxt-dev \
    && rm -rf /var/lib/apt/lists/*

# Install dev httpuv for WebSocket reliability
RUN Rscript -e "install.packages('remotes', repos = 'https://cloud.r-project.org')" \
    && Rscript -e "remotes::install_github('rstudio/httpuv', upgrade = 'never')"

COPY install_packages.R /tmp/install_packages.R
RUN Rscript /tmp/install_packages.R

COPY dashboard_app /srv/shiny-server
COPY shiny-server.conf /etc/shiny-server/shiny-server.conf

RUN chown -R shiny:shiny /srv/shiny-server

EXPOSE 7860

CMD ["/usr/bin/shiny-server"]
