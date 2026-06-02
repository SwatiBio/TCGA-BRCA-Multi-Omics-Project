FROM rocker/shiny-verse:4.3.2

RUN apt-get update && apt-get install -y \
    libfontconfig1-dev \
    libcairo2-dev \
    libxt-dev \
    && rm -rf /var/lib/apt/lists/*

COPY install_packages.R /tmp/install_packages.R
RUN Rscript /tmp/install_packages.R

COPY shiny-server.conf /etc/shiny-server/shiny-server.conf

COPY dashboard_app /srv/shiny-server

EXPOSE 7860

CMD ["/usr/bin/shiny-server"]
