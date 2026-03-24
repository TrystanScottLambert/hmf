FROM rocker/r-ver:4.3.2

# System dependencies needed by some R packages
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libcairo2-dev \
    && apt-get clean all \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages(c( \
    'Cairo', \
    'data.table', \
    'rstan' \
    ), repos='https://cloud.r-project.org')"
RUN R -e 'install.packages("arrow", repos = "https://packagemanager.rstudio.com/all/__linux__/focal/latest")'

# celestial is on CRAN but sometimes needs explicit install
RUN R -e "install.packages('celestial', repos='https://cloud.r-project.org')"

# Copy your script into the image
COPY testing_the_affect_of_a.R /scripts/nessie_scaling.R

WORKDIR /scripts
