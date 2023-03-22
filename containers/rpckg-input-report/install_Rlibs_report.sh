#!/bin/bash

# script to install R packages required for generation of reports by CRISPR-pooled-RSL pipeline
# based on intstall_tidyverse.sh from https://github.com/rocker-org/rocker-versioned2
# used in images based on rocker/tidyverse:4.2.3

set -e

## build ARGs
NCPUS=${NCPUS:--1}

# a function to install apt packages only if they are not installed
function apt_install() {
    if ! dpkg -s "$@" >/dev/null 2>&1; then
        if [ "$(find /var/lib/apt/lists/* | wc -l)" = "0" ]; then
            apt-get update
        fi
        apt-get install -y --no-install-recommends "$@"
    fi
}


#packages for report from CRAN
install2.r --error --skipinstalled -n "$NCPUS" \
    DT \
    plotly\
    viridis \
    knitr \
    ggrepel \
    reshape2 \
    htmltools \
    magrittr \
    MASS \
    cowplot \
    bookdown \
    ggnewscale \
    kableExtra

# Clean up
#rm -rf /var/lib/apt/lists/*
rm -rf /tmp/downloaded_packages

## Strip binary installed lybraries from RSPM
## https://github.com/rocker-org/rocker-versioned2/issues/340
strip /usr/local/lib/R/site-library/*/libs/*.so

# Check the tidyverse core packages' version
echo -e "Check the CRAN packages...\n"

R -q -e "library(DT)"
R -q -e "library(plotly)"
R -q -e "library(viridis)"
R -q -e "library(kableExtra)"

echo -e "\nInstall CRAN packages, done!"

