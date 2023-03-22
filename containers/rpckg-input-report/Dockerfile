FROM rocker/tidyverse@sha256:0c62e932eeb3c5dd2529ed73274e670ffae169d7368073a59bc2282f28e99a6f
LABEL rpckg-input-report based on rocker/tidyverse 4.2.3

WORKDIR .

RUN mkdir ../opt/lib

RUN apt-get update \
 && apt-get upgrade -y \
 && apt-get install -y \
 git \
 git-lfs \
 make \
 gzip \
 rename \
 libglpk-dev \
 libfftw3-dev \
 libxt6 \
 gcc \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

COPY install_pandoc.sh ../opt/lib
RUN bash ../opt/lib/install_pandoc.sh

RUN /rocker_scripts/install_tidyverse.sh

COPY install_Rlibs_report.sh ../opt/lib
RUN bash ../opt/lib/install_Rlibs_report.sh

#for igraph
RUN Rscript -e 'install.packages(c("mvtnorm", "matrixcalc", "igraph", "gplots", "Matrix"))'

#packages from Bioconductor
COPY install_Bioconductorlibs_report.R ../opt/lib
RUN Rscript ../opt/lib/install_Bioconductorlibs_report.R

ENV SHELL=bash

WORKDIR /usr/src/report


