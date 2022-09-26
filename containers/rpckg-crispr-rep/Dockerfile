FROM rocker/tidyverse:4.2.1
LABEL rpckg-crispr-rep based on rocker/tidyverse

WORKDIR .

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
 gcc \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

COPY install_pandoc.sh .
RUN bash install_pandoc.sh

RUN /rocker_scripts/install_tidyverse.sh

COPY install_Rlibs_report.sh .
RUN bash install_Rlibs_report.sh

#for igraph
RUN Rscript -e 'install.packages(c("mvtnorm", "matrixcalc", "igraph", "gplots", "Matrix"))'

#packages from Bioconductor
COPY install_Bioconductorlibs_report.R .
RUN Rscript install_Bioconductorlibs_report.R
