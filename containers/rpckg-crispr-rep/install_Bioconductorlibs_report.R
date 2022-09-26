#!/usr/bin/env Rscript

# script to install Bioconductor R packages required for generation of reports by CRISPR-pooled-RSL pipeline
# based on https://github.com/casbap/ncRNA/blob/main/docker/rpkgs.R
# used in images based on rocker/tidyverse:4.2.1

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.15")

BiocManager::install(c('MAGeCKFlute','edgeR','fgsea','clusterProfiler','enrichplot','DOSE','biomaRt','ReactomePA','reactome.db','org.Hs.eg.db','org.Mm.eg.db'),version = '3.15')

print("Install Bioconductor packages, done!")



