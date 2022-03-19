#!/usr/bin/env Rscript

rm(list=ls())

library(knitr)
args <- commandArgs(TRUE)

if (length(args) < 3) stop("Not all args are set; required: indir name type organism (hs or mm)")

proj.dir <- args[1]
proj.name.pref <- args[2]
data.type <- args[3]
organism <- args[4]

#wrk.dir=file.path(proj.dir,"results",data.type,"report")
wrk.dir=file.path(paste("report",data.type,sep="."))
dir.create(wrk.dir, recursive = TRUE)



rmarkdown::render('crispr_pipeline_report_v0.4.Rmd', output_file = file.path(wrk.dir,paste("CRISPR_Screen_report",proj.name.pref,data.type, 'html', sep=".")))
