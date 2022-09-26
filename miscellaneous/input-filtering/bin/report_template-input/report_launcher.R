#!/usr/bin/env Rscript

rm(list=ls())

library(knitr)
args <- commandArgs(TRUE)

if (length(args) < 3) stop("Not all args are set; required: indir name type organism (hs or mm)")

proj.dir <- args[1]
proj.name.pref <- args[2]
data.type <- args[3]
organism <- args[4]

inputlib.name

refdata.name

refdata.pref

# for the code execution
#path to count table for reference data
refdata.path="/Users/agata.smialowska/NBISproj/5351_CRISPR_pipeline/data/ctables/Calle_Heldin.UMIcounts.csv"

#root dir for filtered input libraries
inputdata.path="/Users/agata.smialowska/NBISproj/5351_CRISPR_pipeline/gh/rackham_tsts/input/Brunello_v2022/results/input_filtered"



#wrk.dir=file.path(proj.dir,"results",data.type,"report")
wrk.dir=file.path(paste("report",inputlib.name,sep="."))
dir.create(wrk.dir, recursive = TRUE)



rmarkdown::render('report-input-filt.Rmd', output_file = file.path(wrk.dir,paste("Input_library_report",proj.name.pref, 'html', sep=".")))
