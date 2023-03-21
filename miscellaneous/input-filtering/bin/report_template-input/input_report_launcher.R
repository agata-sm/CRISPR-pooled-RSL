#!/usr/bin/env Rscript

rm(list=ls())

library(knitr)
args <- commandArgs(TRUE)

if (length(args) < 4) stop("Not all args are set; required: refdata.path inputdata.path inputlib.name refdata.pref")

inputdata.path <- args[1]
inputlib.name <- args[2]
refdata.path <- args[3]
refdata.pref <- args[4]
refdata.defined <- args[5]




# for the code execution
#path to count table for reference data
#refdata.path="/Users/agata.smialowska/NBISproj/5351_CRISPR_pipeline/data/ctables/Calle_Heldin.UMIcounts.csv"

#root dir for filtered input libraries
#inputdata.path="/Users/agata.smialowska/NBISproj/5351_CRISPR_pipeline/gh/rackham_tsts/input/Brunello_v2022/results/input_filtered"


wrk.dir=file.path(paste("report",inputlib.name,sep="."))
dir.create(wrk.dir, recursive = TRUE)


rmarkdown::render('report-input-filt.Rmd', output_file = file.path(wrk.dir,paste("Input_library_report",inputlib.name, 'html', sep=".")))
