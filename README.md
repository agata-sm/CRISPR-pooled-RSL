# CRISPR-pooled-RSL
Pipeline for processing data from pooled CRISPR screens with RSL barcodes.

## About

This pipeline performs two types of analyses of pooled CRISPR screens:

(1) based on all reads, using MAGeCK workflow for Robust Rank Aggregation (RRA);

(2) based on RSL counts, using custom processing scripts and MAGeCK workflow for GSEA to select enriched / depleted guides.

The results of (1) and (2) are used in Gene Set Enrichment Analysis (GSEA) to produce functional annotation.


### Processing of RSL data

This part of the pipeline contains scripts for processing data from RSL barcoded screens.

These scripts form a workflow for filtering RSL-guide reads which likely are artefacts due to sequencing errors. The procedure relies on (i) filtering Input libraries (sequenced in technical replicates) to remove RSL-guide combinations present too few times or only in one replicate. This prefiltered count table is then used for (ii) filtering the sample count tables to retain only reads present in the filtered Input libraries. An additional filter may be applied after this Input-based filtering.



## Installation


### Dependencies

* Nextflow > 0.20 (tested on version 21.10.6)

* perl v5.18.4

* perl modules: 

* R >= 4.1.1

* R packages: please see the script `bin/report_template/pipeline_report.vN.R`

* pandoc >= 2.10.1

* MAGeCK


These dependencies will be included in a singularity container when it's ready.

## Usage



## Misc tools



## To do

* singularity container

* fastqc

* contamination detection


