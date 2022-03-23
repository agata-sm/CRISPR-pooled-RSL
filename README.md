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

To install the pipeline

```
git clone https://github.com/agata-sm/CRISPR-pooled-RSL.git
```

To use github from command line you need to setup [Personal Access Token PTA](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token). For the purpose of this task, the token is used in the same way as the password would be.




### Dependencies

* Nextflow >= 21 (tested on version 21.10.6)

* perl v5.18.4

* perl modules: `Getopt::Long`, `List::Util`, `File::Basename`, `File::Path`, `Sort::Rank`. They are parts of most standard perl distributions (with the exception of `Sort::Rank`, which may need to be installed from CPAN, see [instructions](https://www.uppmax.uu.se/support/faq/software-faq/installing-local-perl-packages/)).

* R >= 4.1.1

* R packages: please see the script `bin/report_template/pipeline_report.vN.R`

* pandoc >= 2.10.1

* MAGeCK


These dependencies will be included in a singularity container when it's ready.

## Usage

The pipeline can be run in several ways, depending on resources available.


***using SLURM queue on Rackham***

This is the preferred way to run the pipeline.

The file `nextflow.config` in *the directory where the pipeline is installed* contains profiles for pipeline execution. This is different than `nextflow.config` in the project directory where the code is run.

To run the pipeline in this mode, each process is submitted as a separate job to SLURM queue. This saves time (tasks are parallelised). The process executing the main pipeline code is run in the login node until the run is completed.

It is more practical to run the pipeline in the background. This protects the run from accidental session interruption - for example when you connect remotely to the server and the session disconnects.

You can use several programs to achieve this, in this example we use [screen](https://linux.die.net/man/1/screen), which is usually already installed in your Linux distribtion.

First, start the program by typing (on login node):

```
screen 
```

A new terminal appears. You can start a process in it, disconnect from it, then reconnect at any time.

To start a new screen press `Ctrl-a`, then `c`. To run the pipeline:

```
nextflow run /proj/sllstore2017103/nbis5351/CRISPR-pooled-RSL/crispr-pooled-rsl.nf -profile cluster
```

To disconnect (detach): `Ctrl+a` `d` (you can do it when the nextflow command is running)

To reconnect: `screen -r`

To view the progress of the jobs in the queue: `jobinfo -u $USER`


***interactive session on Rackam***


We start the interactive session:

```
interactive -A allocation_ID -t 2:00:00 -p core -n 4
```

Nextflow can easily be installed as described in [Nextflow documentation](https://www.nextflow.io/docs/latest/getstarted.html). 


We load the dependencies as modules:

```
module load perl/5.18.4
module load bioinfo-tools
module load MAGeCK/0.5.9.4
module load pandoc/2.17.1.1
module load R_packages/4.1.1
module load FastQC/0.11.9
```

To run the pipeline several configuration and metadata files need to present in the working directory. They are described in detail below, and the examples are given in `config-files`.

To run the pipeline simply type:

```
nextflow run /path/to/CRISPR-pooled-RSL/crispr-pooled-rsl.nf
```

Where `/path/to` is the directory where the pipeline is installed.


To resume the run (if it was interrupted or you realised some arguments need to be changed):

```
nextflow run /path/to/CRISPR-pooled-RSL/crispr-pooled-rsl.nf -resume
```


***locally on a laptop***

(reads analysis only)

Nextflow and the dependencies need to be installed. The procedure is the same as for the interactive session on Rackham.


```
nextflow run /path/to/CRISPR-pooled-RSL/crispr-pooled-rsl.nf
```


## Misc tools



## To do

* singularity container

* contamination detection


