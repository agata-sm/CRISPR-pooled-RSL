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

* java (tested on OracleJDK_11.0.9)

* Nextflow >= 21 (tested on version 21.10.6)

* perl v5.18.4

* perl modules: `Getopt::Long`, `List::Util`, `File::Basename`, `File::Path`, `Sort::Rank`. They are parts of most standard perl distributions (with the exception of `Sort::Rank`, which may need to be installed from CPAN, see [instructions](https://www.uppmax.uu.se/support/faq/software-faq/installing-local-perl-packages/)).

* R >= 4.1.1

* R packages: please see the script `bin/report_template/pipeline_report.vN.R`

* pandoc >= 2.10.1

* MAGeCK


These dependencies are included in a singularity containers, please see section *Containers*

## Usage

The pipeline can be run in several ways, depending on resources available.


Tu run the main workflow for analysis based on **reads**

```
nextflow run /proj/sllstore2017103/nbis5351/CRISPR-pooled-RSL/crispr-pooled-rsl.nf 
```

To run the alternative workflow for analysis based on **RSL counts**

```
nextflow run /proj/sllstore2017103/nbis5351/CRISPR-pooled-RSL/crispr-pooled-rsl.nf -entry RSL
```

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

reads:

```
module load java/OracleJDK_11.0.9
nextflow run /proj/sllstore2017103/nbis5351/CRISPR-pooled-RSL/crispr-pooled-rsl.nf -profile cluster,singularity
```

RSL:

```
module load java/OracleJDK_11.0.9
nextflow run /proj/sllstore2017103/nbis5351/CRISPR-pooled-RSL/crispr-pooled-rsl.nf -entry RSL -profile cluster,singularity
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
module load java/OracleJDK_11.0.9
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

*Currently it is not recommended to run the workflow in this manner. To do this, you need to comment out the commands to load modules, as these are going to produce error on a local system*

(reads analysis only)

Nextflow and the dependencies need to be installed. The procedure is the same as for the interactive session on Rackham (all dependencies need to be installed).

```
nextflow run /path/to/CRISPR-pooled-RSL/crispr-pooled-rsl.nf
```


### Config files

There are two configuration files used by the pipeline. While they both are called `nextflow.config`, they serve different purposes.

* `nextflow.config` in pipeline directory - contains information related to running the pipeline on Rackham, profiles and paths to container library and to `CrisprCounter.jar`; this file should not be modified;

* `nextflow.config` in project working directory - contains information related to the project: paths to input library, input fastq files, metadata and analysis parameters. Please see directory `config-files` for examples and more detailed explanation.


### Metadata files

Two files describing the samples and their relationships are required by the pipeline. They should be named as indicated and follow the tab-delimited format. Column order needs to be preserved.

* `metadata.txt` -  information on file name for each sample and experimental groups (treatments, etc);

* `comparisons.txt` - information on comparisons to be analysed.

Please see directory `config-files` for examples and more detailed explanation.

## Misc tools

*Still to be added*

The following tools are included in this repository. While not part of the analysis pipeline, they prepare the input libraries for analyses.

* tools to filter input libraries sequenced in technical replicates to detect *bona fide* RSL-sgRNA combinations under miscellaneous/input-filtering;

* script to generate `gmt` file from library design file;

* script to process fastq files to remove the set constant part in the middle of the read.


## Containers

Four containers are used by the pipeline:

* `staphb/fastqc` obtained from Dockerhub;

* `mageck-perl` created for this pipeline; obtained from Dockerhub or built locally;

* `rpckg-crispr-rep` created for this pipeline; obtained from Dockerhub or built locally;

* `crisprcounter-perl` created for this pipeline; can only be built locally;

Custom containers `mageck-perl`  and `rpckg-crispr-rep` can be obtained from Dockerhub or built locally using Docker and transferred to Rackham. `crisprcounter-perl` can only be built locally and transferred to Rackham, as it contains `CrisprCounter.jar` which cannot be distributed with this pipeline. 

These images (with the exception of `crisprcounter-perl`) can be fetched from their remote repositories when the pipeline is executed for the first time. However, the process of image building may take some time, and in case of `rpckg-crispr-rep` also resources, it is recommended to fetch them before running the pipeline and save them in a "container library" - a directory accessible to each pipeline run. The default configuration is that the pipeline searches for container images in the "container library".

Dockerfiles and other necessary files for custom containers are included in subdirectories under *containers*.
Instructions on how to build the images are in the readme file under *containers*.


<!-- 
## To do

* contamination detection


 -->