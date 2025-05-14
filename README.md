# CRISPR-pooled-RSL
Pipeline for processing data from pooled CRISPR screens with RSL barcodes.

## About

This pipeline performs two types of analyses of pooled CRISPR screens:

(1) based on total read counts per guide, using MAGeCK workflow for Robust Rank Aggregation (RRA);

(2) based on RSL counts, using custom processing scripts and MAGeCK workflow for pathway analysis2 to select enriched / depleted guides.

The results of (1) and (2) are used in Gene Set Enrichment Analysis (GSEA) to produce functional annotation.


**Note March 2025**: The pipeline has been adjusted for running on Dardel (PDC). The former Rackham specifc config file is still available, for reference.

**Note May 2025**: The config files have been updated for running on Dardel (PDC). The former Rackham specifc config files are still available, for reference.


### Processing of RSL data

This part of the pipeline contains scripts for processing data from RSL barcoded screens.

These scripts form a workflow for filtering RSL-guide reads which likely are artefacts due to sequencing errors. The procedure relies on (i) filtering Input libraries (sequenced in technical replicates) to remove RSL-guide combinations present too few times or only in one replicate. This prefiltered count table is then used for (ii) filtering the sample count tables to retain only reads present in the filtered Input libraries. An additional filter may be applied after this Input-based filtering.



## Installation

To install the pipeline

```
git clone https://github.com/agata-sm/CRISPR-pooled-RSL.git
```

To update the pipeline, when standing in the pipeline directory

```
git pull
```

**to run (May 2025)**

```
cd /cfs/klemming/projects/supr/sllstore2017103/software/CRISPR-pooled-RSL
git pull
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


These dependencies are included in singularity containers, please see section *Containers*.

## Usage

The pipeline can be run in several ways, depending on resources available.


Two main pipeline modes are analysis based on reads (default) or analysis based on RSL tags specified by `-entry RSL`.

Tu run the main workflow for analysis based on **reads**

```
export pipelineDir="/path/to/CRISPR-pooled-RSL/"
nextflow run $pipelineDir/crispr-pooled-rsl.nf 
```

To run the alternative workflow for analysis based on **RSL counts**

```
nextflow run $pipelineDir/crispr-pooled-rsl.nf -entry RSL
```

The pipeline can be run on Rackham (recommended) or locally.

### using SLURM queue on Dardel

**This is the preferred way to run the pipeline.** All other ways have not been recently tested.

To run the pipeline in this mode, each process is submitted as a separate job to SLURM queue. This saves time (tasks are parallelised). The process executing the main pipeline code is run in the login node until the run is completed.

It is more practical to run the pipeline in the background. This protects the run from accidental session interruption - for example when you connect remotely to the server and the session disconnects.

You can use several programs to achieve this, in this example we use [screen](https://linux.die.net/man/1/screen), which is usually already installed in your Linux distribtion. 

All the commands are run on the login node.

To access `screen` on Dardel, you need to load some modules:

```
module load PDC/23.12
module load screen
```

You can then start `screen`. Depending on the settings of your local system, you may need to add arguments to the `screen` command. On MacOS:

```
screen -T xterm
```

This command opens a new shell, which later can be detached from the current session on the login node. The nextflow process, which controls the pipeline execution, will run in this new shell.

To start, we need to load necessary modules:


```
module load nextflow/24.04.2
module load singularity/4.1.1-cpeGNU-23.12
```

We can / should set up nextflow and singularity cache paths, to avoid overcrowding user home directory with temporary files:

```
export NXF_HOME="/cfs/klemming/projects/supr/sllstore2017103/nbis_dardel/nxf"
export NXF_SINGULARITY_CACHEDIR="/cfs/klemming/projects/supr/sllstore2017103/nbis_dardel/nxf"
```

The pipeline in the test runs have been run from another directory:

```
export pipelineDir="/cfs/klemming/projects/supr/sllstore2017103/nbis_dardel/CRISPR-pooled-RSL"
```

**to run after updating pipeline as indicated above (May 2025)**

**reads workflow**

```
export pipelineDir="/cfs/klemming/projects/supr/sllstore2017103/software/CRISPR-pooled-RSL"
nextflow run ${pipelineDir}/crispr-pooled-rsl.nf -profile cluster,singularity
```

**RSL workflow**

```
export pipelineDir="/cfs/klemming/projects/supr/sllstore2017103/software/CRISPR-pooled-RSL"
nextflow run ${pipelineDir}/crispr-pooled-rsl.nf -profile cluster,singularity -entry RSL
```

**OBS!** To detach from `screen` session, press `Ctrl-a` (together), then `d`.


#### Config files on Dardel

Project specific `nextflow.config`, which resides in the pipeline execution directory contains paths on Dardel.


`nextflow.config`:

```
params {
        librarydesign = "/cfs/klemming/projects/supr/sllstore2017103/nbis_dardel/tst1_Xiaonan_HEK/Brunello_Library_USE_THIS_ONLY.csv"
        libraryinputfilt = "/cfs/klemming/projects/supr/sllstore2017103/software/tests/input-lib-test/Brunello_x2022/results/input_filtered/Brunello_x2022.3/Brunello_x2022.filtered.csv"

        project = "naiss2025-22-65"

        projname = "Xiaonan_HEK_tst2"

        fastqdir = "/cfs/klemming/projects/supr/sllstore2017103/B.Schmierer_24_03_Hakan_Xiaonan_Maarten_Rong_Mohan/files/P30654/FASTQ/Xiaonan_FASTQ/"
        sampleinfo = "/cfs/klemming/projects/supr/sllstore2017103/nbis_dardel/tst1_Xiaonan_HEK/metadata_XZ.txt"
        comparisons = "/cfs/klemming/projects/supr/sllstore2017103/nbis_dardel/tst1_Xiaonan_HEK/comparisons_XZ.txt"
        scatters = "/cfs/klemming/projects/supr/sllstore2017103/nbis_dardel/tst1_Xiaonan_HEK/scatters.txt"


        organism = "hs"

        //////////////////////////////
        // RSL filtering params
        filtRowSums = "8"


        /////////////////////////////
        // names of the paramteres and available values are in commented headers

         //mageck count params
        // norm-method (total, median, control)
        mageckCountNorm = "control"
        // type of control list (sgRNA, gene)
        mageckCountCtrl = "gene"
        // file with control genes; if path not given (see below), features marked CON* in librarydesign file will be used
        //control_file = "/proj/software/tests/crispr-screen-test/gprc_ctrl.txt"
        control_file = "/cfs/klemming/projects/supr/sllstore2017103/nbis_dardel/tst1_Xiaonan_HEK/Control_Olfac_gene.txt"
        // no custom file, use CON* features from file librarydesign
        control_file = ""
}
```



### using SLURM queue on Rackham

**This used to be the preferred way to run the pipeline.**  Dardel is now the national HPC resource, and hence the preferred running environment.

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
module load bioinfo-tools
module load Nextflow/22.10.1

export NXF_HOME=/proj/sllstore2017103/software/nextflow_tmp
export APPTAINERENV_TMPDIR="/proj/sllstore2017103/software/containers"

export pipelineDir="/proj/sllstore2017103/software/CRISPR-pooled-RSL/"

nextflow run $pipelineDir/crispr-pooled-rsl.nf -profile cluster,singularity
```

RSL:

```
module load java/OracleJDK_11.0.9
module load bioinfo-tools
module load Nextflow/22.10.1

export NXF_HOME=/proj/sllstore2017103/software/nextflow_tmp
export APPTAINERENV_TMPDIR="/proj/sllstore2017103/software/containers"

export pipelineDir="/proj/sllstore2017103/software/CRISPR-pooled-RSL/"

nextflow run $pipelineDir/crispr-pooled-rsl.nf -entry RSL -profile cluster,singularity
```


To disconnect (detach): `Ctrl+a` `d` (you can do it when the nextflow command is running)

To reconnect: `screen -r`

To view the progress of the jobs in the queue: `jobinfo -u $USER`


### interactive session on Rackam


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

To run the pipeline several configuration and metadata files need to be present in the working directory. They are described in detail below, and the examples are given in `config-files`.

To run the pipeline simply type:

```
nextflow run $pipelineDir/crispr-pooled-rsl.nf
```

Where `$pipelineDir` is the directory where the pipeline is installed (see above).


To resume the run (if it was interrupted or you realised some arguments need to be changed):

```
nextflow run $pipelineDir/crispr-pooled-rsl.nf -resume
```


### locally 

*Currently it is not recommended to run the workflow in this manner. To do this, you need to make sure all dependencies are installed and use the default profile or have singularity installed and use profile singularity; you may need to edit nextflow.config to do this.*

(reads analysis only due to high memory requirements for RSL-sgRNA summarisation)

Nextflow and the dependencies need to be installed. The procedure is the same as for the interactive session on Rackham (all dependencies need to be installed).

```
nextflow run $pipelineDir/crispr-pooled-rsl.nf
```


### Config files

There are two configuration files used by the pipeline. While they both are called `nextflow.config`, they serve different purposes.

* `nextflow.config` in pipeline directory - contains information related to running the pipeline on Rackham, profiles and paths to container library and to `CrisprCounter.jar`; this file should not be modified;

* `nextflow.config` in project working directory - contains information related to the project: paths to input library, input fastq files, metadata and analysis parameters. Please see directory `config-files` for examples and more detailed explanation.


Project specific config file `myProj.config` can be used:

```
nextflow run $pipelineDir/crispr-pooled-rsl.nf -entry RSL -profile cluster,singularity -c myProj.config
```



### Metadata files

Two files describing the samples and their relationships are required by the pipeline. The paths to these files are given in the run specific config file. The files have to follow the tab-delimited format. Column order needs to be preserved.

* `metadata.txt` -  information on file name for each sample and experimental groups (treatments, etc);

* `comparisons.txt` - information on comparisons to be analysed.

Please see directory `config-files` for examples and more detailed explanation.


## Misc tools

The following tools are included in this repository. While not part of the analysis pipeline, they prepare the input libraries for analyses.

* tools to filter input libraries sequenced in technical replicates to detect *bona fide* RSL-sgRNA combinations under miscellaneous/input-filtering;

* script to process fastq files to remove the set constant part in the middle of the read.


## Containers

Four containers are used by the pipeline:

* `staphb/fastqc` obtained from Dockerhub;

* `mageck-perl` created for this pipeline; obtained from Dockerhub or built locally;

* `rpckg-crispr-rep` created for this pipeline; obtained from Dockerhub or built locally;

* `crisprcounter-perl` created for this pipeline; can only be built locally;

* `rpckg-input-report` created for Input library processing pipeline;

* `perl518_list_someutils` created for Input library processing pipeline;

Custom containers `mageck-perl`  and `rpckg-crispr-rep` can be obtained from Dockerhub or built locally using Docker and transferred to Rackham. `crisprcounter-perl` can only be built locally and transferred to Rackham, as it contains `CrisprCounter.jar` which cannot be distributed with this pipeline. 

These images (with the exception of `crisprcounter-perl`) can be fetched from their remote repositories when the pipeline is executed for the first time. However, the process of image building may take some time, and in case of `rpckg-crispr-rep` also resources, it is recommended to fetch them before running the pipeline and save them in a "container library" - a directory accessible to each pipeline run. The default configuration is that the pipeline searches for container images in the "container library".

Dockerfiles and other necessary files for custom containers are included in subdirectories under *containers*.
Instructions on how to build the images are in the readme file under *containers*.


<!-- 
## To do

* contamination detection


 -->
