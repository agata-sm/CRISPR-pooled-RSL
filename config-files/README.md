# Configuration and metadata files

## Config files

There are two configuration files used by the pipeline. While they both are called `nextflow.config`, they are located in different directories and serve different purposes.

* `nextflow.config` in pipeline directory - contains information related to running the pipeline on Rackham, profiles and paths to container directory and `CrisprCounter.jar`; usually this file should not be modified;

* `nextflow.config` in project working directory - contains information related to the project: paths to input library, input fastq files, metadata and analysis parameters. 

File below can be used to perform analysis on the dataset that:

* used Brunello library, the necessary files are in `/proj/AB1234/libraries/Brunello`; please make sure the file has **unix compatible line endings**;

* allocation for computing resources is `snic2021-NN-XXX`;

* project name is `MyProject.UMIcounts`; this will be appended to many output files and used throughout the report;

* fastq files are saved in `/proj/AB1234/tst_data/fastq/`;

* sample info and comparisons are saved in `/proj/AB1234/tst_data/metadata`;

* organism is human `hs`; the other choice is mouse `mm`;


Analysis parameters that can be modified are:

* cutoff for RowSums when filtering off RSL-sgRNA artefacts `filtRowSums` when performing signal filtration based on input library sequencing; this is applied after the input-based sequencing and is designed to remove RSL-sgRNA combinations present at very low abundance throughout the experiment;

* normalisation used by `MAGeCK` when read counting; one of `total`, `median`;

```
params {
        librarydesign = "/proj/AB1234/libraries/Brunello/Brunello_Library_USE_THIS_ONLY.eol.csv"
        libraryinputfilt = "/proj/AB1234/libraries/input_filtered/Brunello.filtered.csv"

        project = "snic2021-NN-XXX"

        projname = "MyProject.UMIcounts"

        fastqdir = "/proj/AB1234/tst_data/fastq/"
        sampleinfo = "/proj/AB1234/tst_data/metadata/metadata.txt"
        comparisons = "/proj/AB1234/tst_data/metadata/comparisons.txt"
        organism = "hs"
        
        //////////////////////////////
        // RSL filtering params
        filtRowSums = "0"


        /////////////////////////////
        // names of the paramteres and available values are in commented headers

        //mageck count params
        // norm-method (total, median)
        mageckCountNorm = "total"
}

```



## Metadata files

Two files describing the samples and their relationships are required by the pipeline. Format requirements:

* Should be named as indicated and follow the tab-delimited format;

* Column names and order needs to be preserved;

* No trailing empty lines are allowed (last line end-of-line character is ok);

* **Sample names cannot start with a digit**.


The following files describe dataset consising of four `fastq` files. The comparisons we are interested in are `high vs. low GFP`, within each replicate. Please make sure the sample names in `reference` and `treatment` columns in `comparisons.txt` are identical to `sample` in `metadata.txt`.
`name` in `comparisons.txt` refers to comparison name, and will be used in naming output and throughout the report.

* `metadata.txt` -  information on file name for each sample and experimental groups (treatments, etc);

```
file	sample	condition	group
LowGFP_REP1_R1_001.fastq.gz	LowGFP_REP1	LowGFP	rep1
LowGFP_REP2_R1_001.fastq.gz	LowGFP_REP2	LowGFP	rep2
HighGFP_REP1_R1_001.fastq.gz	HighGFP_REP1	HighGFP	rep1
HighGFP_REP2_R1_001.fastq.gz	HighGFP_REP2	HighGFP	rep2
```


* `comparisons.txt` - information on comparisons to be analysed.


```
name	reference	treatment
replicate-1	LowGFP_REP1	HighGFP_REP1
replicate-2	LowGFP_REP2	HighGFP_REP2
```
