# Configuration and metadata files

## Config files

There are two configuration files used by the pipeline. While they both are called `nextflow.config`, they serve different purposes.

* `nextflow.config` in pipeline directory - contains information related to running the pipeline on Rackham, profiles and path to `CrisprCounter.jar`; this file should not be modified;

* `nextflow.config` in project working directory - contains information related to the project: paths to inpt library, input fastq files, metadata and analysis parameters. 

File below can be used to perform analysis on the dataset that:

* used Brunello library, the necessary files are in `/proj/AB1234/libraries/Brunello`

* allocation for computing resources is `snic2021-NN-XXX`

* project name is `MyProject.UMIcounts`

* fastq files are saved in `/proj/AB1234/tst_data/fastq/`

* sample info and comparisons are saved in `/proj/AB1234/tst_data/metadata`

* organism is human `hs`; the other choice is mouse `mm`


Analysis parameters that can be modified are:

* cutoff for RowSums when filtering off RSL-sgRNA artefacts `filtRowSums`


```
params {
        librarydesignRSL = "/proj/AB1234/libraries/Brunello/UMICount_Brunello_Library_USE_THIS_ONLY.tsv"
        librarydesign = "/proj/AB1234/libraries/Brunello/Brunello_Library_USE_THIS_ONLY.csv"
        libraryinputfilt = "/proj/AB1234/libraries/Brunello/Input_Brunello_mod.filtered.csv"
        libraryGMT = "/proj/AB1234/libraries/Brunello/Brunello.gmt"

        project = "snic2021-NN-XXX"

        projname = "MyProject.UMIcounts"

        fastqdir = "/proj/AB1234/tst_data/fastq/"
        sampleinfo = "/proj/AB1234/tst_data/metadata/metadata.txt"
        comparisons = "/proj/AB1234/tst_data/metadata/comparisons.txt"
        organism = "hs"
        
        filtRowSums = "0"
}

```



## Metadata files

Two files describing the samples and their relationships are required by the pipeline. Format requirements:

* Sßhould be named as indicated and follow the tab-delimited format;

* Column names and order needs to be preserved;

* No trailing empty lines are allowed (last line end-of-line character is ok);

* Sample names cannot start with a digit.


The following files describe dataset consising of four fastq files. The comparisons we are interested in are `high vs. low GFP`, within each replicate. Please make sire the sample names in `reference` and `treatment` columns in `comparisons.txt` are identical to `sample` in `metadata.txt`.


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


