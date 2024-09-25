# Configuration and metadata files

## Config files

There are two configuration files used by the pipeline. While they both are called `nextflow.config`, they are located in different directories and serve different purposes.

* `nextflow.config` in pipeline directory - contains information related to running the pipeline on Rackham / Dardel, profiles and paths to container directory and `CrisprCounter.jar`; usually this file should not be modified;

* `nextflow.config` in project working directory - contains information related to the project: paths to input library, input fastq files, metadata and analysis parameters. 

File below can be used to perform analysis on the dataset that:

* used Brunello library, the necessary files are in `/proj/AB1234/libraries/Brunello`; please make sure the file has **unix compatible line endings**;

* allocation for computing resources is `snic2021-NN-XXX`;

* project name is `MyProject.UMIcounts`; this will be appended to many output files and used throughout the report;

* fastq files are saved in `/proj/AB1234/tst_data/fastq/`;

* sample info and comparisons are saved in `/proj/AB1234/tst_data/metadata`;

* organism is human `hs`; the other accepted choice is mouse `mm`;


Analysis parameters that can be modified are:

* cutoff for **RowSums** when filtering off RSL-sgRNA artefacts `filtRowSums` when performing signal filtration based on input library sequencing; this is applied after the input-based sequencing and is designed to remove RSL-sgRNA combinations present at very low abundance throughout the experiment; this setting has no effect on read-based analysis;

* normalisation used by `MAGeCK` when read counting; one of `total`, `median`, `control`;

Normalisation method `control` requires a file with control sgRNAs or genes. This file is generated from the library definition file `librarydesign` and contains features labeled as `CON*`; alternatively a custom file with feature list (one per line) can be provided. Please note that if no custom control file is used, it's entry should be nonetheless populated with `control_file = ""`, such as in the example below.



```
params {
        librarydesign = "/proj/AB1234/libraries/Brunello/Brunello_Library_USE_THIS_ONLY.eol.csv"
        libraryinputfilt = "/proj/AB1234/libraries/input_filtered/Brunello.filtered.csv"

        project = "snic2021-NN-XXX"

        projname = "MyProject.UMIcounts"

        fastqdir = "/proj/AB1234/tst_data/fastq/"
        sampleinfo = "/proj/AB1234/tst_data/metadata/metadata.txt"
        comparisons = "/proj/AB1234/tst_data/metadata/comparisons.txt"
        scatters = "/proj/AB1234/tst_data/metadata/scatters.txt"

        organism = "hs"
        
        //////////////////////////////
        // RSL filtering params
        filtRowSums = "0"


        /////////////////////////////
        // names of the paramteres and available values are in commented headers

         //mageck count params
        // norm-method (total, median, control)
        mageckCountNorm = "control"
        // type of control list (sgRNA, gene)
        mageckCountCtrl = "gene"
        // file with control genes; if path not given (see below), features marked CON* in librarydesign file will be used
        //control_file = "/proj/software/tests/crispr-screen-test/gprc_ctrl.txt"
        // no custom file, use CON* features from file librarydesign
        control_file = ""

}

```



## Metadata files

Two files describing the samples and their relationships are required by the pipeline. Format requirements:

* Should follow the tab-delimited format;

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

**Important** Sample name (column `sample`) requirements: 

        * alphanumeric characters (letters, digits and `_`) are allowed

        * must not start with a digit

        * must not contain hyphen `-`


* `comparisons.txt` - information on comparisons to be analysed.


```
name	reference	treatment
replicate-1	LowGFP_REP1	HighGFP_REP1
replicate-2	LowGFP_REP2	HighGFP_REP2
```

**name** in this file refers to comparison name (and will be referenced in the report).


## Output parameters

It is possible to output a document containing a collection of **scatter plots of log2FC values** in a series of *comparison pairs* of choice. This is set by parameter `scatters` by providing a path to file `scatters.txt` which lists *comparisons* to be plotted on scatter plots (by their comparison **name** used in file `comparisons.txt`). Arbitrary number of pairs for scatter plots can be defined, just be mindufl of file size and its memory footprint when displaying, as these plots are interactive.

If no scatters are to be plotted, please use `scatters = "none"` as the value of the `scatters` parameter (the parameter should be set).


