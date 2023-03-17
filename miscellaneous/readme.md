# Miscallaneous tools distributed with CRISPR-pooled-RSL pipeline

The following tools are included in this repository. While not part of the analysis pipeline, they prepare the input libraries for analyses.

* tools to filter input libraries sequenced in technical replicates to detect *bona fide* RSL-sgRNA combinations under `miscellaneous/input-filtering`;

* script to process fastq files to remove the set constant part in the middle of the read `fastq_trim_mid_v2.pl`.


## Input filtering

This is a pipeline to filter **input libraries** sequenced in two technical replicates to assess the effect of different filtering cutoffs (range 1 to 10) on data. The output is a report with various descriptive statistics and plots.

The pipeline uses apptainer / singularty containers.

Usage:

```
module load java/OracleJDK_11.0.9
module load bioinfo-tools
module load Nextflow/22.10.1

export NXF_HOME="/proj/sllstore2017103/software/nextflow_tmp"
export APPTAINERENV_TMPDIR="/proj/sllstore2017103/software/containers"

nextflow run /path/to/Input-filter.nf -profile cluster, singularity
```

### Config file

As with the main pipeline, this pipeline uses a config file to set the run specific variables.

**nextflow.config** in pipeline run working directory:

```
params {

	// path to file proj.properties
	properties = "/proj/sllstore2017103/software/tests_2023/runs/input_fastq/test_run.properties"
  
  // computing project allocation
  project = "snic2022-22-1180"

	projname = "test_run"

  librarydesign = "/proj/sllstore2017103/software/tests/library_files/Brunello_Library_USE_THIS_ONLY.eol.csv"

  // whether to use reference data (recommended if available)
	usereference = "TRUE" // or "FALSE"

	// path to count table with evaluation data set
	refdatacnttable = "/proj/sllstore2017103/software/tests/ref_data/refdata.UMIcounts.csv"

	refdatapref = "refdata.UMIcounts" //

	// additional filtering by RowSums to remove low abundance RSL-sgRNA combinations
	filtRowSums = 1
  
}

```

Paths to fastq files of technical replicates if Input library are given in `test_run.properties`. Example of this file in given in the CRISPR facility repository or at `/proj/sllstore2017103/software/tests_2023/config_examples`


## fastq_trim_mid_v2.pl

Script to process fastq files: removes middle part of the fastq read (based on the position); works with fastq.gz.

Usage:

```
perl fastq_trim_mid_v2.pl --infile in.fastq --pos 24 --len 20
```

OBS! This script was tested on perl 5.26.x. It may not work on perl 5.30.x and above.

