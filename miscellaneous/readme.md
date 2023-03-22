# Miscallaneous tools distributed with CRISPR-pooled-RSL pipeline

The following tools are included in this repository. While not part of the analysis pipeline, they prepare the input libraries for analyses.

* tools to filter input libraries sequenced in technical replicates to detect *bona fide* RSL-sgRNA combinations under `miscellaneous/input-filtering`;

* script to process fastq files to remove the set constant part in the middle of the read `fastq_trim_mid_v2.pl`.


## Input filtering

This is a pipeline to filter **input libraries** sequenced in two technical replicates to assess the effect of different filtering cutoffs (range 1 to 10) on data. The output is a report with various descriptive statistics and plots.

The pipeline uses Singularity / Apptainer containers.


Usage:

```
module load java/OracleJDK_11.0.9
module load bioinfo-tools
module load Nextflow/22.10.1

nextflow run /path/to/Input-filter.nf -profile cluster,singularity
```

### Containers

Three containers are used in this pipeline:

* `crisprcounter-perl.sif` described in the main pipeline;

* `agatasm/perl518_list_someutils` perl 5.18.4 with several packages;

* `agatasm/rpckg-input-report` R packages and other required for report generation;


Docker files for these containers are included with the main pipeline.


### Config files

Run specifc `nextflow.config` file should contain the following:

```
params {
	// allocation
	project = "snic2022-22-1180"

	// run name
	projname = "inputrun"

	// path to file proj.properties
	properties = "/path/to/inputrun.properties"

	// library design file
	librarydesign = "/path/to/Brunello_Library_USE_THIS_ONLY.eol.csv"

	// use reference data
	usereference = "TRUE" // or "FALSE"

	// path to count table with evaluation data set
	refdatacnttable = "/proj/sllstore2017103/software/tests/ref_data/Heldin2020/Calle_Heldin.UMIcounts.csv"

	refdatapref = "Calle_Heldin.UMIcounts"

	// additional filtering by RowSums to remove low abundance RSL-sgRNA combinations in reference data
	filtRowSums = 1
}
```

* `inputrun.properties` has format appropriate for `CrisprCounter.jar` and is described in detail elsewhere;

* `refdatacnttable` is a path to count table with reference data set;

The following prefixes should be identical: `projname` , `projname.properties` and the output file given in  `projname.properties`.


If desired, the run may be performed without using reference data (`usereference = "FALSE"`).



## fastq_trim_mid_v2.pl

Script to process fastq files: removes middle part of the fastq read (based on the position); works with fastq.gz.

Usage:

```
perl fastq_trim_mid_v2.pl --infile in.fastq --pos 24 --len 20
```

OBS! This script was tested on perl 5.26.x. It may not work on perl 5.30.x and above.

