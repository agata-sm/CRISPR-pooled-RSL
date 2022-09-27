# Miscallaneous tools distributed with CRISPR-pooled-RSL pipeline

The following tools are included in this repository. While not part of the analysis pipeline, they prepare the input libraries for analyses.

* tools to filter input libraries sequenced in technical replicates to detect *bona fide* RSL-sgRNA combinations under `miscellaneous/input-filtering`;

* script to process fastq files to remove the set constant part in the middle of the read `fastq_trim_mid_v2.pl`.


## Input filtering

This is a pipeline to filter **input libraries** sequenced in two technical replicates to assess the effect of different filtering cutoffs (range 1 to 10) on data. The output is a report with various descriptive statistics and plots.

The pipeline uses module system on Rackham.

Usage:

```
module load java/OracleJDK_11.0.9
nextflow run /path/to/Input-filter.nf -profile cluster
```


## fastq_trim_mid_v2.pl

Script to process fastq files: removes middle part of the fastq read (based on the position); works with fastq.gz.

Usage:

```
perl fastq_trim_mid_v2.pl --infile in.fastq --pos 24 --len 20
```



