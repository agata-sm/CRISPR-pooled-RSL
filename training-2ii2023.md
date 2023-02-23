# Training instructions 

3 Feb 2023
updated 23 Feb 2023

## Before we start

The current version of the pipeline is on Rackham at `/proj/sllstore2017103/software/CRISPR-pooled-RSL`.

This pipeline is written in `nextflow`, the workflow manager for high performance computing (HPC) servers, which is compatible with major software container environments, of which we use `Singularity / Apptainer` together with containers created by `Docker`.


The file whcih defines the workflow is `crispr-pooled-rsl.nf`.


Today we will work in directory `/proj/sllstore2017103/training_4ii2023`. Upon login to Uppmax please type:

```
cd /proj/sllstore2017103/training_4ii2023
```

Each person will work in their own directory:

```
cd your_dir
```



## Config files

There are several files present in `/proj/sllstore2017103/training_4ii2023/your_dir`.

* `metadata.txt` defines files, sample names and conditions;

* `comparisons.txt` defines comparisons;

* `gprc_ctrl.txt` is a list of control genes for normalisation. Please note that this is **for demonstration only**. Normalisation by control genes should be made with lists of at least 100 genes according to MAGeCK documentation;

* `nextflow.config` files formatted for the following normalisation menhods: **total** and **control** with *gene* (two variants: with ctrl genes from the library file and custom gene list gprc_ctrl.txt) and *sgRNAs* (only the library defined controls)



OBS! You will need to change the allocation in these files to your current one. 

```
# open the file in the text editor
nano nextflow.config_total

# go to line "project", remove the current allocation, paste the correct one

# press Ctrl-x
# press Y (you want to save changes)
# modify file name by removing the _total
# press enter
```

Now you saved the modified config file as `nextflow.config`.

This is the example:

```
cat nextflow.config_ctrl_sgRNA 
params {
	librarydesign = "/proj/sllstore2017103/software/tests/library_files/Brunello_Library_USE_THIS_ONLY.eol.csv"
	libraryinputfilt = "/proj/sllstore2017103/software/tests/input-lib-test/Brunello_x2022/results/input_filtered/Brunello_x2022.3/Brunello_x2022.filtered.csv"

	project = "snic2022-22-634"

	projname = "Yumeng_ctrl_sgRNA"

	fastqdir = "/proj/sllstore2017103/software/tests/ref_data/Yumeng2021/FASTQ_Yumeng"
	sampleinfo = "/proj/sllstore2017103/software/tests/ref_data/Yumeng2021/metadata.txt"
	comparisons = "/proj/sllstore2017103/software/tests/ref_data/Yumeng2021/comparisons.txt"
	organism = "hs"

	//////////////////////////////
	// RSL filtering params
	filtRowSums = "1"


	/////////////////////////////
	// names of the paramteres and available values are in commented headers

	//mageck count params
	// norm-method (total, median, control)
	mageckCountNorm = "control"
	// type of control list (sgRNA, gene)
	mageckCountCtrl = "sgRNA"
	//control_file = "/proj/sllstore2017103/software/tests/crispr-screen-test/gprc_ctrl.txt"
	control_file = ""

}
```




## Running the pipeline

First we run a program called `screen` which would allow us to disconnect the session when the pipeline is running.

```
screen
```

A new bash session starts (one which can be detached after the pipeline is started).

We load necessary modules and export path for Nextflow temporary files:

```
module load java/OracleJDK_11.0.9
module load bioinfo-tools
module load Nextflow/22.10.1

export NXF_HOME=/proj/sllstore2017103/software/nextflow_tmp
```


To start a new screen press Ctrl-a, then c. To run the pipeline:

reads:

```
module load java/OracleJDK_11.0.9
module load bioinfo-tools
module load Nextflow/22.10.1

export NXF_HOME=/proj/sllstore2017103/software/nextflow_tmp

nextflow run /proj/sllstore2017103/software/CRISPR-pooled-RSL/crispr-pooled-rsl.nf -profile cluster,singularity
```

RSL:

```
module load java/OracleJDK_11.0.9
module load bioinfo-tools
module load Nextflow/22.10.1

export NXF_HOME=/proj/sllstore2017103/software/nextflow_tmp

nextflow run /proj/sllstore2017103/software/CRISPR-pooled-RSL/crispr-pooled-rsl.nf -entry RSL -profile cluster,singularity
```


To disconnect (detach): `Ctrl+a d `(you can do it when the nextflow command is running)

To reconnect: `screen -r`

To view the progress of the jobs in the queue: `jobinfo -u $USER`



## To think about before the run

Common pitfals when running the pipeline:

1. Incorrect allocation in `nextflow.config` - has it expired? The error message says `Failed to submit process to grid scheduler for execution`

2. Library definition file e.g. `Brunello_Library_USE_THIS_ONLY.eol.csv` (see `nextflow.config` in the test directory for complete path) has incorrect line endings; the expected line endings are Linux style `\n`. If this file was open on any Windows / MacOS spreadsheet application, it is likely that the line endings have been converted to the Win or MacOS style, and this results in scripts not being able to read this file properly. The error messages may be less obvious, though. To run the pipeline the line endings need to be converted to Linux style.

A discussion on line endings and how to convert them between formats please check [here]https://kuantingchen04.github.io/line-endings/.


## Report recompilation


