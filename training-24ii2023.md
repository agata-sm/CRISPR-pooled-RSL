# Training 24 February 2023

1. Local report recompilation.
2. Running pipeline for processing Input libraries.

## Local report recompilation

Sometimes (often, always) we want to introduce project specific changes to the standard report. This can be done by modifying the report template and recompiling the report using this modified project specific template. 

Report reproducibility is accomplished by running the recompilation inside a *container*. Software container is a solution to bundle the OS and all required software and libraries and run it in isolation from the native local system. This means that by using the same container on Rackham and locally, the software environment is idential.

Locally (on our laptops) we will use Docker (the client needs to be installed from ). This is not possible on HPCs due to security issues, and Apptainer (formerly Singularity) is used instead; however, the image to run is equivalent as it is built using the same instructions.


### Docker image 

First we need to build the image. This is done only once, the same image can be used for all subsequent runs.

1. Start the Docker app. The app needs some time to configure itself. When ready to use, the command `docker -v` should result in printing the docker version.

2. We now can pull the image from Docker Hub:

```
docker pull agatasm/rpckg-crispr-rep
```

This process takes time, as the image is being created from scratch.


### Data

Report recompilation relies on data processed and collected during pipeline run. We need to copy it to our local system. It is important to preserve the structure of the directory, otherwise the script to produce report won't be able to find necessary files.

First, we locate the directory which contains the standard report from the pipeline run. Nextflow links all files used by a process to its own subdirectory located under `work`. Let's find the directory for `report.reads`. We will use the report from my most recent test run.

On Rackham:


on Rackham when in the pipeline run directory `/crex/proj/sllstore2017103/software/tests_2023/runs/Yumeng_1`

```
cd work
find . -name "report.reads"
```

The output of this comamnd is:

```
./be/81a4f78b2879ac36f00437eebf1d4d/report.reads
./ce/cd666308f533d2e0d8477a5d6d34e0/Yumeng_total_23ii2023/report.reads
./80/cf1115f1e532a718af19d0917a6016/Yumeng_total_23ii2023/report.reads
./1c/c22bb5f678c1c19e0b8aa2c45932a2/Yumeng_total_23ii2023/report.reads
./a9/f7cc52bbd0ce68aefea6d2967c60c9/report.reads
./a9/f7cc52bbd0ce68aefea6d2967c60c9/Yumeng_total_23ii2023/report.reads
```

This is because I run the pipeline twice (`report.reads` x2), and the RSL version three times (`Yumeng_total_23ii2023/report.reads` x3). Which is the correct one? Let's check time of file creation and use the newest one.

```
ll ./be/81a4f78b2879ac36f00437eebf1d4d/report.reads

ll ./a9/f7cc52bbd0ce68aefea6d2967c60c9/report.reads
```

The output:

```
agata@rackham2:work$ ll ./be/81a4f78b2879ac36f00437eebf1d4d/report.reads
total 24396
drwxrwsr-x 3 agata sllstore2017103     4096 Feb 23 22:43 CRISPR_Screen_report.Yumeng_total_23ii2023.reads_files
-rw-rw-r-- 1 agata sllstore2017103 24977353 Feb 23 23:03 CRISPR_Screen_report.Yumeng_total_23ii2023.reads.html
agata@rackham2:work$ ll ./a9/f7cc52bbd0ce68aefea6d2967c60c9/report.reads
total 24524
drwxrwsr-x 3 agata sllstore2017103     4096 Feb 24 01:11 CRISPR_Screen_report.Yumeng_total_23ii2023.reads_files
-rw-rw-r-- 1 agata sllstore2017103 25104956 Feb 24 01:33 CRISPR_Screen_report.Yumeng_total_23ii2023.reads.html
```

Let's work with

```
./a9/f7cc52bbd0ce68aefea6d2967c60c9/report.reads
```

We need to copy the entire directory to our local system:

```
scp -r <USER>@rackham.uppmax.uu.se:/proj/sllstore2017103/software/tests_2023/runs/Yumeng_1/work/a9/f7cc52bbd0ce68aefea6d2967c60c9/report.reads
```

**OBS! This is a large chunk of data!** For this test data set it's 2.5 GB.


To avoid confusing old report and formatted results, you can remove (or copy to a different location) the following directories, which will be created during the recompilation:

* ``report.reads`` -  the report itself

* ``results/reads/rra_annotation`` - formatted results and plots


```
rm -r report.RSL
rm -r results/RSL/rra_annotation
```

### Recompilation

We are going to work *locally* in the directory we have just copied.

```
cd f7cc52bbd0ce68aefea6d2967c60c9/report.reads
```

File ´crispr_pipeline_report_v0.4.Rmd` contains the main template for the report. The text
can be edited, following syntax guidelines for markdown (e.g. https://www.markdownguide.org/basicsyntax/).
Please do not change text in the code cunks and in other `Rmd`  files present in
this directory. When ready, the report may be recompiled.
The container is then started with a shell, which can be used to run commands.

To initialise Docker container:

```
docker run --rm -it -v $(pwd):/usr/src/report agatasm/rpckg-crispr-rep bash
```

The console prompt will be something like:

```
root@a82aacd301b1:/usr/src/report# 
```

You can now recompile the report.


The comand is:

```
Rscript report_launcher.R PROJNAME PROJNAME data-type org
```


In our case PROJNAME is `Yumeng_total_23ii2023`, data type is reads and organism is hs:

```
Rscript report_launcher.R Yumeng_total_23ii2023 Yumeng_total_23ii2023 reads hs
```


Time for a cup of your beverage of choice, as it takes a while to complete.



## Processing Input libraries


This pipeline was tested recently and the last run can be found at

```
/proj/sllstore2017103/software/tests/input-lib-test
```


The command to run this pipeline is:

```
nextflow run /proj/sllstore2017103/software/tests_2023/CRISPR-pooled-RSL/miscellaneous/input-filtering/Input-filter.nf -profile cluster
```


