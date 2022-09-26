# Containers for CRISPR-pooled-RSL pipeline

This directory contains Dockerfiles and auxiliary files for generating the following images:


* `mageck-perl` created for this pipeline; obtained from Dockerhub or built locally;

* `rpckg-crispr-rep` created for this pipeline; obtained from Dockerhub or built locally;

* `crisprcounter-perl` created for this pipeline; can only be built locally;

Custom containers `mageck-perl`  and `rpckg-crispr-rep` can be obtained from Dockerhub or built locally using Docker and transferred to Rackham. `crisprcounter-perl` can only be built locally and transferred to Rackham, as it contains `CrisprCounter.jar` which cannot be distributed with this pipeline. 

The images are Docker images, but as we cannot use Docker on Rackham, they are automatically converted to Singularity images by the pipeline.

## Retrieving the image from Dockerhub

The images can be downloaded to the "container library" prior to pipeline first run, to speed the task execution.

Let's consider that the container library is at `/crex/proj/MYPROJ/library/containers`.

1. set the directory for singularity cache files

```
export SINGULARITY_CACHEDIR="/crex/proj/MYPROJ/library/containers"
```

2. pull the images from Dockerhub:

```
singularity pull --name agatasm-mageck-perl.img docker://agatasm/mageck-perl

singularity pull --name agatasm-rpckg-crispr-rep.img docker://agatasm/rpckg-crispr-rep

singularity pull --name staphb-fastqc.img docker://staphb/fastqc
```


Please note that installing `rpckg-crispr-rep` requires more CPU usage than is normally accepted on the login node, therefore should be done on an interactive session instead.


## Building the image from Dockerfile locally

Image building is best inititiated in a directory which contains *only* the files required for the image (i.e. the content of subdirectories in this directory).

Installation of a `Docker` client is required, please see https://www.docker.com/.

If working on MacOS, initialisation of the Docker Desktop client (obtained from the link above) is required prior to following the instructions below. All steps require internet connection.

First, we build the image using the information in the Docker file:

```
docker build --no-cache -t crisprcounter-perl .
```

This image can be used to spin a container working in the local system. To move it to another location (such as Rackham):

1. we need to find out the tag of the image to package it:

```
docker images
```

the output (the actual IMAGE ID will be different):

```
REPOSITORY           TAG       IMAGE ID       CREATED         SIZE
crisprcounter-perl   latest    c233203f7b16   3 minutes ago   1.24GB
```

2. we save the image as tar archive:

```
docker save c233203f7b16 -o crisprcounter-perl.tar 
```

3. we copy the image to Rackham (using the singularity image library as the destination):

```
scp crisprcounter-perl.tar agata@rackham.uppmax.uu.se:/crex/proj/MYPROJ/library/containers
```

4. we build the singularity image on Rackham:

```
cd /crex/proj/MYPROJ/library/containers
singularity build crisprcounter-perl.sif docker-archive://crisprcounter-perl.tar
```

**For building `crisprcounter-perl` image `CrisprCounter.jar` needs to be copied to the directory with the Dockerfile.**
