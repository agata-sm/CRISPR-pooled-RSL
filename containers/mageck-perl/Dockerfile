FROM rocker/tidyverse:4.2.1
LABEL mageck-perl based on rocker/tidyverse:4.2.1
LABEL use for running mageck 0.5.9.4 as part of the CRISPR-pooled-RSL pipeline

ENV SHELL=bash
ENV DEBIAN_FRONTEND noninteractive

WORKDIR .

RUN apt-get update && apt-get install perlbrew --no-install-recommends -y
RUN mkdir -p root/perl5/perlbrew/dists
RUN perlbrew install-patchperl
RUN perlbrew --force install perl-5.18.4

RUN wget -O - http://cpanmin.us | perl - --self-upgrade
RUN cpanm install Getopt::Long
RUN cpanm install List::Util
RUN cpanm install File::Basename
RUN cpanm install Sort::Rank
RUN cpanm install File::Path

RUN Rscript  -e 'install.packages("xtable",repos="http://cran.fhcrc.org");' -e 'install.packages("gplots",repos="http://cran.fhcrc.org");' -e 'install.packages("ggplot2",repos="http://cran.fhcrc.org");'


RUN apt-get update && apt-get install g++ make texinfo texlive-latex-base -y
RUN apt-get install wget texlive texlive-latex-extra --no-install-recommends -y

RUN apt-get install python3 python3-pip -y

RUN pip3 install numpy && pip3 install scipy

RUN wget http://downloads.sourceforge.net/project/mageck/0.5/mageck-0.5.9.4.tar.gz && tar xvzf mageck-0.5.9.4.tar.gz && cd mageck-0.5.9.4 && python3 setup.py install 

RUN apt-get autoclean -y && apt-get autoremove  -y &&  apt-get clean -y
