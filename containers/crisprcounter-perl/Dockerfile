FROM ubuntu:5.26.2

LABEL crisprcounter-perl based on perl:5.26.2
LABEL use for running CrisprCounter.jar as part of the CRISPR-pooled-RSL pipeline

RUN cpanm install Getopt::Long
RUN cpanm install List::Util
RUN cpanm install File::Basename
RUN cpanm install Sort::Rank
RUN cpanm install File::Path


# Install OpenJDK-8
RUN echo 'deb http://ftp.debian.org/debian stretch-backports main' | tee /etc/apt/sources.list.d/stretch-backports.list
	# Install OpenJDK-8
	RUN apt-get update && \
	    apt-get install -y openjdk-8-jdk && \
	    apt-get install -y ant && \
	    apt-get clean;
	    
	# Fix certificate issues
	RUN apt-get update && \
	    apt-get install ca-certificates-java && \
	    apt-get clean && \
	    update-ca-certificates -f;

# Setup JAVA_HOME -- useful for docker commandline
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
RUN export JAVA_HOME

RUN mkdir ../opt/myjar
COPY CrisprCounter.jar ../opt/myjar

WORKDIR ./app
