# Set the base image to Ubuntu
FROM ubuntu:20.04
ENV DEBIAN_FRONTEND=noninteractive

# File Author / Maintainer
MAINTAINER Samantha Zarate

# System packages
RUN apt-get update -q
RUN apt-get upgrade -y -q

# Install dependencies for conda
RUN apt-get install -y \
    curl

# Install conda
RUN curl -LO http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh && bash Miniconda-latest-Linux-x86_64.sh -p /miniconda -b && rm Miniconda-latest-Linux-x86_64.sh
ENV PATH=/miniconda/bin:${PATH}
RUN conda update -y conda

# Install tools using conda
RUN conda config --add channels bioconda
RUN conda install -c bioconda \
    bwa \
    mosdepth \
    picard \
    samtools \
    openssl=1.0 \
    bcftools

RUN apt-get install -y \
    software-properties-common \
    unzip \
    wget

# Install Java 8
RUN wget -qO - https://adoptopenjdk.jfrog.io/adoptopenjdk/api/gpg/key/public | apt-key add -
RUN add-apt-repository --yes https://adoptopenjdk.jfrog.io/adoptopenjdk/deb/
RUN apt-get update 
RUN apt-get install -y \
    adoptopenjdk-8-hotspot

# Install GATK
RUN wget -q https://github.com/broadinstitute/gatk/releases/download/4.1.9.0/gatk-4.1.9.0.zip -O /gatk-4.1.9.0.zip
RUN unzip -q /gatk-4.1.9.0.zip -d /
ENV PATH=/gatk-4.1.9.0:${PATH}

# Install bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary
RUN mv bedtools.static.binary bin/bedtools
RUN chmod a+x bin/bedtools