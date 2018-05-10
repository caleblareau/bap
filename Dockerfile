#########################
# Dockerfile to build bap
#########################

From r-base:3.5.0

MAINTAINER Caleb Lareau; clareau@broadinstitute.org
ENV SHELL bash

### Installing necessary R packages
RUN Rscript -e 'install.packages(c("data.table", "dplyr", "ggplot2"), repos = "http://cran.us.r-project.org"); \
				source("https://bioconductor.org/biocLite.R");\
				biocLite("GenomicAlignments"); \
				biocLite("GenomicRanges"); \
				biocLite("Rsamtools"); \ 
				biocLite("BiocParallel")'
				
RUN apt-get update && apt-get install -y python3-pip
RUN cp /usr/bin/python3 /usr/bin/python

COPY . /bap
WORKDIR /bap

RUN pip3 install .

RUN bap bam -i data/test.small.bam -z -bt CB -ji 0.002 
