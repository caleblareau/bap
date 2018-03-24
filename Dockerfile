############################################################
# Dockerfile to build bap-- IN PROGRESS
############################################################

FROM r-base:3.4.4

MAINTAINER Caleb Lareau; caleblareau@g.harvard.edu

RUN apt-get update && apt-get install -y python3 
RUN cp /usr/bin/python3 /usr/bin/python

### Installing necessary R packages
RUN Rscript -e 'install.packages("devtools", repos = "http://cran.us.r-project.org"); \
				devtools::install_github("caleblareau/BuenColors")'
	
### Installing bowtie2
ENV VERSION 2.3.0
ENV NAME bowtie2
ENV URL "https://github.com/BenLangmead/bowtie2/archive/v${VERSION}.tar.gz"

RUN wget -q -O - $URL | tar -zxv && \
    cd ${NAME}-${VERSION} && \
    make -j 4 && \
    cd .. && \
    cp ./${NAME}-${VERSION}/${NAME} /usr/local/bin/ && \
    cp ./${NAME}-${VERSION}/${NAME}-* /usr/local/bin/ && \
    strip /usr/local/bin/*; true && \
	rm -rf ./${NAME}-${VERSION}/

