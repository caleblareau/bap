############################################################
# Dockerfile to build proatac
############################################################

FROM ubuntu:16.04

MAINTAINER Caleb Lareau; caleblareau@g.harvard.edu

# apt update and install global requirements
RUN apt-get clean all &&\
    apt-get update &&\
    apt-get upgrade -y && \
    apt-get install -y -qq  \
        autotools-dev   \
        automake        \
        cmake           \
        curl            \
        fuse            \
        git             \
        wget            \
        zip             \
        unzip           \
        libtbb-dev      \
        openjdk-8-jdk   \
        build-essential \
        r-base          \
        pkg-config      \
        python          \
        python-dev      \
        python-pip      \
        python3         \
        python3-dev     \
        python3-pip     \
        zlib1g-dev      \
        libncurses5-dev \
        libbz2-dev      \
        liblzma-dev     \
        libmagic-dev  &&\
    apt-get clean && \
    apt-get purge && \
    
    # Install python packages
    pip install python-magic && \
    pip install numpy && \
    pip install wheel && \
    pip install setuptools && \
	pip install macs2 && \
	pip3 install proatac && \
	
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /usr/share/doc/*

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



### Installing samtools/htslib/tabix/bgzip

ENV VERSIONH 1.4
ENV NAMEH htslib
ENV URLH "https://github.com/samtools/${NAMEH}/releases/download/${VERSIONH}/${NAMEH}-${VERSIONH}.tar.bz2"

ENV VERSION 1.4
ENV NAME "samtools"
ENV URL "https://github.com/samtools/${NAME}/releases/download/${VERSION}/${NAME}-${VERSION}.tar.bz2"

RUN wget -q $URLH && \
	bzip2 -d ${NAMEH}-${VERSIONH}.tar.bz2 && \
	tar -xf ${NAMEH}-${VERSIONH}.tar && \
	cd ${NAMEH}-${VERSIONH} && \
	./configure && \
	make -j 4 && \
	cd .. && \
	cp ./${NAMEH}-${VERSIONH}/tabix /usr/local/bin/ && \
	cp ./${NAMEH}-${VERSIONH}/bgzip /usr/local/bin/ && \
	cp ./${NAMEH}-${VERSIONH}/htsfile /usr/local/bin/ && \
	strip /usr/local/bin/tabix; true && \
	strip /usr/local/bin/bgzip; true && \
	strip /usr/local/bin/htsfile; true && \

	wget -q $URL && \
	bzip2 -d ${NAME}-${VERSION}.tar.bz2 && \
	tar -xf ${NAME}-${VERSION}.tar && \
	cd ${NAME}-${VERSION} && \
	./configure && \
	make -j 4 && \
	cd .. && \
	cp ./${NAME}-${VERSION}/${NAME} /usr/local/bin/ && \
	strip /usr/local/bin/${NAME}; true && \

	rm -rf ./${NAME}-${VERSION}/ && \
	rm -rf ./${NAME}-${VERSION}.tar && \
	rm -rf ./${NAMEH}-${VERSIONH}/ && \
	rm -rf ./${NAMEH}-${VERSIONH}.tar



### Installing bedtools

ENV VERSION 2.26.0
ENV NAME bedtools2
ENV URL "https://github.com/arq5x/bedtools2/releases/download/v${VERSION}/bedtools-${VERSION}.tar.gz"

RUN wget -q -O - $URL | tar -zxv && \
	cd ${NAME} && \
	make -j 4 && \
	cd .. && \
	cp ./${NAME}/bin/bedtools /usr/local/bin/ && \
	strip /usr/local/bin/*; true && \
	rm -rf ./${NAME}/
	
# Shouldn't need to update the working directory or copy anything. Maybe do a test like in Ariadne