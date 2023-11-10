FROM ubuntu:rolling

RUN apt-get update \
    && apt-get install --no-install-recommends -y \
    curl \
    python3 \
    python3-pip \ 
    default-jre \
    && apt-get clean

# Nextflow
RUN curl -s https://get.nextflow.io | bash \
    && mv nextflow /usr/local/bin/

# Bioinfo tools 
RUN apt-get install --no-install-recommends -y \
    samtools \
    bedtools \
    bowtie2 \
    && apt-get clean