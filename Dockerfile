FROM ubuntu:24.04

RUN apt-get update \
    && apt-get install --no-install-recommends -y \
    curl \
    python3 \
    python3-pip \
    python3-virtualenv \ 
    default-jre \
    && apt-get clean

# Nextflow
RUN curl -s https://get.nextflow.io | bash \
    && mv nextflow /usr/local/bin/

# Bioinfo tools 
RUN apt-get install --no-install-recommends -y \
    samtools=1.19.2-1build2 \
    bedtools=2.31.1+dfsg-2 \
    bowtie2=2.5.2-1 \
    muscle=1:5.1.0-1 \
    hmmer=3.4+dfsg-2 \
    ncbi-blast+=2.12.0+ds-4build2 \
    && apt-get clean

# install mview dependence version 1.67
RUN apt-get install --no-install-recommends -y \
    libfile-copy-recursive-perl=0.45-4 \
    && apt-get clean

# install mview version 1.67
COPY thirdPartySoft/mview-1.67.tar.gz .
RUN tar xvzf mview-1.67.tar.gz
WORKDIR /mview-1.67
RUN perl install.pl /usr/bin/



# creating virtual environment 
ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m virtualenv --python=/usr/bin/python3 $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Install python dependencies
COPY requirements.txt .
RUN pip install -r requirements.txt