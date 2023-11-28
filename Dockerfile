FROM ubuntu:rolling

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
    samtools \
    bedtools \
    bowtie2 \
    && apt-get clean

# creating virtual environment 
ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m virtualenv --python=/usr/bin/python3 $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Instal python dependencies
COPY requirements.txt .
RUN pip install -r requirements.txt