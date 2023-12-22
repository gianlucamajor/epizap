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
    samtools=1.17-1 \
    bedtools=2.30.0+dfsg-2ubuntu1 \
    bowtie2=2.5.0-3 \
    muscle=1:5.1.0-1 \
    hmmer=3.3.2+dfsg-1 \
    ncbi-blast+=2.12.0+ds-3build1 \
    && apt-get clean

# creating virtual environment 
ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m virtualenv --python=/usr/bin/python3 $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Instal python dependencies
COPY requirements.txt .
RUN pip install -r requirements.txt