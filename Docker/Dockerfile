FROM jupyter/minimal-notebook:9fe5186aba96

# change user to root to modify .bashrc for e.g. highlighting in the console
USER root
RUN cat /etc/skel/.bashrc >> /etc/bash.bashrc

# install bc for counting reads
RUN apt-get update && apt-get install -y bc

# change back to notebook user
USER $NB_UID

# remove auto-generated work directory
RUN rm -r /home/jovyan/work

# intall required software with conda
RUN /opt/conda/bin/conda install -y --quiet -c bioconda \
python=3.8 \
samtools=1.12 \
bwa=0.7.17 \
fastqc=0.11.8 \
sra-tools=2.11.0 \
cutadapt=3.4 \
bowtie2=2.4.2 \
hisat2=2.2.1 \
subread=2.0.1 \
entrez-direct=13.9
