#!/usr/bin/env bash

# change directory to project
cd ~/project
# create a directory for reads
mkdir reads
# change directory to reads
cd reads
# download reads from SRA
prefetch SRR519926
fastq-dump --split-files SRR519926
