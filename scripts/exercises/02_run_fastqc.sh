#!/usr/bin/env bash

# change directory to the reads directory
cd ~/project/reads

# run fastqc on all fastq files
fastqc *.fastq
