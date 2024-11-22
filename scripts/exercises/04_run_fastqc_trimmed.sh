#!/usr/bin/env bash

# change directory to the trimmed directory
cd ~/project/results/trimmed

# run fastqc on all fastq files
fastqc trimmed*.fastq
