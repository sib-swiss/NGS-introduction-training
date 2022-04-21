#!/usr/bin/env bash

cd ~/workdir
mkdir reads
cd reads
prefetch SRR519926
fastq-dump --split-files SRR519926
