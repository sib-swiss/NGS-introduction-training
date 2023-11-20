#!/usr/bin/env bash

cd ~/project
mkdir reads
cd reads
prefetch SRR519926
fastq-dump --split-files SRR519926
