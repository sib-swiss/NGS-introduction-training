#!/usr/bin/env bash

WORKDIR=/config/workdir/projects/project3
cd "$WORKDIR"/data/fastq

fastqc *.fastq.gz
