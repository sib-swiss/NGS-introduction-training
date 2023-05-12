#!/usr/bin/env bash

WORKDIR=/config/workdir/projects/project1
cd "$WORKDIR"/data/fastq

fastqc *.fastq.gz
