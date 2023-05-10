#!/usr/bin/env bash

WORKDIR=/config/workdir/project1
cd "$WORKDIR"/data/fastq

fastqc *.fastq.gz
