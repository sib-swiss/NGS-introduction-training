#!/usr/bin/env bash

WORKDIR=/config/workdir/project1
cd "$WORKDIR"/data/fastq

fastqc *_trimmed_R?.fastq.gz
