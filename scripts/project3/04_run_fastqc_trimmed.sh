#!/usr/bin/env bash

WORKDIR=/config/workdir/projects/project3
cd "$WORKDIR"/results/trimmed/

fastqc *_trimmed_R?.fastq.gz
