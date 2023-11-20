#!/usr/bin/env bash

PROJDIR=/config/project/projects/project1
cd "$PROJDIR"/results/trimmed/

fastqc *_trimmed_R?.fastq.gz
