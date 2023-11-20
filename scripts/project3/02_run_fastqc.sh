#!/usr/bin/env bash

PROJDIR=/config/project/projects/project3
cd "$PROJDIR"/data/fastq

fastqc *.fastq.gz
