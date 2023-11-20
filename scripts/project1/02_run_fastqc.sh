#!/usr/bin/env bash

PROJDIR=/config/project/projects/project1
cd "$PROJDIR"/data/fastq

fastqc *.fastq.gz
