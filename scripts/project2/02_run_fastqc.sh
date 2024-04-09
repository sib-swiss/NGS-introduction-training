#!/usr/bin/env bash

PROJDIR=/config/workspace/projects/project2
cd "$PROJDIR"/reads

fastqc *.fastq.gz
