#!/usr/bin/env bash

WORKDIR=/config/workdir/project1

cd "$WORKDIR"/data/alignments/

samtools merge -o merged.bam mother.bam father.bam son.bam

samtools index merged.bam
