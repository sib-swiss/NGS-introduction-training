#!/usr/bin/env bash

WORKDIR=/config/workdir/project1

cd "$WORKDIR"/data/alignments/

samtools merge -o merged.bam mother.md.bam father.md.bam son.md.bam

samtools index merged.bam
