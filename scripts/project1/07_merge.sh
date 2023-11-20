#!/usr/bin/env bash

PROJDIR=/config/project/projects/project1

cd "$PROJDIR"/alignments/

samtools merge -o merged.bam mother.md.bam father.md.bam son.md.bam

samtools index merged.bam
